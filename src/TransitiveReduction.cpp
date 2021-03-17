#include <cmath>
#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "../include/cxxopts.hpp"
#include "../include/TransitiveReductionSR.hpp"
#include "../include/Utils.hpp"

#include <map>
#include <fstream>

/*! Namespace declarations */
using namespace combblas;

/*! Type definitions */
typedef MinPlusBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> MinPlusSR_t;
typedef ReduceMBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ReduceMSR_t;
typedef Bind2ndBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> Bind2ndSR_t;

void TransitiveReduction(PSpMat<dibella::CommonKmers>::MPI_DCCols& B)
{
    PSpMat<dibella::CommonKmers>::MPI_DCCols BT = B;
    BT.Transpose();
    if(!(BT == B))
    {
        B += BT;
    }
    B.PrintInfo();

    uint nnz, prev;
    double timeA2 = 0, timeC = 0, timeI = 0, timeA = 0;

    /* Gonna iterate on B until there are no more transitive edges to remove */
    do
    {
      prev = B.getnnz();

      /* Find two-hops neighbors
       * C = B^2
       */
      double start = MPI_Wtime();
      PSpMat<dibella::CommonKmers>::MPI_DCCols F = B;
      PSpMat<dibella::CommonKmers>::MPI_DCCols C = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, F);
      timeA2 += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix C = B^2: ");
      C.PrintInfo();
  #endif
    
      start = MPI_Wtime();
      FullyDistVec<int64_t, dibella::CommonKmers> vA(B.getcommgrid());

      dibella::CommonKmers id; 
      vA = B.Reduce(Row, ReduceMSR_t(), id);
      vA.Apply(PlusFBiSRing<dibella::CommonKmers, dibella::CommonKmers>());

      F.DimApply(Row, vA, Bind2ndSR_t());
      timeC += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix F = B + FUZZ: ");
      F.PrintInfo();
  #endif

      /* Find transitive edges that can be removed
       * I = F >= C 
       */
      start = MPI_Wtime();
      bool isLogicalNot = false;
      PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(F, C, GreaterBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), isLogicalNot, id);

      /* Prune potential zero-valued nonzeros */
      I.Prune(ZeroUnaryOp<bool>(), true);
      timeI += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix I = F >= B: ");
      I.PrintInfo();
  #endif

      /* Remove transitive edges
       * B = B .* not(I)
       */ 
      start = MPI_Wtime();
      isLogicalNot = true;
      B = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, I, EWiseMulOp<dibella::CommonKmers, bool>(), isLogicalNot, true);

      /* Prune zero-valued overhang */
      B.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
      timeA += MPI_Wtime() - start;

  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix B = B .* not(I): ");
      B.PrintInfo();
  #endif
      nnz = B.getnnz();
       
    } while (nnz != prev);

    tu.print_str("Matrix B, i.e AAt after transitive reduction: ");
    B.PrintInfo();

 #ifdef DIBELLA_DEBUG
    double maxtimeA2, maxtimeC, maxtimeI, maxtimeA;
    
    MPI_Reduce(&timeA2, &maxtimeA2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeC,  &maxtimeC,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeI,  &maxtimeI,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeA,  &maxtimeA,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myrank == 0)
    {
      std::cout << "TransitiveReduction:TimeA2 = " << maxtimeA2 << std::endl;
      std::cout << "TransitiveReduction:TimeC  = " <<  maxtimeC << std::endl;
      std::cout << "TransitiveReduction:TimeI  = " <<  maxtimeI << std::endl;
      std::cout << "TransitiveReduction:TimeA  = " <<  maxtimeA << std::endl;
    }
 #endif
}