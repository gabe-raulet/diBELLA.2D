// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "CC.h"

// #define MATRIXPOWER
#define APSP
#define MAXPATHLEN 5000

#ifdef APSP // All-Pairs Shortest Path (not great for big graphs, I wanna use Seidel's algorithm but modify it for sparse matrices)

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>

using namespace boost::numeric::ublas;
#define INF 10000

typedef struct
{
    int len;
    int dir; // I can use short (not ushort, -1 is undefined)
    std::string seq;
} Overlap;

void toBinary(ushort n, int* arr) 
{ 
    int nbit = 2;
    for(int i = 0; i < nbit; i++)
    { 
        arr[i] = n % 2; 
        n = n / 2; 
    }
}

// If false, dir (ij.dir) isn't modified
bool DirectionalityCheck(int dir1, int dir2, int& dir)
{
    int rbit, lbit;
    int start, end;

    int mybin1[2] = {0, 0}; // 0 1
    int mybin2[2] = {0, 0}; // 0 0

    if(dir1 != 0) toBinary(dir1, mybin1);
    if(dir2 != 0) toBinary(dir2, mybin2);

    rbit = mybin1[0]; // 1 
    lbit = mybin2[1]; // 0

    if(rbit != lbit)
    {
        start = mybin1[1]; // 11
        end   = mybin2[0]; // 00

        if(start == 0)
        {
            if(end == 0) dir = 0;
            else dir = 1;
        }
        else
        {
            if(end == 0) dir = 2;
            else dir = 3;      
        }

        return true;
    }
    else return false;
}

// If one of the two is INF, the condition ij > ik + kj is not gonna be satified
bool isDirOk(int dir1, int dir2, int& dir)
{
    if(DirectionalityCheck(dir1, dir2, dir)) return true;
    else return false;
}

void PrintAPSP(mapped_matrix<Overlap>& D, int nvertex)
{
    for (unsigned i = 0; i < D.size1(); ++i)
        for (unsigned j = 0; j < D.size2(); ++j)
        {
            Overlap entry = D(i, j);

            if(entry.len == INF) printf("%7s\t", "INF");
            else printf("%s\t", entry.seq.c_str());
            // else printf("%7d\t", entry.len);

            if(j == nvertex - 1) printf("\n");
        }
    printf("\n");
}

void APSP(mapped_matrix<Overlap>& S, int nvertex, int nnz)
{
    mapped_matrix<Overlap> D(S); // the distances matrix (copy constructor)
    unsigned i, j, k;
 
    for (k = 0; k < nvertex; k++)
        for (i = 0; i < nvertex; i++)
            for (j = 0; j < nvertex; j++)
            {
                Overlap ij = D(i, j);
                Overlap kj = D(k, j);
                Overlap ik = D(i, k);

                // The direction check is very simple right now, the bidirectional one hasn't been integrated/tested yet
                if((ij.len > (ik.len + kj.len)) & isDirOk(ik.dir, kj.dir, ij.dir))
                {
                    ij.len = ik.len + kj.len;
                    ij.seq = ik.seq + kj.seq;
                }

                D(i, j) = ij;
            }
    PrintAPSP(D, nvertex);
}

#endif

/*! Namespace declarations */
using namespace combblas;
// typedef ContigSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ContigSRing_t;

std::vector<std::string> 
CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu, 
    PSpMat<dibella::CommonKmers>::DCCols* spSeq, std::shared_ptr<DistributedFastaData> dfd, int64_t nreads)
{    

    float balance = S.LoadImbalance();
    int64_t nnz   = S.getnnz();

    std::ostringstream outs;
    outs.str("");
    outs.clear();
    outs << "CreateContig::LoadBalance: " << balance << endl;
    outs << "CreateContig::nonzeros: "    << nnz     << endl;
    SpParHelper::Print(outs.str());

    int64_t nCC = 0;
    FullyDistVec<int64_t, int64_t> myLabelCC = CC(S, nCC);

	std::string myccoutput = myoutput + ".cc";
	myLabelCC.ParallelWrite(myccoutput, 1);

    uint64_t nContig = myLabelCC.Reduce(maximum<int64_t>(), (int64_t)0);
    nContig++; // because of zero based indexing for cluster

    std::stringstream ncc;
    ncc << "nContig: " << nContig << endl;
    SpParHelper::Print(ncc.str());

    First4Clust(myLabelCC);
    HistCC(myLabelCC, nCC);

// GGGG: This must be executed sequentially for now (prototyping on E. coli CCS single node)
#ifdef APSP

    // 1) From CombBLAS to Boost matrix format

    int nvertex = nreads;
    // int nnz = S.getnnz();

    // The string matrix boostS
    mapped_matrix<Overlap> boostS(nvertex, nvertex, nnz); // This is not compressed (I need to modify Seidel's to get a sparse version)

    // Local sequences (entire dataset since it's sequential for now)
	uint64_t z = 0;
	auto dcsc = spSeq->GetDCSC();

    // I wanna have a vector of tuples <row idx, col idx, direction, len, suffix> to fill my boostS matrix
    VecTupType mattuples(nnz); // nnz in the vector

	for (uint64_t i = 0; i < dcsc->nzc; ++i)
	{
		for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
		{
			std::get<0>(mattuples[z]) = dcsc->ir[j]; // row idx
			std::get<1>(mattuples[z]) = dcsc->jc[i]; // col idx

            // This is the nonzero, meaning CommonKmers
            // I wanna extract the direction from here and the compute the string which is not currently stored there because CombBLAS doesn't like string)
            // To get the string I need direction, start/end position
			dibella::CommonKmers* cks = &(dcsc->numx[j]);
            int dir = cks->overhang & 3;
            std::get<2>(mattuples[z]) = dir;

            // Get row/col sequences paying attention to the consistency of V/H
            seqan::Dna5String seqH = *(dfd->row_seq(dcsc->ir[j])); // extract row sequence
            seqan::Dna5String seqV = *(dfd->row_seq(dcsc->jc[i])); // extract col sequence

            seqan::Dna5String overlap_suffix;

            // These should have already been updated according to overhang/overhangT during TR
            ushort begpV = cks->first.first;  // Updated post-alignment, need to update in TR semiring when transposing
			ushort begpH = cks->first.second; // Check correctness of H/V

            ushort endpV = cks->first.second;  // Updated post-alignment, need to update in TR semiring when transposing
			ushort endpH = cks->second.second; // Check correctness of H/V

            uint rlenV = length(seqV);
            uint rlenH = length(seqH);

            // Get suffix substring
			if(dir == 1 || dir == 2) // !reverse complement
			{
				if(begpH > begpV)
				{
					assert(dir == 1);

					// suffix = rlenV - endpV;
                	// suffix = seqV.substr(endpV, length(seqV)); 
                    overlap_suffix = seqan::suffix(seqV, endpV); // seqH entering in seqV so we extract the ending part of seqV (fwd strand)
				}	
				else
				{
                    assert(dir == 2);

					// suffix  = rlenH - endpH;
					// suffix = seqH.substr(endpH, length(seqH)); 
                    overlap_suffix = seqan::suffix(seqH, endpH); // seqH entering in seqV so we extract the ending part of seqV (fwd strand)

				} 
			}
			else
			{
				if((begpV > 0) & (begpH > 0) & (rlenV-endpV == 0) & (rlenH-endpH == 0))
				{
					assert(dir == 0);

					// suffix = begpV;
					// suffix = seqV.substr(0, begpV); 
                    // >----> <----<, I want to get the reverse complement of this substring because seqV is on the reverse strand in this case
                    overlap_suffix = seqan::prefix(seqV, begpV+1); // The 2nd par is excluding end position
                    
                    // RevComplement(overlap_suffix); // Declare this function!
				}
				else
				{
					assert(dir == 3);

					// suffix  = rlenV - endpV;	
                    // suffix = seqV.substr(endpV,  length(seqV)); 
                    // <----< >---->, I don't want get the reverse complement of this substring because seqV is on the fwd strand in this case
                    overlap_suffix = seqan::suffix(seqV, endpV);
				}
			}

            std::get<3>(mattuples[z]) = overlap_suffix; // Everything should technically be consistent, if correct
			++z;
		}
	}

	assert(z == nnz);
	std::cout << "Local nnz count: " << nnz << std::endl;

    // for (uint64_t i = 0; i < dcsc->nz; ++i)
	// {
    //     int64_t lrid = dcsc->ir[i]; // local row idx
    //     seqan::Dna5String rseq = *(dfd->row_seq(lrid)); // extract row sequence
    // }

    // The entire matrix boostS is initialized to INF with 0 on the diagonal
    for (unsigned i = 0; i < boostS.size1(); ++i)
        for (unsigned j = 0; j < boostS.size2(); ++j)
        {
            Overlap entry;
            entry.seq = "";

            entry.dir = -1; // if INF (zeros) direction must be "indefined", 0 is a direction in our bidirected graph

            if(i != j) entry.len = INF;
            else entry.len = 0; // this is important dude (zeros on the diagonal must be actual zeros) ---Question how do I modify this for sparsity?

            boostS(i, j) = entry;
        }
    
    // PrintAPSP(boostS, nvertex);

    // 2) I need to add substring to the nonzero
    // The matrix boostS stores the nonzeros plus the substring suffix

    // PrintAPSP(boostS, nvertex);
      
    // 3) APSP (inefficient for big matrix but let's see if contig makes sense first)
    // PrintAPSP(S, nvertex);
    // APSP(S, nvertex, nnz);

#endif

#ifdef MATRIXPOWER

    // PSpMat<dibella::CommonKmers>::MPI_DCCols T = S; // (T) traversal matrix
    // PSpMat<dibella::CommonKmers>::MPI_DCCols ContigM(S.getcommgrid()); // Contig matrix
    // dibella::CommonKmers defaultBVal; 

    // Read vector is gonna multiply the matrix and create contig from there
    // FullyDistVec<int64_t, std::array<char, MAXCONTIGLEN>> ReadVector(S.getcommgrid());
    // char* nt; // NT initial value

    CustomVectorEntry nt;
    /*
     * A vector of read ids that is my path [0, 1, 3, 67] (paths)
     * A vector of offset [10, 40, 50] (offsets) and offsets.size = paths.size -1 
     * Offset 10 tells me that i have to cut the last 10 bases of read1 and concatenate them to read0, then 40 from read3 and concatenate them to read1 etc. 
    */
    FullyDistVec<int64_t, char*> ReadVector(S.getcommgrid(), nreads, nt);
	auto dcsc = spSeq->GetDCSC();

    // GGGG: fill the read vector with sequences
    // IT * ir; //!< row indices, size nz
	for (uint64_t i = 0; i < dcsc->nz; ++i)
	{
		int64_t lrid = dcsc->ir[i]; // local row idx
        seqan::Dna5String rseq = *(dfd->row_seq(lrid));

        // std::array<char, MAXCONTIGLEN> crseq;
        std::string cpprseq;
        std::copy(begin(rseq), end(rseq), begin(cpprseq)); // @GGGG: this doesnt work

        char* crseq = new char[MAXCONTIGLEN];
        crseq = &cpprseq[0]; // C++14

        std::cout << rseq << std::endl;
    
        ReadVector.SetElement(lrid, crseq); // SetElement only work locally (owner)
	}

    // ReadVector.DebugPrint();
    // FullyDistVec<int64_t, std::array<char, MAXSEQLEN>> ContigVector(S.getcommgrid());
    // FullyDistVec<int64_t, char*> ContigVector(S.getcommgrid());

    // do
    // { 
    //     // ContigSR concatenates entries
    //     ReadVector = ContigVector;
    //     ContigVector = SpMV<ContigSR>(S, ReadVector);
                
    // } while (ReadVector != ContigVector); // Once the two vec are identical we're done

    // // GGGG: we know how long are the contig(s) from the CC so we could just extract those or can I use a FullySpDist vector and get only on contig?

    // ContigM.ParallelWriteMM("contig.miracle.mm", true, dibella::CkOutputMMHandler());

#endif

    std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    return myContigSet;
}

// FullyDistVec<int64_t, dibella::CommonKmers> ReduceV(T.getcommgrid());
// FullyDistVec<int64_t, dibella::CommonKmers> ContigV(T.getcommgrid());

// ContigV =  T.Reduce(Row, ReduceMSR_t(), NullValue);
// ReduceV = nT.Reduce(Row, ReduceMSR_t(), NullValue);

// useExtendedBinOp doesn't seem to be used anywhere, only passed as argument?
// ContigV.EWiseApply(ReduceV, GreaterSR_t(), IsNotEndContigSR_t(), false)