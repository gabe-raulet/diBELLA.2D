// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "ParallelOps.hpp"
#include "CC.h"

// #define MATRIXPOWER
#define MAXPATHLEN 5000

/*! Namespace declarations */
using namespace combblas;

void ReadProcessorAssignment(FullyDistVec<int64_t, int64_t> *assignment, const FullyDistVec<int64_t, int64_t>& vCC, FullyDistVec<int64_t, int64_t>& ccSizes, int maxperproc)
{
    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    /* getting the ids of all contigs with at least 2 reads */
    FullyDistVec<int64_t, int64_t> contigIds = ccSizes.FindInds(bind2nd(std::greater<int64_t>(), 1));

    /* sort ccSizes */
    /* starting from procecssor 0, count off reads from the largest
     * contig until you are out of reads or you hit the maxperproc limit.
     * Continue onto the next largest contig when you run out of reads, and
     * continue onto the next process when you hit the maxperproc limit. Do this
     * until you run out of contigs. Shouldn't run out of processors because
     * maxperproc should be set to be larger than the number of reads divided by
     * the number of processors */
}


void CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu, PSpMat<dibella::CommonKmers>::DCCols* spSeq, std::shared_ptr<DistributedFastaData> dfd, int64_t nreads)
{

    float balance = S.LoadImbalance();
    int64_t nnz = S.getnnz();

    std::ostringstream outs;
    outs.str("");
    outs.clear();
    outs << "CreateContig::LoadBalance: " << balance << endl;
    outs << "CreateContig::nonzeros: "    << nnz     << endl;
    SpParHelper::Print(outs.str());

    /* std::shared_ptr<ParallelOps> parops = ParallelOps::init(NULL, NULL); */
    /* int myrank = parops->world_proc_rank; */
    /* int nprocs = parops->world_procs_count; */

    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    std::shared_ptr<CommGrid> fullWorld = parops->grid;


    /* boolean matrix for degree calculations */
    PSpMat<bool>::MPI_DCCols D1(S.getcommgrid());

    /* cast overlap matrix nonzeros to 1, for degree calculation */
    D1 = S;

    /* @degs1 is vector of degrees before branch adjacent edges are removed */
    FullyDistVec<int64_t, int64_t> degs1(D1.getcommgrid());

    D1.Reduce(degs1, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    FullyDistVec<int64_t, int64_t> branches = degs1.FindInds(bind2nd(std::greater<int64_t>(), 2));

    /* delete branch adjacenct edges */
    S.PruneFull(branches, branches);

    /* boolean matrix for degree calculations, AFTER branch edge deletion */
    PSpMat<bool>::MPI_DCCols D2(S.getcommgrid());

    D2 = S;

    /* @degs2 is vector of degrees after branch adjacent edges are removed */
    FullyDistVec<int64_t, int64_t> degs2(D2.getcommgrid());

    D2.Reduce(degs2, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    /* @roots vector stores indices of each end of each contig, so should have even
     * number of nonzeros */
    FullyDistVec<int64_t, int64_t> roots = degs2.Find(bind2nd(std::equal_to<int64_t>(), 1));

    /* GRGR TODO: assertion with roots length, figure out how to find nnz of roots */

    //PSpMat<bool>::MPI_DCCols CCMat(S.getcommgrid());
    //CCMat = S;

    /* Calculate connected components */
    int64_t nCC = 0;
    FullyDistVec<int64_t, int64_t> vCC = CC(S, nCC);

    /* @LocalCCSizes will store the local counts for connected components, which will be used
     * later to instantiate a distributed vector with the global counts. These vectors are
     * the same size in each process. */
    std::vector<int64_t> LocalCCSizes(nCC, 0);

    /* For each local process, get the locally stored vector of CC membership */
    std::vector<int64_t> localCC = vCC.GetLocVec();
    int64_t localcclen = vCC.LocArrSize();

    /* On each process, compute the local contributions of connected components members and
     * read sizes */
    for (int64_t i = 0; i < localcclen; ++i)
        LocalCCSizes[localCC[i]]++;

    int avelen = nCC / nprocs; /* size of local vector on every process except for last one */
    int lastlen = nCC - (avelen * (nprocs - 1)); /* above except only for last process */

    /* @recvcounts used by MPI_Reduce_scatter to know the local array sizes for each process (??) */
    std::vector<int> recvcounts(nprocs, avelen);
    recvcounts.back() = lastlen;

    int mylen = (myrank != nprocs - 1)? avelen : lastlen; /* local length of this process */

    int64_t *fillarrCC = new int64_t[mylen];

    /* Compute the global CC counts, and then distribute them according to CombBLAS distributed
     * vector constructor requirements */
    MPI_Reduce_scatter(LocalCCSizes.data(), fillarrCC, recvcounts.data(), MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

    /* Copy array into vector */
    std::vector<int64_t> fillvecCC(mylen);
    fillvecCC.assign(fillarrCC, fillarrCC + mylen);

    FullyDistVec<int64_t, int64_t> ccSizes(fillvecCC, fullWorld);

    delete[] fillarrCC;

    //std::string stringm = myoutput;
    //stringm += ".reads.per.contig.txt";
    //ccSizes.ParallelWrite(stringm, true);

    /* !!!! STARTING WITH nprocs == 1 !!!! */


    /* ReadProcessorAssignment(); */

    /* Once we have the read processor assignment (FillyDistVec<int64_t,
     * int64_t> assignment (??)), then we call the proper MPI collective for
     * distributing the reads, using DistributedFastaData and/or FastaData. Look
     * at KmerOps for ideas about this sort of exchange. */

    /* Once we have all the reads properly distributed, run contig assembly
     * locally on each process. If an entire contig is on one process, then the
     * entire contig will be assembled locally. If a contig happens to spill
     * over, then we wait for all local assemblies to complete, and then rerun
     * ReadProcessorAssignment (which maybe should be renamed then?) to get
     * assignment of partially assembled contigs (treated as reads maybe?) to
     * single processors */

}



// typedef ContigSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ContigSRing_t;

// std::vector<std::string>
// CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu,
//    PSpMat<dibella::CommonKmers>::DCCols* spSeq, std::shared_ptr<DistributedFastaData> dfd, int64_t nreads)
// {

   // float balance = S.LoadImbalance();
   // int64_t nnz   = S.getnnz();

   // std::ostringstream outs;
   // outs.str("");
   // outs.clear();
   // outs << "CreateContig::LoadBalance: " << balance << endl;
   // outs << "CreateContig::nonzeros: "    << nnz     << endl;
   // SpParHelper::Print(outs.str());

   // int64_t nCC = 0;
   // FullyDistVec<int64_t, int64_t> myLabelCC = CC(S, nCC);

   //	std::string myccoutput = myoutput + ".cc";
   //	myLabelCC.ParallelWrite(myccoutput, 1);

   // uint64_t nContig = myLabelCC.Reduce(maximum<int64_t>(), (int64_t)0);
   // nContig++; // because of zero based indexing for cluster

   // std::stringstream ncc;
   // ncc << "nContig: " << nContig << endl;
   // SpParHelper::Print(ncc.str());

   // First4Clust(myLabelCC);
   // HistCC(myLabelCC, nCC);

   // PrintCC(myLabelCC, nCC);

#ifdef MATRIXPOWER

    // PSpMat<dibella::CommonKmers>::MPI_DCCols T = S; // (T) traversal matrix
    // PSpMat<dibella::CommonKmers>::MPI_DCCols ContigM(S.getcommgrid()); // Contig matrix
    // dibella::CommonKmers defaultBVal;

    // Read vector is gonna multiply the matrix and create contig from there
    // FullyDistVec<int64_t, std::array<char, MAXCONTIGLEN>> ReadVector(S.getcommgrid());
    // char* nt; // NT initial value

    // CustomVectorEntry nt;
    /*
     * A vector of read ids that is my path [0, 1, 3, 67] (paths)
     * A vector of offset [10, 40, 50] (offsets) and offsets.size = paths.size -1
     * Offset 10 tells me that i have to cut the last 10 bases of read1 and concatenate them to read0, then 40 from read3 and concatenate them to read1 etc.
    */
    // FullyDistVec<int64_t, char*> ReadVector(S.getcommgrid(), nreads, nt);
    //	auto dcsc = spSeq->GetDCSC();

    // GGGG: fill the read vector with sequences
    // IT * ir; //!< row indices, size nz
    //	for (uint64_t i = 0; i < dcsc->nz; ++i)
    //	{
    //		int64_t lrid = dcsc->ir[i]; // local row idx
    //    seqan::Dna5String rseq = *(dfd->row_seq(lrid));

    //    std::array<char, MAXCONTIGLEN> crseq;
    //    std::string cpprseq;
    //    std::copy(begin(rseq), end(rseq), begin(cpprseq)); // @GGGG: this doesnt work

    //    char* crseq = new char[MAXCONTIGLEN];
    //    crseq = &cpprseq[0]; // C++14

    //    std::cout << rseq << std::endl;

    //   ReadVector.SetElement(lrid, crseq); // SetElement only work locally (owner)
    //	}

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

    // std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    // return myContigSet;
// }

// FullyDistVec<int64_t, dibella::CommonKmers> ReduceV(T.getcommgrid());
// FullyDistVec<int64_t, dibella::CommonKmers> ContigV(T.getcommgrid());

// ContigV =  T.Reduce(Row, ReduceMSR_t(), NullValue);
// ReduceV = nT.Reduce(Row, ReduceMSR_t(), NullValue);

