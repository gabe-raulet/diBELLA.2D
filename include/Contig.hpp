// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "CC.h"

// #define MATRIXPOWER
#define APSPDEF
#define MAXPATHLEN 5000

#ifdef APSPDEF // All-Pairs Shortest Path (not great for big graphs, I wanna use Seidel's algorithm but modify it for sparse matrices)

#include <algorithm>
#define INF 10000

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

void PrintAPSP(ContigMatType& D, int nvertex)
{
    for (unsigned i = 0; i < nvertex; ++i)
        for (unsigned j = 0; j < nvertex; ++j)
        {
            Overlap entry = D[i][j];

            if(entry.len == INF) printf("%7s\t", "INF");
            else printf("%7d\t", entry.len); 

            if(j == nvertex - 1) printf("\n");
        }
    printf("\n");
}

void APSP(ContigMatType& CustomS, int nvertex, int nnz) 
{
    ContigMatType D = CustomS; // the distances matrix (copy constructor)
    unsigned i, j, k;
 
    for (k = 0; k < nvertex; k++)
        for (i = 0; i < nvertex; i++)
            for (j = 0; j < nvertex; j++)
            {
                Overlap ij = D[i][j];
                Overlap kj = D[k][j];
                Overlap ik = D[i][k];

                if((ij.len > (ik.len + kj.len)) && DirectionalityCheck(ik.dir, kj.dir, ij.dir))
                {   
                    /* Debug print */
                    // std::cout << "(i, j) = (" << i << ", " << j << ")" << " | (i, k) = (" << i << ", " << k << ")" << " | (k, j) = (" << k << ", " << j << ")" << std::endl;
                    // std::cout << "(i, j) = " << ij.dir << " | (i, k) = " << ik.dir << " | (k, j) = " << kj.dir << std::endl;
                    // std::cout << std::endl;
                    
                    ij.len = ik.len + kj.len;

                    seqan::append(ik.seq, kj.seq);
                    assert(ij.len == length(ik.seq));

                    ij.seq = ik.seq; // + kj.seq;
                }

                D[i][j] = ij;
            }
    PrintAPSP(D, nvertex);
}

#endif

/*! Namespace declarations */
using namespace combblas;

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
#ifdef APSPDEF

    // The string matrix using vec of vec there's a better way to do this using sparsity
    ContigMatType CustomS;
    
    CustomS.resize(nreads);
    for(int i = 0; i < nreads; i++)
    {
        CustomS[i].resize(nreads);
    }

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
            int overhang = cks->overhang >> 2;

            // Get row/col sequences paying attention to the consistency of V/H
            seqan::Dna5String seqH = *(dfd->row_seq(dcsc->ir[j])); // extract row sequence
            seqan::Dna5String seqV = *(dfd->row_seq(dcsc->jc[i])); // extract col sequence

            seqan::Dna5String contigsubstr;

            /* Debug print */
            // std::cout << dcsc->ir[j]+1 << " " << dcsc->jc[i]+1 << " " << dir << " " << overhang << std::endl;

            uint rlenV = length(seqV);
            uint rlenH = length(seqH);

            // Use case switch!
            if (dir == 1)
            {
                int endpV = rlenV - overhang;
                contigsubstr = seqan::suffix(seqV, endpV); // seqH entering in seqV so we extract the ending part of seqV (fwd strand)
            }
            else if (dir == 2)
            {
                int endpH = rlenH - overhang;
                contigsubstr = seqan::suffix(seqH, endpH); // seqH entering in seqV so we extract the ending part of seqV (fwd strand)
            }
            else if (dir == 0)
            {
                int begpV = overhang;
                contigsubstr = seqan::prefix(seqV, begpV+1); // 2nd parameter excludes end position
                seqan::Dna5StringReverseComplement twin(contigsubstr);
                std::get<3>(mattuples[z]) = twin;
            }
            else if (dir == 3)
            {
                int endpV = rlenV - overhang;
                contigsubstr = seqan::suffix(seqV, endpV);
            }
            
            std::get<2>(mattuples[z]) = dir;
            if(dir != 0) std::get<3>(mattuples[z]) = contigsubstr; // Everything should technically be consistent, if correct

            /* Debug print */
            // std::cout << seqH << std::endl;
            // std::cout << seqV << std::endl;
            // std::cout << std::get<3>(mattuples[z]) << std::endl;
            // std::cout << std::endl;

            ++z;
		}
	}

	assert(z == nnz);
	std::cout << "Local nnz count: " << nnz << std::endl;

    // The entire matrix S is initialized to INF with 0 on the diagonal
    for (unsigned i = 0; i < nreads; ++i)
        for (unsigned j = 0; j < nreads; ++j)
        {
            Overlap entry;

            entry.seq = "";
            entry.dir = -1; // if INF (zeros) direction must be "indefined", 0 is a direction in our bidirected graph

            if(i != j) entry.len = INF;
            else entry.len = 0; // this is important dude (zeros on the diagonal must be actual zeros) ---Question how do I modify this for sparsity?

            CustomS[i][j] = entry;
        }

    // The matrix stores the nonzeros
    for(int z = 0; z < nnz; z++)
    {
        Overlap nnzentry;

        uint i = std::get<0>(mattuples[z]);
        uint j = std::get<1>(mattuples[z]);

        assert(i != j);

        nnzentry.dir = std::get<2>(mattuples[z]);
        nnzentry.seq = std::get<3>(mattuples[z]);
        nnzentry.len = length(nnzentry.seq);

        CustomS[i][j] = nnzentry;        
    }
     
    // 3) APSP (inefficient for big matrix but let's see if contig makes sense first)
    PrintAPSP(CustomS, nreads);
    APSP(CustomS,nreads, nnz); 
    
    printf("APSP is done!\n");

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