#include <stdio.h>
#include <iostream>
#include <string>
#include <string.h>
#include <algorithm>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;
#define INF 10000

/* Data struct for the nonzeros */

typedef struct
{
    int len;
    int dir; // I can use short (not ushort, -1 is undefined)
    std::string seq;
} Overlap;

/* Bit encoded directionality function */

void tobinary(ushort n, int* arr) 
{ 
    int nbit = 2;
    for(int i = 0; i < nbit; i++)
    { 
        arr[i] = n % 2; 
        n = n / 2; 
    }
}

// if false, dir (ij.dir) isn't modified
bool DirectionalityCheck(int dir1, int dir2, int& dir)
{
    int rbit, lbit;
    int start, end;

    int mybin1[2] = {0, 0}; // 0 1
    int mybin2[2] = {0, 0}; // 0 0

    if(dir1 != 0) tobinary(dir1, mybin1);
    if(dir2 != 0) tobinary(dir2, mybin2);

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

/* Directionality check main function */

// If one of the two is INF, the condition ij > ik + kj is not gonna be satified
bool isDirOk(int dir1, int dir2, int& dir)
{
    if(DirectionalityCheck(dir1, dir2, dir)) return true;
    else return false;
}

/* Print function for result */

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

/* Flyod-Warshall function */
 
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

/* Dijkstra function */

// Dijkstra might be sufficient (and more efficient) if I know where to start from (i.e., one read per contig)
void Dijkstra()
{
    /* Code here */
}

/* The main function */

// @GGGG: Can I incorporate the string into the output file from ELBA?
int main(int argc, char **argv)
{
    // int nvertex = std::stoi(argv[1]); // number of sequences
    // int nnz = std::stoi(argv[2]);     // number of nonzeros
    // char* input = argv[3];            // name of string matrix file

    int nvertex = 4;
    int nnz = nvertex - 1;

    // the string matrix S
    mapped_matrix<Overlap> S(nvertex, nvertex, nnz);

    // for (auto iter1 = S.begin1(); iter1 != S.end1(); ++iter1)
    //     for (auto iter2 = S.begin2(); iter2 != S.end2(); ++iter2)

    // GGGG: iterate over the string matrix market input and store nozneros into the matrix S
    for (unsigned i = 0; i < S.size1(); ++i)
        for (unsigned j = 0; j < S.size2(); ++j)
        {
            Overlap entry;
            entry.seq = "";

            entry.dir = -1; // if INF (zeros) direction must be "indefined", 0 is a direction in our bidirected graph

            if(i != j) entry.len = INF;
            else entry.len = 0; // this is important dude (zeros on the diagonal must be actual zeros) ---Question how do I modify this for sparsity?

            S(i, j) = entry;
        }

    PrintAPSP(S, nvertex);

    Overlap entry;

    /* Graph */
    entry.len = 5;
    entry.dir = 1; // >--->
    entry.seq = "CANE";
    S(0, 1) = entry;
    entry.dir = 0; // >---<
    entry.len = 3;
    S(1, 2) = entry;
    entry.len = 1;
    entry.dir = 1; // >--->
    S(2, 3) = entry;

    /* Let's add an extra edge to check length check */
    entry.seq = "GATTO";
    entry.len = 10;
    entry.dir = 1; // >--->
    S(0, 3) = entry;

    entry.len = 5;
    entry.dir = 2;
    entry.seq = "ENAC";
    S(1, 0) = entry;
    entry.len = 3;
    entry.dir = 0;
    S(2, 1) = entry;
    entry.len = 1;
    entry.dir = 2;
    S(3, 2) = entry;

    entry.seq = "OTTAG";
    entry.len = 10;
    entry.dir = 2;
    S(3, 0) = entry;

    PrintAPSP(S, nvertex);
    APSP(S, nvertex, nnz);

    return 0;
}

/* Giulia's TODOS */

/* 1. Floyd-Warshall for sequence concatenation (x) and path checking (x)
        - Path checking should only be done on the real nonzeros, INF should be ignored?
        - TODO: Check if string reversion is a thing of think if I'll have everything in place from my output (in terms of string concatenation) <--
        - TODO: Get a toy contig from E. coli HiFi for the meeting on Wed (i.e., I need to figure out the string bit)
        - TODO: Dijkstra version
 * 2. Seidel's algorithm for sequence concatenation and path checking
 * 3. Seidel's algortihm for sparse matrices (1 and 2 assumed dense matrices)
 * 4. Parallel Seidel's algorithm for sparse matrices
*/
