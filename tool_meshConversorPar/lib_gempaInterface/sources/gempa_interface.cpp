
#include <Partitioner.h>
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <algorithm>

using namespace std;

extern "C" {
    void gempa_do_partition(int numNodes, int numRanks, double coord_x[], double coord_y[], double coord_z[], int weights[], int partition[]);

    }

void gempa_do_partition(int numNodes, int numRanks, double coord_x[], double coord_y[], double coord_z[], int weights[], int partition[])
{
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    gempa::Partitioner P1;

    P1.sdim = 3;
    P1.N = numNodes;
    P1.npart = numRanks;
    P1.maxIter = 100;
    P1.isParallel = true;
    P1.isWeighted = true;

    P1.coords.resize(P1.N * P1.sdim);
    P1.weights.resize(P1.N);

    if(!mpi_rank) std::cout << "# Doing GEMPA partitioning..." << std::endl;

    for (int i = 0; i < P1.N; ++i)
    {
        int index = i * P1.sdim;

        P1.coords[index+0] = coord_x[i];
        P1.coords[index+1] = coord_y[i];
        P1.coords[index+2] = coord_z[i];

        P1.weights[i] = weights[i];
    }

    // Partition and check
    P1.part();
    double balance = P1.checkBalance();
    bool isOrdered = P1.checkOrdering();

    for (int i = 0; i < P1.N; ++i)
    {
      //cout << "["<<i<<"] "<<P1.partition[i] << std::endl;
      partition[i] = P1.partition[i];
    }
    if(!mpi_rank) std::cout << "# GEMPA partitioning done!" << std::endl;
}

