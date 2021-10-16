#include "time.h"
#include "random"
#include "mpi.h"
#include <iostream>
#include "gauss.hpp"

void print_matrix(double *matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << matrix[i * n + j] << " ";
        }
        std::cout << std::endl;
    }
}

void print_vector(double *vector, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << vector[i] << " ";
    }
    std::cout << std::endl;
}

void matrix_mul_matrix(double *left, double *right, double *ans, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            ans[i * n + j] = 0;
            for (int k = 0; k < n; k++)
            {
                ans[i * n + j] += left[i * n + k] * right[k * n + j];
            }
        }
    }
}

//returns swap count
void lu_decomposition(double *matrix, int n, int rank, int numprocs)
{
    MPI_Barrier(MPI_COMM_WORLD);
    for (int row_index = 0; row_index < n - 1; row_index++)
    {
        int pseudorank = row_index % numprocs;
        if (pseudorank == rank)
        {
            for (int i = row_index + 1; i < n; i++)
            {
                matrix[i * n + row_index] = matrix[i * n + row_index] / matrix[row_index * n + row_index];

                for (int z = row_index + 1; z < n; z++)
                {
                    matrix[i * n + z] = matrix[i * n + z] - (matrix[i * n + row_index] * matrix[row_index * n + z]);
                }
            }
        }
        MPI_Bcast(matrix, n * n, MPI_DOUBLE, pseudorank, MPI_COMM_WORLD);
    }
}

double *extract_L_matrix(double *matrix, int n)
{
    double *L = new double[n * n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j > i)
            {
                L[i * n + j] = 0;
            }
            else if (i == j)
            {
                L[i * n + j] = 1;
            }
            else
            {
                L[i * n + j] = matrix[i * n + j];
            }
        }
    }

    return L;
}

double *extract_U_matrix(double *matrix, int n)
{
    double *U = new double[n * n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j >= i)
            {
                U[i * n + j] = matrix[i * n + j];
            }
            else
            {
                U[i * n + j] = 0;
            }
        }
    }

    return U;
}

double *solve_y(double *L, double *vector, int n, int rank, int numprocs)
{
    double *solution = new double[n]{0};

    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += L[i * n + j] * solution[j];
        }
        solution[i] = vector[i] - sum;
    }

    return solution;
}

double *solve_x(double *U, double *y, int n, int rank, int numprocs)
{
    double *solution = new double[n]{0};

    for (int i = n-1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += U[i * n + j] * solution[j];
        }
        solution[i] = (y[i] - sum) / U[i * n + i];
    }

    return solution;
}

double *solve_lu_decomposition(double *matrix, double *vector, int n, int rank, int numprocs)
{
    lu_decomposition(matrix, n, rank, numprocs);

    double *L = extract_L_matrix(matrix, n);
    double *U = extract_U_matrix(matrix, n);

    double *lu = new double[n * n];
    matrix_mul_matrix(L, U, lu, n);

    if (rank == 0)
    {
        std::cout << "L:\n";
        print_matrix(L, n);
        std::cout << "U:\n";
        print_matrix(U, n);
        std::cout << "L * U = \n";
        print_matrix(lu, n);
    }

    double *y = solve_y(L, vector, n, rank, numprocs);
    double *solution = solve_x(U, y, n, rank, numprocs);

    return solution;
} //5 = -1.17584 -3.6856 1.52991 -4.08969 2.00373

//mpiexec -n 1 $fileNameWithoutExt
//cd "c:\Users\boral\MPIProjects\lab1\lab5\" && g++ lab5.cpp -o lab5 -fopenmp -l msmpi -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"

int main()
{
    MPI_Init(nullptr, nullptr);
    int rank;
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(32549201);

    int N;
    if (rank == 0)
    {
        std::cout << "Enter matrix size(N): ";
        std::cin >> N;
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //matrix generation
    double *matrix = new double[N * N];
    double sum = 0;
    do
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                matrix[i * N + j] = rand() % 100 - 50;
            }
        }
    } while (std::abs(get_gauss_determinant(matrix, N)) < 1e-5);

    //vector generation
    double *vector = new double[N];
    for (int i = 0; i < N; i++)
    {
        vector[i] = rand() % 100 + 1;
    }
    if (rank == 0)
    {
        std::cout << "Initial matrix: \n";
        print_matrix(matrix, N);
        std::cout << "Initial vector: \n";
        print_vector(vector, N);
    }
    double begin = MPI_Wtime();

    double *answer = solve_lu_decomposition(matrix, vector, N, rank, numprocs);

    double end = MPI_Wtime();
    if (rank == 0)
    {
        if (answer == nullptr)
        {
            std::cout << "Answer is non-existant" << std::endl;
            return 0;
        }
        std::cout << "Answer:\n";
        print_vector(answer, N);
        std::cout << "Time: " << end - begin << " s." << std::endl;
    }

    MPI_Finalize();
    return 0;
}