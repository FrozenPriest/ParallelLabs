#include "time.h"
#include "random"
#include <iostream>
#include <cstring>
#include "mpi.h"

void print_matrix(long double *matrix, int n)
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

void print_vector(long double *vector, int n)
{
    for (int i = 0; i < n; i++)
    {
        std::cout << vector[i] << " ";
    }
    std::cout << std::endl;
}
//returns vector
void matrix_mul_vector(long double *matrix_left, long double *vector, int n, long double *result)
{
    for (int row_index = 0; row_index < n; row_index += 1)
    {
        result[row_index] = 0;
        for (int column_index = 0; column_index < n; column_index++)
        {
            auto matrix_index = row_index * n + column_index;
            result[row_index] += matrix_left[matrix_index] * vector[column_index];
        }
    }
}
void vector_sub_vector(long double *vector_left, long double *vector_right, int n, long double *result)
{
    for (int index = 0; index < n; index += 1)
    {
        result[index] = vector_left[index] - vector_right[index];
    }
}

long double vector_norm(long double *vector, int n)
{
    long double sum = 0;
    for (int row_index = 0; row_index < n; row_index += 1)
    {
        sum += vector[row_index] * vector[row_index];
    }
    return sqrt(sum);
}

long double vector_max(long double *vector, int n)
{
    long double max = vector[0];

    for (int row_index = 1; row_index < n; row_index += 1)
    {
        if (max < abs(vector[row_index]))
            max = abs(vector[row_index]);
    }
    return max;
}
//mpiexec -n 1 $fileNameWithoutExt
//cd "c:\Users\boral\MPIProjects\lab1\lab3\" && g++ lab3.cpp -o lab3 -fopenmp -l msmpi -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
int main(int argc, char *argv[])
{
    //init mpi
    if (int rc = MPI_Init(&argc, &argv))
    {
        std::cout << "Ошибка запуска, выполнение остановлено " << std::endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    int rank;
    int numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(134932190);

    long double convergence_criteria = 1e-5;
    long double omega = 1.9;

    std::cout << "Enter matrix size(N): ";
    int N;
    std::cin >> N;

    long double *matrix = new long double[N * N];
    long double sum = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i < j)
            {
                matrix[i * N + j] = rand() % 58 + 1;
            }
            else if (i == j)
            {
                matrix[i * N + j] = 0;
            }
            else
            {
                matrix[i * N + j] = matrix[j * N + i];
            }
            sum += matrix[i * N + j];
        }
    }
    for (int i = 0; i < N; i++)
    {
        matrix[i * N + i] = sum + rand() % 100 + 1;
    }
    long double *vector = new long double[N];
    for (int i = 0; i < N; i++)
    {
        vector[i] = rand() % 100 + 1;
    }
    //example 1
    //  N = 4;
    //  matrix = new long double[N * N]{4, -1, -6, 0, -5, -4, 10, 8, 0, 9, 4, -2, 1, 0, -7, 5};
    //  vector = new long double[N]{2, 21, -12, -6};
    //  N = 2;
    //example 2
    // matrix = new long double[N * N]{1, 1, 1, 3};
    //vector = new long double[N]{3, 7};
    //example 3
    //  N = 6;
    // matrix = new long double[N * N]{9, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 7, -1, 0, 1, 1, 1, 1, 1, 1, 1, 7, 0, 1, 12, 1, 1 };
    // vector = new long double[N]{ 90, 15, 16, 19, -20, 10 };

    std::cout << "Initial matrix: \n";
    //print_matrix(matrix, N);
    std::cout << "Initial vector: \n";
    //print_vector(vector, N);

    ///////////////////////////////////
    double startwtime = 0;
    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }

    int step = 0;
    long double *phi = new long double[N]{0};
    long double *phi_old = new long double[N]{0};

    long double *phi_diff = new long double[N]{0};
    //resudal
    long double residal = 1e100;
    while (residal > convergence_criteria)
    {
        for (int row_index = 0; row_index < N; row_index += 1)
        {
            long double sigma = 0;
            for (int column_index = 0; column_index < row_index; ++column_index)
            {
                sigma += matrix[row_index * N + column_index] * phi[column_index] * omega;
            }
            for (int column_index = row_index + 1; column_index < N; ++column_index)
            {
                sigma += matrix[row_index * N + column_index] * phi_old[column_index] * omega;
            }
            phi[row_index] = (vector[row_index] * omega - sigma) / matrix[row_index * N + row_index] - phi_old[row_index] * (omega - 1);
        }
        vector_sub_vector(phi, phi_old, N, phi_diff);
        residal = vector_norm(phi_diff, N);

        memcpy(phi_old, phi, N * sizeof(long double));
        step += 1;
        // std::cout << "Error vector: ";
        // print_vector(vector_sub_vector(matrix_mul_vector(matrix, phi, N), vector, N), N);
        //std::cout << "Phi: ";
        //print_vector(phi, N);
        // std::cout << "Step " << step  << ", residal = " << residal << std::endl;
    }

    if (rank == 0)
    {
        double endwtime = MPI_Wtime();

        std::cout << "Result: ";
        print_vector(phi, N);
        std::cout << "Error: " << residal << std::endl;
        std::cout << (endwtime - startwtime) * 1000 << " ms." << std::endl;
    }

    delete phi;
    delete phi_old;
    delete phi_diff;
    delete matrix;
    delete vector;

    MPI_Finalize();
    return 0;
}