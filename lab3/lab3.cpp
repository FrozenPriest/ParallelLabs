#include "time.h"
#include "random"
#include <iostream>
#include <cstring>
#include "mpi.h"

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
//returns vector
void matrix_mul_vector(double *matrix_left, double *vector, int n, double *result)
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
//5000
//Error: 9.53317e-06
//352.359 ms. 310.76 ms. 457.948 ms. 1698.79 ms.
void vector_sub_vector(double *vector_left, double *vector_right, int n, double *result)
{
    for (int index = 0; index < n; index += 1)
    {
        result[index] = vector_left[index] - vector_right[index];
    }
}

double vector_norm(double *vector, int n)
{
    double sum = 0;
    for (int row_index = 0; row_index < n; row_index += 1)
    {
        sum += vector[row_index] * vector[row_index];
    }
    return sqrt(sum);
}

double vector_max(double *vector, int n)
{
    double max = vector[0];

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

    double convergence_criteria = 1e-8;
    double omega = 1.9;

    int N;

    if (rank == 0)
    {
        std::cout << "Enter matrix size(N): ";
        //std::cin >> N; 333.505 ms. 358.683 ms. 776.442 ms. 1057.27 ms.
        N = 1000; //3622.26 ms
    }
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double *matrix = new double[N * N];
    double sum = 0;
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
    double *vector = new double[N];
    for (int i = 0; i < N; i++)
    {
        vector[i] = rand() % 100 + 1;
    }

    if (rank == 0)
    {
    //    std::cout << "Initial matrix: \n";
    //    print_matrix(matrix, N);
    //    std::cout << "Initial vector: \n";
    //    print_vector(vector, N);
    }
    MPI_Bcast(matrix, N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vector, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    ///////////////////////////////////
    double startwtime = 0;
    if (rank == 0)
    {
        startwtime = MPI_Wtime();
    }

    int step = 0;
    double *phi = new double[N]{0};
    double *phi_old = new double[N]{0};

    double *phi_diff = new double[N]{0};
    //resudal
    double residal = 1e100;

    while (residal > convergence_criteria)
    {
        for (int row_index = 0; row_index < N; row_index += 1)
        {
            double sigma = 0;
            for (int column_index = 0 + rank; column_index < row_index; column_index += numprocs)
            {
                sigma += matrix[row_index * N + column_index] * phi[column_index] * omega;
            }
            for (int column_index = row_index + 1 + rank; column_index < N; column_index += numprocs)
            {
                sigma += matrix[row_index * N + column_index] * phi_old[column_index] * omega;
            }
            double sigma_reduced = 0;

            MPI_Reduce(&sigma, &sigma_reduced, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            
            if (rank == 0)
            {
                phi[row_index] = (vector[row_index] * omega - sigma_reduced) / matrix[row_index * N + row_index] - phi_old[row_index] * (omega - 1);
            }
            MPI_Bcast(phi + row_index, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        vector_sub_vector(phi, phi_old, N, phi_diff);
        residal = vector_norm(phi_diff, N);

        memcpy(phi_old, phi, N * sizeof(double));
        step += 1;
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