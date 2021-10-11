#include "time.h"
#include "random"
#include <iostream>

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
double *matrix_mul_vector(double *matrix_left, double *vector, int n)
{
    double *result = new double[n];

    for (int row_index = 0; row_index < n; row_index += 1)
    {
        for (int column_index = 0; column_index < n; column_index++)
        {
            auto matrix_index = row_index * n + column_index;
            result[row_index] += matrix_left[matrix_index] * vector[column_index];
        }
    }
    return result;
}
double *vector_sub_vector(double *vector_left, double *vector_right, int n)
{
    double *result = new double[n];
    for (int index = 0; index < n; index += 1)
    {
        result[index] = vector_left[index] - vector_right[index];
    }
    return result;
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
//mpiexec -n 1 $fileNameWithoutExt
//cd "c:\Users\boral\MPIProjects\lab1\lab3\" && g++ lab3.cpp -o lab3 -fopenmp -l msmpi -L "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64" -I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include"
int main()
{
    srand(time(nullptr));

    double convergence_criteria = 1e-7;
    double omega = 1.5;

    std::cout << "Enter matrix size(N): ";
    int N;
    std::cin >> N;

    double *matrix = new double[N * N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i <= j)
            {
                matrix[i * N + j] = rand() % 100;
            }
            else
            {
                matrix[i * N + j] = matrix[j * N + i];
            }
        }
    }
    double *vector = new double[N];
    for (int i = 0; i < N; i++)
    {
        vector[i] = rand() % 100;
    }
    //example 1
   // matrix = new double[N * N]{4, -1, -6, 0, -5, -4, 10, 8, 0, 9, 4, -2, 1, 0, -7, 5};
   // vector = new double[N]{2, 21, -12, -6};
    //example 2
   // matrix = new double[N * N]{1, 1, 1, 3};
   // vector = new double[N]{3, 7};

    std::cout << "Initial matrix: \n";
    print_matrix(matrix, N);
    std::cout << "Initial vector: \n";
    print_vector(vector, N);

    int step = 0;
    double *phi = new double[N]{0};
    //resudal
    double residal = abs(vector_norm(vector_sub_vector(matrix_mul_vector(matrix, phi, N), vector, N), N));
    while (residal > convergence_criteria)
    {
        for (int row_index = 0; row_index < N; row_index += 1)
        {
            double sigma = 0;
            for (int column_index = 0; column_index < N; column_index++)
            {
                auto matrix_index = row_index * N + column_index;
                if (row_index != column_index)
                {
                    sigma += matrix[matrix_index] * phi[column_index];
                }
            }
            phi[row_index] = (1 - omega) * phi[row_index] + (omega / matrix[row_index * N + row_index]) * (vector[row_index] - sigma);
        }
        residal = vector_norm(vector_sub_vector(matrix_mul_vector(matrix, phi, N), vector, N), N);
        step += 1;
        std::cout << "Step " << step << ", residal = " << residal << ", phi = ";
        print_vector(phi, N);
    }
    std::cout << "Result: ";
    print_vector(phi, N);
    return 0;
}