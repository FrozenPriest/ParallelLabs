#include "time.h"
#include "random"
#include <omp.h>
#include <iostream>

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

//замена столбца матрицы вектором свободных элементов
void swap_column(long double *matrix, long double *free_term_column, int i, int n)
{
    for (int col_element_index = 0; col_element_index < n; ++col_element_index)
    {
        long double temp = matrix[col_element_index * n + i];
        matrix[col_element_index * n + i] = free_term_column[col_element_index];
        free_term_column[col_element_index] = temp;
    }
}

//замена двух строк
void swap_row(long double *matrix, int i, int j, int n)
{
    for (int row_element_index = 0; row_element_index < n; ++row_element_index)
    {
        long double temp = matrix[i * n + row_element_index];
        matrix[i * n + row_element_index] = matrix[j * n + row_element_index];
        matrix[i * n + row_element_index] = temp;
    }
}

//поиск максимального элемента в столбце
int col_max(long double *matrix, int col, int n)
{
    long double max = std::abs(matrix[col * n + col]);
    int maxPos = col;
    for (int i = col + 1; i < n; ++i)
    {
        long double element = std::abs(matrix[i * n + col]);
        if (element > max)
        {
            max = element;
            maxPos = i;
        }
    }
    return maxPos;
}

//returns swap count
int triangulation(long double *matrix, int n)
{
    int swap_count = 0;
    for (int i = 0; i < n - 1; ++i)
    {
        int imax = col_max(matrix, i, n);
        if (i != imax)
        {
            swap_row(matrix, i, imax, n);
            ++swap_count;
        }
        //#pragma omp parallel for
        for (int j = i + 1; j < n; ++j)
        {
            long double mul = -matrix[j * n + i] / matrix[i * n + i];
            for (int k = i; k < n; ++k)
            {
                matrix[j * n + k] += matrix[i * n + k] * mul;
            }
        }
    }
    return swap_count;
}
//4.17148  2.38426  1.42846  1.16603
long double *copy_array(long double *array, int n)
{
    long double *copy = new long double[n];
    for (size_t i = 0; i < n; i++)
    {
        copy[i] = array[i];
    }
    return copy;
}

//определитель по гауссу
long double get_gauss_determinant(long double *matrix, int n)
{
    long double *matrix_copy = copy_array(matrix, n * n);

    int swap_count = triangulation(matrix_copy, n);
    long double determinant = 1;

    if (swap_count % 2 == 1)
    {
        determinant = -1;
    }
    for (int i = 0; i < n; i++)
    {
        determinant *= matrix_copy[i * n + i];
    }

    delete matrix_copy;
    return determinant;
}

long double *solve_crammer(long double *matrix, long double *vector, int n)
{
    long double determinant = get_gauss_determinant(matrix, n);
    if (std::abs(determinant) < 1e-5)
    {
        return nullptr; //no answer
    }
    long double *solution = new long double[n]{0};
#pragma omp parallel
    {
        long double *matrix_copy = copy_array(matrix, n * n);
        long double *vector_copy = copy_array(vector, n);

#pragma omp for
        for (int col_index = 0; col_index < n; col_index++)
        {
            swap_column(matrix_copy, vector_copy, col_index, n);
            solution[col_index] = get_gauss_determinant(matrix_copy, n) / determinant;
            swap_column(matrix_copy, vector_copy, col_index, n);
        }

        delete matrix_copy;
        delete vector_copy;
    }

    return solution;
}

int main()
{
    srand(32549201);

    int N;
    std::cout << "Enter matrix size(N): ";
    std::cin >> N;

    //matrix generation
    long double *matrix = new long double[N * N];
    long double sum = 0;
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
    long double *vector = new long double[N];
    for (int i = 0; i < N; i++)
    {
        vector[i] = rand() % 100 + 1;
    }

    std::cout << "Initial matrix: \n";
    //print_matrix(matrix, N);
    std::cout << "Initial vector: \n";
    //print_vector(vector, N);

    int thread_count;
    std::cout << "Enter thread count: ";
    std::cin >> thread_count;

    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);

    double begin = omp_get_wtime();

    long double *answer = solve_crammer(matrix, vector, N);

    double end = omp_get_wtime();

    if (answer == nullptr)
    {
        std::cout << "Answer is non-existant" << std::endl;
        return 0;
    }
    print_vector(answer, N);
    std::cout << "Time: " << end - begin << " ms." << std::endl;

    return 0;
}