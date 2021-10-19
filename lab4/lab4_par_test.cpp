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

void vector_sub_vector(long double *vector_left, long double *vector_right, long double *result, int n)
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

void matrix_mul_vector(long double *matrix_left, long double *vector, long double *result, int n)
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

//returns swap count
void triangulation(long double *matrix, int n)
{
    for (int i = 0; i < n - 1; ++i)
    {
        //#pragma omp parallel for
        for (int j = i + 1; j < n; ++j)
        {
            long double mul = -matrix[j * n + i] / matrix[i * n + i];
            matrix[j * n + i] = 0;
            for (int k = i + 1; k < n; ++k)
            {
                matrix[j * n + k] += matrix[i * n + k] * mul;
            }
        }
    }
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

    triangulation(matrix_copy, n);
    
    long double determinant = 1;

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

    for (int i = 1; i < 40; i++)
    {
        int N = i * 20;
        for (int thread_count = 1; thread_count < 9; thread_count++)
        {

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
            //print_vector(answer, N);
            //std::cout << "Time: " << end - begin << " s." << std::endl;

            long double *check_result = new long double[N]{0};
            matrix_mul_vector(matrix, answer, check_result, N);
            vector_sub_vector(vector, check_result, check_result, N);
            //print_vector(check_result, N);
            long double error = vector_norm(check_result, N);
            //std::cout << "Error: " << error << std::endl;

            std::cout << N << ";" << thread_count << ";" << end - begin << ";" << error << std::endl;

            delete check_result;
            delete vector;
            delete matrix;
        }
    }
    return 0;
}