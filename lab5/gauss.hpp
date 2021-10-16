#include <iostream>

//замена столбца матрицы вектором свободных элементов
void swap_column(double *matrix, double *free_term_column, int i, int n)
{
    for (int col_element_index = 0; col_element_index < n; ++col_element_index)
    {
        double temp = matrix[col_element_index * n + i];
        matrix[col_element_index * n + i] = free_term_column[col_element_index];
        free_term_column[col_element_index] = temp;
    }
}

//замена двух строк
void swap_row(double *matrix, int i, int j, int n)
{
    for (int row_element_index = 0; row_element_index < n; ++row_element_index)
    {
        double temp = matrix[i * n + row_element_index];
        matrix[i * n + row_element_index] = matrix[j * n + row_element_index];
        matrix[i * n + row_element_index] = temp;
    }
}

//поиск максимального элемента в столбце
int col_max(double *matrix, int col, int n)
{
    double max = std::abs(matrix[col * n + col]);
    int maxPos = col;
    for (int i = col + 1; i < n; ++i)
    {
        double element = std::abs(matrix[i * n + col]);
        if (element > max)
        {
            max = element;
            maxPos = i;
        }
    }
    return maxPos;
}

//returns swap count
int triangulation(double *matrix, int n)
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
            double mul = -matrix[j * n + i] / matrix[i * n + i];
            for (int k = i; k < n; ++k)
            {
                matrix[j * n + k] += matrix[i * n + k] * mul;
            }
        }
    }
    return swap_count;
}
//4.17148  2.38426  1.42846  1.16603
double *copy_array(double *array, int n)
{
    double *copy = new double[n];
    for (size_t i = 0; i < n; i++)
    {
        copy[i] = array[i];
    }
    return copy;
}

//определитель по гауссу
double get_gauss_determinant(double *matrix, int n)
{
    double *matrix_copy = copy_array(matrix, n * n);

    int swap_count = triangulation(matrix_copy, n);
    double determinant = 1;

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