#include "time.h"
#include "random"
#include <omp.h>
#include <iostream>

void print_array(int *array, int size)
{
    std::cout << "Array:";
    for (size_t i = 0; i < size; i++)
    {
        std::cout << " " << array[i];
    }
    std::cout << "." << std::endl;
}

#define TASK_SIZE 100

//Lomuto partition scheme
int partition(int *array, int low, int high)
{
    auto pivot = array[high];
    int pivotIndex = low - 1;

    for (int cur = low; cur <= high; cur++)
    {
        if (array[cur] <= pivot)
        {
            pivotIndex++;
            //swap a[pivotIndex] a[cur]
            auto temp = array[pivotIndex];
            array[pivotIndex] = array[cur];
            array[cur] = temp;
        }
    }
    return pivotIndex;
}

void quicksort(int *array, int low, int high)
{
    if (low >= 0 && low < high)
    {
        auto pivotIndex = partition(array, low, high);
//std::cout << "Info: " << low << ", " << high << "." << std::endl;
#pragma omp task shared(array) if (high - pivotIndex > TASK_SIZE)
        quicksort(array, low, pivotIndex - 1);

#pragma omp task shared(array) if (high - pivotIndex > TASK_SIZE)
        quicksort(array, pivotIndex + 1, high);
    }
}

bool checkSorted(int *array, int size)
{
    for (int i = 1; i < size; i++)
    {
        if (array[i] < array[i - 1])
        {
            return false;
        }
    }
    return true;
}

int main()
{
    srand(9373937);
    std::cout << "Enter array size: ";
    int size;
    std::cin >> size;

    std::cout << "Enter thread count: ";
    int thread_count;
    std::cin >> thread_count;

    int *array = new int[size];
    for (int i = 0; i < size; i++)
    {
        array[i] = rand() % 100;
    }
    // print_array(array, size);

    //openmp settings
    omp_set_dynamic(0);
    omp_set_num_threads(thread_count);

    double begin = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp single
        quicksort(array, 0, size - 1);
    }
    double end = omp_get_wtime();

    //print_array(array, size);
    std::cout << "Time: " << end - begin << " ms." << std::endl;
    std::cout << "Sort is confirmed!" << std::endl;

    delete array;
    return 0;
}

/*
    Check when making report:
    https://stackoverflow.com/questions/16007640/openmp-parallel-quicksort
    https://pro-prof.com/archives/1220
    
*/