#include "time.h"
#include "random"
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
        std::cout << "Info: " << low << ", " << high << "." << std::endl;
        quicksort(array, low, pivotIndex - 1);
        quicksort(array, pivotIndex + 1, high);
    }
}

int main()
{
    srand(time(nullptr));
    std::cout << "Enter array size: ";
    int size;
    std::cin >> size;

    int *array = new int[size];
    for (int i = 0; i < size; i++)
    {
        array[i] = rand() % 100;
    }
    print_array(array, size);
    quicksort(array, 0, size - 1);
    print_array(array, size);

    return 0;
}