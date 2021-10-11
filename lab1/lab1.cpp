#include "mpi.h"
#include <iostream>
#include <windows.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  //init mpi
  if (int rc = MPI_Init(&argc, &argv))
  {
    cout << "Ошибка запуска, выполнение остановлено " << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }

  //matrix & vector
  size_t matrix_width = 1000;
  size_t matrix_height = 325;

  double *matrix = new double[matrix_width * matrix_height];
  double *mul_vector = new double[matrix_height];

  //mockup data
  for (size_t i = 0; i < matrix_height * matrix_width; i++)
  {
    matrix[i] = i;
  }
  for (size_t i = 0; i < matrix_height; i++)
  {
    mul_vector[i] = i * 4 + 7;
  }

  int rank;
  int numprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*  if (rank == 0)
  {
    cout << "Matrix:\n";
    for (int i = 0; i < matrix_height; i++)
    {
      for (int j = 0; j < matrix_width; j++)
      {
        cout << matrix[i * matrix_width + j] << " ";
      }
      cout << endl;
    }
    cout << "Vector:" << endl;
    for (int j = 0; j < matrix_height; j++)
    {
      cout << mul_vector[j] << endl;
    }
    cout << "The number of processes: " << numprocs << endl;
  } */

  //timer ready
  double startwtime;
  if (rank == 0)
  {
    startwtime = MPI_Wtime();
  }

  double *t_res = new double[matrix_height];

  //main action
  for (int row_index = rank; row_index < matrix_height; row_index += numprocs)
  {
    for (int column_index = 0; column_index < matrix_width; column_index++)
    {
      auto matrix_index = row_index * matrix_width + column_index;
      t_res[row_index] += matrix[matrix_index] * mul_vector[column_index];
    }
  }
  //reduce to main
  double *result = new double[matrix_height];
  MPI_Reduce(t_res, result, matrix_height, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  //print results
  if (rank == 0)
  {
    double endwtime = MPI_Wtime();

    /* cout << "Result:" << endl;
    for (int i = 0; i < matrix_height; i++)
    {
      cout << result[i] << endl;
    } */
    cout << (endwtime - startwtime) * 1000 << endl;
  }

  MPI_Finalize();
  return 0;
}
