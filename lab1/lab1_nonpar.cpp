#include <iostream>
#include <windows.h>
#include <vector>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
{
  if (int rc = MPI_Init(&argc, &argv))
  {
    cout << "Ошибка запуска, выполнение остановлено " << endl;
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  for (int z = 290; z < 345; z += 5)
  {
    double time = 0.0;
    for (int jk = 0; jk < 100; jk++)
    {
      //matrix & vector
      const size_t matrix_width = 1000;
      const size_t matrix_height = 325;

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

      //timer ready
      auto startwtime = MPI_Wtime();

      double *t_res = new double[matrix_height];

      //main action
      for (int row_index = 0; row_index < matrix_height; row_index += 1)
      {
        for (int column_index = 0; column_index < matrix_width; column_index++)
        {
          auto matrix_index = row_index * matrix_width + column_index;
          t_res[row_index] += matrix[matrix_index] * mul_vector[column_index];
        }
      }

      double endwtime = MPI_Wtime();
      time += (endwtime - startwtime) * 1000;
      //cout << "Time: " << (endwtime - startwtime) * 1000 << " ms" << endl;
      delete matrix;
      delete mul_vector;
      delete t_res;
    }
    time /= 100.0;
    cout << time << endl;
  }
  MPI_Finalize();
  return 0;
}
