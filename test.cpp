/*************************************************************************
	> File Name: test.cpp
	> Author: Binbin Li
	> Mail:
	> Created Time: Sun 02 Oct 2016 03:20:27 PM EDT
 ************************************************************************/
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

void communicate(
  const int& bot, const int& top, const int& right, const int& left,
  const int& squareRowId, const int& squareColId, double **localArray,
  const int& width, const int& rank, const int& size);
void update(int& rank, double **localArray, double **prevArray, const int& top,
    const int& bot, const int& left, const int& right, const int& m, const int& n);

int main(int argc, char** argv) {
  int rank, size;
  MPI_Init(0, 0);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Compute the size of each square
  int m = atoi(argv[1]), n = atoi(argv[2]);
  int width = (int) sqrt(size);
  int rowSize = m / width, colSize = n / width;
  if (rowSize * width < m) ++rowSize;
  if (colSize * width < n) ++colSize;
  int squareRowId = rank / width, squareColId = rank % width;

  // Compute upperleft index and bottomright index
  int left = squareColId * colSize;
  int right = min((squareColId + 1) * colSize - 1, n - 1);
  int top = squareRowId * rowSize;
  int bot = min((squareRowId + 1) * rowSize - 1, m - 1);

  // Initialize the array
  double **localArray = new double*[bot - top + 3];
  double **prevArray = new double*[bot - top + 3];
  for (int i = 0; i < bot - top + 3; ++i) {
    localArray[i] = new double [right - left + 3];
    prevArray[i] = new double [right - left + 3];
    for (int j = 0; j < right - left + 3; ++j) {
      localArray[i][j] = 0;
    }
    memcpy(prevArray[i], localArray[i], sizeof(double) * (right - left + 3));
  }
  vector<double> jVal(right - left + 1, 0);
  int iIdx, jIdx;
  for (int i = top; i <= bot; ++i) {
    double a = i * sin(i);
    iIdx = i - top + 1;
    if (i == top) {
      for (int j = left; j <= right; ++j) {
        jVal[j - left] = j * cos(j);
        localArray[iIdx][j - left + 1] = a + jVal[j - left] + sqrt(i + j);
      }
    } else {
      for (int j = left; j <= right; ++j) {
        localArray[iIdx][j - left + 1] = a + jVal[j - left] + sqrt(i + j);
      }
    }
  }
  double startTime = MPI_Wtime();

  // Update in 10 iterations
  for (int iteration = 1; iteration <= 10; ++iteration) {
    // Communicate with neighbors
    communicate(bot, top, right, left, squareRowId, squareColId, localArray, width, rank, size);
    // Update the array
    update(rank, localArray, prevArray, top, bot, left, right, m, n);
    //MPI_Barrier(MPI_COMM_WORLD);
  }

  double sum = 0;
  double squareSum = 0;
  for (int i = 1; i < bot - top + 2; ++i) {
    for (int j = 1; j < right - left + 2; ++j) {
      sum += localArray[i][j];
      squareSum += localArray[i][j] * localArray[i][j];
    }
  }
  for (int i = 0; i < bot - top + 3; ++i) {
    delete [] localArray[i];
    delete [] prevArray[i];
  }
  delete [] localArray;
  delete [] prevArray;

  // Send result to processor 0
  if (rank != 0) {
    MPI_Request request;
    MPI_Isend(&sum, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &request);
    MPI_Isend(&squareSum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &request);
  }
  else {
    double totalSum = sum, totalSquareSum = squareSum;
    MPI_Status status;
    for (int i = 1; i < size; ++i) {
      MPI_Recv(&sum, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(&squareSum, 1, MPI_DOUBLE, i , 1, MPI_COMM_WORLD, &status);
      totalSum += sum;
      totalSquareSum += squareSum;
    }
    cout << "Sum: " << totalSum << endl;
    cout << "Square Sum: " << totalSquareSum << endl;
  }
  if (rank == 0) cout << "Time: " << MPI_Wtime() - startTime << endl << endl;;

  MPI_Finalize();
  return 0;
}

void communicate(
  const int& bot, const int& top, const int& right, const int& left,
  const int& squareRowId, const int& squareColId, double **localArray,
  const int& width, const int& rank, const int& size) {
  MPI_Status status;
  int height = bot - top + 3, length = right - left + 3;
  // Send to below in the even row
  if (squareRowId % 2 == 0 && rank + width < size) {
    MPI_Ssend(localArray[bot - top + 1], length, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD);
  }
  // Receive from upper in the odd row
  if (squareRowId % 2 == 1) {
    MPI_Recv(localArray[0], length, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &status);
  }
  // Send to below in the odd row
  if (squareRowId % 2 == 1 && squareRowId != width - 1) {
    MPI_Ssend(localArray[bot - top + 1], length, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD);
  }
  // Receive from upper in the even row
  if (squareRowId % 2 == 0 && squareRowId != 0) {
    MPI_Recv(localArray[0], length, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &status);
  }
  // Send to upper in the even row
  if (squareRowId % 2 == 0 && squareRowId != 0) {
    MPI_Ssend(localArray[1], length, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD);
  }
  // Receive from below in the odd row
  if (squareRowId % 2 == 1 && squareRowId != width - 1) {
    MPI_Recv(localArray[height - 1], length, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &status);
  }
  // Send to upper in the odd row
  if (squareRowId % 2 == 1) {
    MPI_Ssend(localArray[1], length, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD);
  }
  // Receive from below in the even row
  if (squareRowId % 2 == 0 && squareRowId != width - 1) {
    MPI_Recv(localArray[height - 1], length, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &status);
  }

  double *rowVec = new double [height];
  // Send to right in the even column
  if (squareColId % 2 == 0 && squareColId != width - 1) {
    for (int i = 0; i < height; ++i) rowVec[i] = localArray[i][length - 2];
    MPI_Ssend(rowVec, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  }
  // Receive from left in the odd column
  if (squareColId % 2 == 1) {
    MPI_Recv(rowVec, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
    for (int i = 0; i < height; ++i) localArray[i][0] = rowVec[i];
  }
  // Send to right in the odd column
  if (squareColId % 2 == 1 && squareColId != width - 1) {
    for (int i = 0; i < height; ++i) rowVec[i] = localArray[i][length - 2];
    MPI_Ssend(rowVec, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  }
  // Receive from left in the even column
  if (squareColId % 2 == 0 && squareColId != 0) {
    MPI_Recv(rowVec, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
    for (int i = 0; i < height; ++i) localArray[i][0] = rowVec[i];
  }
  // Send to left in the even column
  if (squareColId % 2 == 0 && squareColId != 0) {
    for (int i = 0; i < height; ++i) rowVec[i] = localArray[i][1];
    MPI_Ssend(rowVec, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  }
  // Receive from right in the odd column
  if (squareColId % 2 == 1 && squareColId != width - 1) {
    MPI_Recv(rowVec, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    for (int i = 0; i < height; ++i) localArray[i][length - 1] = rowVec[i];
  }
  // Send to left in the odd column
  if (squareColId % 2 == 1) {
    for (int i = 0; i < height; ++i) rowVec[i] = localArray[i][1];
    MPI_Ssend(rowVec, height, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
  }
  // Receive from right in the even column
  if (squareColId % 2 == 0 && squareColId != width - 1) {
    MPI_Recv(rowVec, height, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
    for (int i = 0; i < height; ++i) localArray[i][length - 1] = rowVec[i];
  }
  delete [] rowVec;
}

void update(
  int& rank, double **localArray, double **prevArray, const int& top,
  const int& bot, const int& left, const int& right, const int& m, const int& n) {
  int topBorder = (top == 0) ? 2 : 1;
  int botBorder = (bot == m - 1) ? bot - top : bot - top + 1;
  int leftBorder = (left == 0) ? 2 : 1;
  int rightBorder = (right == n - 1) ? right - left : right - left + 1;
  for (int i = 0; i < bot - top + 3; ++i) {
    memcpy(prevArray[i], localArray[i], sizeof(double) * (right - left + 3));
  }
  for (int i = topBorder - 1; i <= botBorder + 1; ++i) {
    for (int j = leftBorder - 1; j <= rightBorder + 1; ++j) {
      double x = prevArray[i][j];
      double y = x;
      for (int k = 1; k <= 10; ++k) {
        y += sin(x * k) / pow(2.0, k);
      }
      prevArray[i][j] = y;
    }
  }
  for (int i = topBorder; i <= botBorder; ++i) {
    for (int j = leftBorder; j <= rightBorder; ++j) {
      double z = (prevArray[i-1][j] + prevArray[i][j-1] + prevArray[i][j] + prevArray[i][j+1] + prevArray[i+1][j]) / 5.0;
      localArray[i][j] = max(-100.0, min(100.0, z));
    }
  }
}


