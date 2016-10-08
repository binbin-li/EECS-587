/*************************************************************************
	> File Name: test.c
	> Author:
	> Mail:
	> Created Time: Sun 02 Oct 2016 03:20:27 PM EDT
 ************************************************************************/
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

double f(const double& x);
void send(const int& bot, const int& top, const int& right, const int& left, const int& squareRowId, const int& squareColId, double **localArray, const int& width, const int& rank, const int& size);
void receive(const int& bot, const int& top, const int& right, const int& left, const int& squareRowId, const int& squareColId, double **localArray, const int& width, const int& rank, const int& size);
void update(int& rank, double **localArray, double **prevArray, const int& top, const int& bot, const int& left, const int& right, const int& m, const int& n);

int main(int argc, char** argv) {
  int rank, size;
  double startTime, endTime;
  startTime = MPI_Wtime();
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

  // Update in 10 iterations
  for (int iteration = 1; iteration <= 10; ++iteration) {
    // Send to neighbors
    send(bot, top, right, left, squareRowId, squareColId, localArray, width, rank, size);
    // Receive from neighbors
    receive(bot, top, right, left, squareRowId, squareColId, localArray, width, rank, size);
    // Update the array
    update(rank, localArray, prevArray, top, bot, left, right, m, n);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  double sum = 0;
  double squareSum = 0;
  for (int i = 1; i < bot - top + 2; ++i) {
    for (int j = 1; j < right - left + 2; ++j) {
      sum += localArray[i][j];
      squareSum += localArray[i][j] * localArray[i][j];
    }
  }

  MPI_Request request;
  if (rank != 0) {
    MPI_Isend(&sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(&squareSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
  }
  else {
    double *sums = new double [size];
    //double *squareSums = new double [size];
    sums[0] = sum;
    MPI_Status status;
    for (int i = 1; i < size; ++i) {
      MPI_Recv(sums+i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
    }
    double val = 0;
    for (int i = 0; i < size; ++i) {
      val += sums[i];
    }
    cout << val << endl;
  }

  MPI_Finalize();
  endTime = MPI_Wtime();
  if (rank == 0) cout << "Time: " << endTime - startTime << endl;
  return 0;
}

void send(const int& bot, const int& top, const int& right, const int& left, const int& squareRowId, const int& squareColId, double **localArray, const int& width, const int& rank, const int& size) {
  MPI_Request request;
  if (squareRowId == 0) {
    // send to below
    if (rank + width < size) {
      MPI_Isend(localArray[bot - top + 1], right - left + 3, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &request);
    }
  }
  else if (squareRowId == width - 1) {
    // send to upper
    if (rank - width >= 0) {
      MPI_Isend(localArray[1], right - left + 3, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &request);
    }
  }
  else {
    // send to upper and below
    MPI_Isend(localArray[bot - top + 1], right - left + 3, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &request);
    MPI_Isend(localArray[1], right - left + 3, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &request);
  }
  if (squareColId == 0) {
    // send to right
    if (rank + 1 < size) {
      double *rowVec = new double [bot - top + 3];
      for (int i = 0; i < bot - top + 3; ++i) rowVec[i] = localArray[i][right - left + 1];
      MPI_Isend(rowVec, bot - top + 3, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request);
      delete [] rowVec;
    }
  }
  else if (squareColId == width - 1) {
    // send to left
    if (rank - 1 >= 0) {
      double *rowVec = new double [bot - top + 3];
      for (int i = 0; i < bot - top + 3; ++i) rowVec[i] = localArray[i][1];
      MPI_Isend(rowVec, bot - top + 3, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
      delete [] rowVec;
    }
  }
  else {
    // send to right and left
    double *rowVec = new double [bot - top + 3];
    for (int i = 0; i < bot - top + 3; ++i) rowVec[i] = localArray[i][1];
    MPI_Isend(rowVec, bot - top + 3, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request);
    for (int i = 0; i < bot - top + 3; ++i) rowVec[i] = localArray[i][right - left + 1];
    MPI_Isend(rowVec, bot - top + 3, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &request);
    delete [] rowVec;
  }
}

void receive(const int& bot, const int& top, const int& right, const int& left, const int& squareRowId, const int& squareColId, double **localArray, const int& width, const int& rank, const int& size) {
  MPI_Status status;
  int height = bot - top + 3, length = right - left + 3;
  if (squareRowId == 0) {
    // Receive from below
    if (rank + width < size) {
      MPI_Recv(localArray[height - 1], length, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &status);
    }
  }
  else if (squareRowId == width - 1) {
    // Receive from upper
    if (rank - width >= 0) {
      MPI_Recv(localArray[0], length, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &status);
    }
  }
  else {
    // Receive from upper and below
    MPI_Recv(localArray[bot - top + 2], right - left + 3, MPI_DOUBLE, rank + width, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(localArray[0], right - left + 3, MPI_DOUBLE, rank - width, 0, MPI_COMM_WORLD, &status);
  }
  if (squareColId == 0) {
    // Receive from right
    if (rank + 1 < size) {
      double *rowVec = new double [bot - top + 3];
      MPI_Recv(rowVec, bot - top + 3, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < bot - top + 3; ++i) localArray[i][right - left + 2] = rowVec[i];
      delete [] rowVec;
    }
  }
  else if (squareColId == width - 1) {
    // Receive from left
    if (rank - 1 >= 0) {
      double *rowVec = new double [bot - top + 3];
      MPI_Recv(rowVec, bot - top + 3, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < bot - top + 3; ++i) localArray[i][0] = rowVec[i];
      delete [] rowVec;
    }
  }
  else {
    // Receive from right and left
      double *rowVec = new double [bot - top + 3];
      MPI_Recv(rowVec, bot - top + 3, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < bot - top + 3; ++i) localArray[i][right - left + 2] = rowVec[i];
      MPI_Recv(rowVec, bot - top + 3, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
      for (int i = 0; i < bot - top + 3; ++i) localArray[i][0] = rowVec[i];
      delete [] rowVec;
  }
}

void update(int& rank, double **localArray, double **prevArray, const int& top, const int& bot, const int& left, const int& right, const int& m, const int& n) {
  int topBorder = (top == 0) ? 2 : 1;
  int botBorder = (bot == m - 1) ? bot - top : bot - top + 1;
  int leftBorder = (left == 0) ? 2 : 1;
  int rightBorder = (right == n - 1) ? right - left : right - left + 1;
  for (int i = 0; i < bot - top + 3; ++i) {
    memcpy(prevArray[i], localArray[i], sizeof(double) * (right - left + 3));
  }
  for (int i = topBorder - 1; i <= botBorder + 1; ++i) {
    for (int j = leftBorder - 1; j <= rightBorder + 1; ++j) {
      prevArray[i][j] = f(prevArray[i][j]);
    }
  }
  for (int i = topBorder; i <= botBorder; ++i) {
    for (int j = leftBorder; j <= rightBorder; ++j) {
      double z = (prevArray[i-1][j] + prevArray[i+1][j] + prevArray[i][j-1] + prevArray[i][j+1] + prevArray[i][j]) / 5.0;
      localArray[i][j] = max(-100.0, min(100.0, z));
    }
  }
}

double f(const double& x) {
  double y = x;
  double divisor = x;
  int divider = 2;
  for (int i = 1; i <= 10; ++i) {
    y += sin(divisor) / divider;
    divider <<= 1;
    divisor += x;
  }
  return y;
}
