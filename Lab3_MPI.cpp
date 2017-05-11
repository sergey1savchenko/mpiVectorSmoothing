// by Serhii Savchenko

#include "header.h"

using namespace std;

const int N = 7;										// Size of input vector
const int M = 3;										// Number of iterations to do

// Create Input Vector
double* generateVector() {

	double *vector = new double[N];

	for (int i = 0; i < N; ++i) {

		double randElem = rand() % 100;					// Random value from 0 to 99
		if (randElem <= 66) {							// ~ 66% of coordinates = 0
			vector[i] = 0.0;
		}
		else {
			vector[i] = randElem;
		}
	}

	return vector;
}

// Writing to file
void writeVector(double* vector, int newSize, string fileName) {
	ofstream fout;
	fout.open(fileName, ios_base::trunc);				// Replace all file content
	for (int i = 0; i < newSize; ++i) {
		fout << "[" << i << "]" << " = " << vector[i] << "; " << endl;
	}
	fout.close();
}

int main(int argc, char**argv){

	int num_proc;
	int proc_rank;

	// Start time
	clock_t t;
	t = clock();

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	int partSize;
	if (N % num_proc == 0) {
		partSize = N / num_proc;
	}
	else {
		partSize = N / num_proc + 1;
	}

	int newSize = partSize * num_proc;
	double *inputVector = new double[newSize];
	double *vector = new double[newSize];
	double *result = new double[newSize];

	int extraCoordinates = newSize - N;

	double *tempValues = new double[newSize];

	// First initialization
	if (proc_rank == 0) {

		// Create new Input Vector
		inputVector = generateVector();
		if (extraCoordinates > 0) {
			for (int e = N; e < newSize; e++) {
				inputVector[e] = -1.0;
			}
		}

		writeVector(inputVector, newSize, "InputVector.txt");

		// Create copies
		for (int i = 0; i < newSize; ++i) {
			vector[i] = inputVector[i];
			result[i] = inputVector[i];
			tempValues[i] = tempValues[i];
		}

	}

	// Vectors
	double *temp_initialVector = new double[partSize];
	double *temp_vector = new double[partSize];
	double *temp_result = new double[partSize];

	MPI_Scatter(inputVector, partSize, MPI_DOUBLE, temp_initialVector, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(vector, partSize, MPI_DOUBLE, temp_vector, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(result, partSize, MPI_DOUBLE, temp_result, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(tempValues, newSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Smoothing
		for (int i = 0; i <= M; i++) {

			for (int j = 0; j < partSize; ++j) {
				if (temp_initialVector[j] == 0) {
					if((proc_rank*partSize + j - 1) >= 0 && (proc_rank*partSize + j + 1) < newSize){
						if (tempValues[proc_rank*partSize + j - 1] >= 0 && tempValues[proc_rank*partSize + j + 1] >= 0 &&
							tempValues[proc_rank*partSize + j - 1] <= 100 && tempValues[proc_rank*partSize + j + 1] <= 100) {
								temp_result[j] = (tempValues[proc_rank*partSize + j - 1] + tempValues[proc_rank*partSize + j + 1]) / 2;
						}
					}
				}else {
					temp_result[j] = temp_initialVector[j];
				}
				
			}

			temp_vector = temp_result;

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Gather(temp_vector, partSize, MPI_DOUBLE, vector, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			for (int b = 0; b < newSize; b++) {
				if(proc_rank*partSize + b < newSize){
					tempValues[proc_rank*partSize + b] = vector[b];
				}
			}
			MPI_Bcast(tempValues, newSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			
			if (i == 0 && proc_rank == 0) {
				writeVector(tempValues, newSize, "1stIter.txt");
			}

			MPI_Barrier(MPI_COMM_WORLD);

		}
	
		// Result Gathering
		MPI_Gather(temp_vector, partSize, MPI_DOUBLE, vector, partSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Finalize();

		if (proc_rank == 0) {
			// End time
			t = clock() - t;
			float time = ((float)t) / CLOCKS_PER_SEC;

			// Results
			cout << "Vector size: " << N << "\n" <<
				"Iterations: " << M << "\n" <<
				"Threads: " << num_proc << "\n" <<
				"Time: " << time << " seconds" << endl;

			// ResultVector
			writeVector(vector, newSize, "ResultWithExtras.txt");
			writeVector(vector, N, "ResultGeneral.txt");
		}
	

	system("Pause");
	return 0;

}