#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>


void inv_double_gs(double* a, int n, double* u, double* b);

// VECTOR ROUTINES: user must allocate memory for return param w
void vectorSubtract(int n, double* u, double* v, double* w); 	// u - v = w
void vectorAdd(int n, double* u, double* v, double* w); 		// u + v = w
void vectorScale(int n, double* u, double c, double* w);		// cu = w
void dotProduct(int n, double* u, double* v, double* w);		// u * v = w
void vectorComp(int n, double* u, double* v, double* w);		// comp(u,v) = w
void vectorNorm(int n, double* u, double* w);					// ||u|| = w
void vectorMagnitude(int n, double* u, double* w);				// |u| = w

// MATRIX ROUTINES: user must allocate memory for retrun param w
void matrixTranspose(int n, double** u, double** w);			// uT = w
void matrixMultiply(int n, double** u, double** v, double** w); // uv = w

// MATRIX UTILITIES: user must allocate memory for return param w
void matrixExpand(int n, double* aFlat, double** w);
void matrixFlatten(int n, double** aExpanded, double* w);


int main(){

}

void inv_double_gs(double* a, int n, double* u, double* b){

}

void vectorSubtract(int n, double* u, double* v, double* w){
	int i;

	for(i = 0; i < n; i++){
		w[i] = u[i] - v[i];
	}

	return;
}

void vectorAdd(int n, double* u, double* v, double* w){
	int i;

	for(i = 0; i < n; i++){
		w[i] = u[i] + v[i];
	}

	return;
}

void vectorScale(int n, double* u, double c, double* w){
	int i;

	for(i = 0; i < n; i++){
		w[i] = c * u[i];
	}

	return;
}

void dotProduct(int n, double* u, double* v, double* w){
	int i;
	double product;

	product = 0;

	for(i = 0; i < n; i++){
		product += u[i] * v[i];
	}

	*w = product;
	return;
}

void vectorComp(int n, double* u, double* v, double* w){
	double product[1];

	dotProduct(n, u, v, product);

	vectorScale(n, u, *product, w);

	return;
}

void vectorNormalize(int n, double* u, double* w){
	int i;
	double magnitude[1];

	// potential segfault. might have to initialize.
	vectorMagnitude(n, u, magnitude);

	for(int i; i < n; i++){
		w[i] /= *magnitude;
	}

	return;
}

void vectorMagnitude(int n, double* u, double* w){	
	int i;

	double sum = 0;
	for(i = 0; i < n; i++){
		sum += pow(u[i], 2);
	}

	*w = sqrt(sum);
	return;
}

void matrixTranspose(int n, double** u, double** w){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[i][j] = u[j][i];
		}
	}

	return;
}

void matrixMultiply(int n, double** u, double** v, double** w){
	int i, j, k;
	double vectorSum;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			vectorSum = 0;
			for(k = 0; k < n; k++){
				vectorSum += u[i][k] * v[k][j];
			}
			w[i][j] = vectorSum;
		}
	}

	return;
}

void matrixExpand(int n, double* aFlat, double** w){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[i][j] = aFlat[(i*n)+j];
		}
	}

	return;
}

void matrixFlatten(int n, double** aExpanded, double* w){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[(i*n)+j] = aExpanded[i][j];
		}
	}

	return;
}