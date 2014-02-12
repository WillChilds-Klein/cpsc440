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
void vectorNormalize(int n, double* u, double* w);				// ||u|| = w
void vectorMagnitude(int n, double* u, double* w);				// |u| = w

// MATRIX ROUTINES: user must allocate memory for retrun param w
void matrixTranspose(int n, double u[][n], double w[][n]);		// uT = w
void matrixMultiply(int n, double u[][n] , double v[][n], double w[][n]); // uv = w

// MATRIX UTILITIES: user must allocate memory for return param w
void matrixExpand(int n, double* aFlat, double w[][n]);
void matrixFlatten(int n, double aExpanded[][n], double* w);

// PRINTING UTILITIES
void vectorPrint(int n, double* u);
void matrixPrint(int n, double u[][n]);
void flatPrint(int n, double* u);


int main(){
	int i, j, sum = 0, n = 3;
	double a[3][3], u[n][n], b[3][3], x[n], y[1];

	for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a[i][j] = sum++;
		}
	}

	// Test subroutines.

	printf("u:\n");
	vectorPrint(n, a[1]);	
	printf("v:\n");
	vectorPrint(n, a[0]);
	printf("u - v:\n");
	vectorSubtract(n, a[1], a[0], x);
	vectorPrint(n, x);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);	
	printf("v:\n");
	vectorPrint(n, a[2]);
	printf("u + v:\n");
	vectorAdd(n, a[1], a[2], x);
	vectorPrint(n, x);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);
	printf("9u:\n");
	vectorScale(n, a[1], 9, x);
	vectorPrint(n, x);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);	
	printf("v:\n");
	vectorPrint(n, a[0]);
	printf("u * v:\n");
	dotProduct(n, a[1], a[0], y);
	printf("%f\n", *y);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);	
	printf("v:\n");
	vectorPrint(n, a[0]);
	printf("comp(u,v:)\n");
	vectorComp(n, a[1], a[0], x);
	vectorPrint(n, x);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);
	printf("||u||:\n");
	vectorNormalize(n, a[1], x);
	vectorPrint(n, x);
	printf("\n\n");

	printf("u:\n");
	vectorPrint(n, a[1]);
	printf("|u|:\n");
	vectorMagnitude(n, a[1], y);
	printf("%f\n", *y);
	printf("\n\n");

	printf("a:\n");
	matrixPrint(n, a);
	printf("aT:\n");
	matrixTranspose(n, a, b);
	matrixPrint(n, b);
	printf("\n\n");

	printf("a:\n");
	matrixPrint(n, a);
	printf("a^2:\n");
	matrixMultiply(n, a, a, b);
	matrixPrint(n, b);
	printf("\n\n");

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
	printf("mag: %f\n", *magnitude);

	for(i = 0; i < n; i++){
		w[i] = u[i] / (*magnitude);
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

void matrixTranspose(int n, double u[][n], double w[][n]){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[i][j] = u[j][i];
		}
	}

	return;
}

void matrixMultiply(int n, double u[][n], double v[][n], double w[][n]){
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

void matrixExpand(int n, double* aFlat, double w[][n]){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[i][j] = aFlat[(i*n)+j];
		}
	}

	return;
}

void matrixFlatten(int n, double aExpanded[][n], double* w){
	int i, j;

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			w[(i*n)+j] = aExpanded[i][j];
		}
	}

	return;
}

void vectorPrint(int n, double* u){
	int i;
	printf("[");
	for(i = 0; i < n; i++){
		if(i == n-1)
			printf(" %f", u[i]);
		else
			printf(" %f,", u[i]);
	}
	printf("]\n");

	return;
}

void matrixPrint(int n, double u[][n]){
	int i, j;

	for(i = 0; i < n; i++){
		printf("[");
		for(j = 0; j < n; j++){
			if(j == n-1)
				printf(" %f]\n", u[i][j]);
			else
				printf(" %f,", u[i][j]);
		}
	}
	printf("\n");

	return;
}

void flatPrint(int n, double* u){
	int i;

	for(i = 0; i < n*n; i++){
		printf("%f, ", u[i]);
	}
	printf("\n");

	return;
}