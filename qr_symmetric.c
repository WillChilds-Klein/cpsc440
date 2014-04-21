#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#define EPS (.001)

void qr_symmetric(double *a, int n, double *b);
void apply_rotations(int n, double b[][n], int i, double w[][n]);
void shift(int n, double b[][n], int i, double mu, double bShifted[][n]);
void unshift(int n, double bShifted[][n], int i, double mu, double b[][n]);
int converged(int n, double b[][n], int i);

// UPPERHES ROUTINES: routines used for finding upper hessenberg form
void upperhes(int n, double *a, double *u, double *b);
void fast_left_multiply(int n, double b[][n], int i, double phi, double w[][n]);
void fast_right_multiply(int n, double b[][n], int i, double phi, double w[][n]);

// MATRIX ROUTINES: user must allocate memory for retrun param w
void matrixCopy(int n, double u[][n], double w[][n]); 						// w := u
void matrixTranspose(int n, double u[][n], double w[][n]);					// uT = w
void matrixMultiply(int n, double u[][n] , double v[][n], double w[][n]); 	// uv = w
void identity(int n, double u[][n]); 					 // inits u to n x n identity

// MATRIX UTILITIES: user must allocate memory for return param w
void matrixExpand(int n, double *aFlat, double w[][n]);
void matrixFlatten(int n, double aExpanded[][n], double *w);

// PRINTING UTILITIES
void vectorPrint(int n, double *u);
void matrixPrint(int n, double u[][n]);
void flatPrint(int n, double *u);

void qr_symmetric(double *aFlat, int n, double *bFlat){
	int i, j, nCurr;
	double a[n][n], b[n][n], uFlat[n*n], bUpperHesFlat[n*n], bShifted[n][n], 
		spectrum[n], mu, x, y, phi;

	// make symmetric a tridiagonal/upperhes
	upperhes(n, aFlat, uFlat, bUpperHesFlat);
	matrixExpand(n, bUpperHesFlat, b);

	// init nCurr to n
	nCurr = n;

	for(i = 0; i < n; i++){
		j=0;
		printf("i = %d\n", i);
		while(j < 7){//!converged(n, b, i)){

			printf("i = %d, j = %d, b:\n", i, j);
			matrixPrint(n, b);

			// perform shift
			mu = b[i][i];
			shift(n, b, i, mu, bShifted);

			// apply rotations
			apply_rotations(n, bShifted, i, bShifted);

			// add shift back to b
			unshift(n, bShifted, i, mu, b);
			j++;
		}
	}

	matrixFlatten(n, b, bFlat);

	return;
}

// no mem write overlap
void apply_rotations(int n, double bOriginal[][n], int i, double w[][n]){
	int j, k;
	double b[n][n], temp1[2][n], temp2[n][2], phiArr[n], x, y, phi; // might need array of phi's...
	double Qj[n][n], QjT[n][n];

	matrixCopy(n, bOriginal, b);

	// make b lower triangular
	for(j = n-2; j >= i; j--){ // j is as used in Step 3 of spec
		// find x and y (zero'ing x out), then calculate phi
		x = b[j][j+1];
		y = b[j+1][j+1];
		phi = atan(x / y);
		phiArr[j] = phi;
		
		// create Qj
		identity(n, Qj);
		Qj[j][j] = cos(phi);
		Qj[j][j+1] = sin(phi);
		Qj[j+1][j] = -1 * sin(phi);
		Qj[j+1][j+1] = cos(phi);

		printf("Q%d\n", j);
		matrixPrint(n, Qj);

		// simulate left multiplication
		for(k = i; k < n; k++){
			temp1[0][k] = cos(phi)*b[j][k] - sin(phi)*b[j+1][k];
			temp1[1][k] = sin(phi)*b[j][k] + cos(phi)*b[j+1][k];
		}

		// update b
		for(k = i; k < n; k++){
			b[j][k] = temp1[0][k];
			b[j+1][k] = temp1[1][k];
		}
	}

	// printf("lower triangular?\n");
	// matrixPrint(n, b);

	// "apply adjoints"
	for(j = n-2; j >= i; j--){ // j is as used in Step 3 of spec
		// retrieve appropriate phi
		phi = phiArr[j];

		// create QjT
		identity(n, QjT);
		QjT[j][j] = cos(phi);
		QjT[j][j+1] = -1 * sin(phi);
		QjT[j+1][j] = sin(phi);
		QjT[j+1][j+1] = cos(phi);

		printf("Q%dT\n", j);
		matrixPrint(n, QjT);
		
		// simulate left multiplication
		for(k = i; k < n; k++){
			temp2[k][0] = cos(phi)*b[k][j] - sin(phi)*b[k][j+1];
			temp2[k][1] = sin(phi)*b[k][j] + cos(phi)*b[k][j+1];
		}

		// update b
		for(k = i; k < n; k++){
			b[k][j] = temp2[k][0];
			b[k][j+1] = temp2[k][1];
		}
	}

	// printf("tridiagonal?\n");
	// matrixPrint(n, b);

	for(k = 0; k < n; k++){
		memcpy(w[i], b[i], n*sizeof(double));
	}

	return;
}

void shift(int n, double b[][n], int i, double mu, double bShifted[][n]){
	int k;

	matrixCopy(n, b, bShifted);

	for(k = i; k < n; k++){
		bShifted[k][k] = b[k][k] - mu;
	}

	return;
}

void unshift(int n, double bShifted[][n], int i, double mu, double b[][n]){
	int k;

	matrixCopy(n, bShifted, b);

	for(k = i; k < n; k++){
		b[k][k] = bShifted[k][k] + mu;
	}

	return;
}

int converged(int n, double b[][n], int i){
	return (abs(b[i][i+1]) <= EPS);
}

void upperhes(int n, double *aFlat, double *uFlat, double *bFlat){
	int i, j, k;
	double a[n][n], u[n][n], b[n][n], left[n][n], right[n][n], 
		uInv[n][n], x = 0.0, y = 0.0, phi = 0.0;

	// expand a 
	matrixExpand(n, aFlat, a);

	// init b to a
	matrixCopy(n, a, b);

	// init u and uInv to I_n
	identity(n, u);
	identity(n, uInv);

	// perform itrations
	for(j = 0; j < n; j++){
		for(i = n-1; i > j+1; i--){
			// find x and y, then calculate phi
			x = b[i-1][j];
			y = b[i][j];
			phi = atan(-1 * (y / x));

			// kill offenders of b.
			fast_left_multiply(n, b, i, phi, b);
			fast_right_multiply(n, b, i, phi, b);

			// kill corresponding "offenders" of u.
			fast_left_multiply(n, u, i, phi, u);
			fast_right_multiply(n, uInv, i, phi, uInv);
		}
	}

	// write to outputs
	matrixFlatten(n, u, uFlat);
	matrixFlatten(n, b, bFlat);

	return;
}

void fast_left_multiply(int n, double b[][n], int i, double phi, double w[][n]){
	int k;
	double temp[2][n];

	// simulate left multiplication
	for(k = 0; k < n; k++){
		temp[0][k] = cos(phi)*b[i-1][k] - sin(phi)*b[i][k];
		temp[1][k] = sin(phi)*b[i-1][k] + cos(phi)*b[i][k];
	}

	memcpy(w[i-1], temp[0], n*sizeof(double));
	memcpy(w[i], temp[1], n*sizeof(double));

	return;
}

void fast_right_multiply(int n, double b[][n], int i, double phi, double w[][n]){
	int k;
	double temp[n][2];

	// simulate right multiplication
	for(k = 0; k < n; k++){
		temp[k][0] = cos(phi)*b[k][i-1] - sin(phi)*b[k][i];
		temp[k][1] = sin(phi)*b[k][i-1] + cos(phi)*b[k][i];
	}

	for(k = 0; k < n; k++){
		w[k][i-1] = temp[k][0];
		w[k][i] = temp[k][1];
	}

	return;
}

void matrixCopy(int n, double u[][n], double w[][n]){
	int i;

	for(i = 0; i < n; i++){
		memcpy(w[i], u[i], n*sizeof(double));
	}

	return;
}

// fixed for memory overlap issue btwn u and w
void matrixTranspose(int n, double u[][n], double w[][n]){
	int i, j;
	double temp[n][n];

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			temp[i][j] = u[j][i];
		}
	}

	// copy over transpose
	for(i = 0; i < n; i++){
		memcpy(w[i], temp[i], n*sizeof(double));
	}

	return;
}

// no mem overlap btwn u, v, and w
void matrixMultiply(int n, double u[][n], double v[][n], double w[][n]){
	int i, j, k;
	double vectorSum, temp[n][n];

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			vectorSum = 0.0;
			for(k = 0; k < n; k++){
				vectorSum += u[i][k] * v[k][j];
			}
			temp[i][j] = vectorSum;
		}
	}

	// copy over temp
	for(i = 0; i < n; i++){
		memcpy(w[i], temp[i], n*sizeof(double));
	}

	return;
}

void identity(int n, double u[][n]){
	int i, j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == j)
				u[i][j] = 1.0;
			else
				u[i][j] = 0.0;
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

int main(){
	int i, j, sum = 0, n = 4;
	double a[n][n], u[n][n], b[n][n],
		aFlat[n*n], uFlat[n*n], bFlat[n*n],
		uInv[n][n], bShifted[n][n], temp;

	// symmetric matrix should create tridiagonal
	for (i = 0; i < n; i++){
		for (j = 0; j < i; j++){
			temp = (rand() % 10) + 1;
			a[i][j] = temp;
			a[j][i] = temp;
		}
		a[i][i] = i+1;
	}

	matrixFlatten(n, a, aFlat);
	upperhes(n, aFlat, uFlat, bFlat);
	matrixExpand(n, uFlat, u);
	matrixExpand(n, bFlat, b);

	// check outputs.
	printf("symmetric a:\n");
	matrixPrint(n, a);
	printf("u:\n");
	matrixPrint(n, u);
	printf("tridiagonal b:\n");
	matrixPrint(n, b);

	matrixFlatten(n, b, bFlat);
	qr_symmetric(bFlat, n, aFlat);
	matrixExpand(n, aFlat, a);

	printf("new a is tridiagonal b above.\n");
	printf("spectrum of a:\n");
	matrixPrint(n, a);

	return 1;
}