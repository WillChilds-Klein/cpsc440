#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#define EPS (.0000001)

void upperhes(int n, double *a, double *u, double *b);
void kill_ze_offenders(int n, double b[][n], int i, int j, double phi);
void upperhes_slow(int n, double *aFlat, double *uFlat, double *bFlat);
void assignVals(int n, double left[][n], double right[][n], double phi, int i, int j);


// MATRIX ROUTINES: user must allocate memory for retrun param w
void matrixCopy(int n, double u[][n], double w[][n]); // w = u
void matrixTranspose(int n, double u[][n], double w[][n]);		// uT = w
void matrixMultiply(int n, double u[][n] , double v[][n], double w[][n]); // uv = w
void identity(int n, double u[][n]); // inits u to n x n identity

// MATRIX UTILITIES: user must allocate memory for return param w
void matrixExpand(int n, double *aFlat, double w[][n]);
void matrixFlatten(int n, double aExpanded[][n], double *w);

// PRINTING UTILITIES
void vectorPrint(int n, double *u);
void matrixPrint(int n, double u[][n]);
void flatPrint(int n, double *u);

// v0 - compose sparse u1,u2,...,un and actually multiply them out
// v1 - all of this multiplication ends up being too slow, so instead
// of actually composing the full sparse matrices, just apply the pertinent
// trig functions w vals, and do the same to the identity to iteratively 
// compose u.

void upperhes(int n, double *aFlat, double *uFlat, double *bFlat){
	int i, j, k;
	double a[n][n], u[n][n], b[n][n], left[n][n], right[n][n], 
		uInv[n][n], x = 0.0, y = 0.0, phi = 0.0;

	// expand a 
	matrixExpand(n, aFlat, a);

	// init b to a
	matrixCopy(n, a, b);

	// init u and uInv to Inxn
	identity(n, u);
	identity(n, uInv);

	// perform ze itrations
	for(j = 0; j < n; j++){
		for(i = n-1; i > j+1; i--){
			// init left and right to I each at start of each iter
			identity(n, left);
			identity(n, right);

			// find x and y, then calculate phi
			x = b[i-1][j];
			y = b[i][j];
			phi = atan(-1 * (y / x));

			// kill offenders of b.
			kill_ze_offenders(n, b, i, j, phi);

			// symmetrically kill "offenders" of u.
			kill_ze_offenders(n, u, i, j, phi);
		}
	}

	// write to outputs
	matrixFlatten(n, u, uFlat);
	matrixFlatten(n, b, bFlat);

	return;
}

void kill_ze_offenders(int n, double b[][n], int i, int j, double phi){

}

void upperhes_slow(int n, double *aFlat, double *uFlat, double *bFlat){
	int i, j;
	double a[n][n], u[n][n], b[n][n], left[n][n], right[n][n], 
		uInv[n][n], x = 0.0, y = 0.0, phi = 0.0;

	// expand a 
	matrixExpand(n, aFlat, a);

	// init b to a
	matrixCopy(n, a, b);

	// init u and uInv to Inxn
	identity(n, u);
	identity(n, uInv);

	// perform ze itrations
	for(j = 0; j < n; j++){
		for(i = n-1; i > j+1; i--){
			// init left and right to I each at start of each iter
			identity(n, left);
			identity(n, right);

			// find x and y, then calculate phi
			x = b[i-1][j];
			y = b[i][j];
			phi = atan(-1 * (y / x));

			// assign values to sparse matrices appropriately
			assignVals(n, left, right, phi, i, j);

			// keel ze y
			matrixMultiply(n, left, b, b);
			matrixMultiply(n, b, right, b);

			// compose u, uInv
			matrixMultiply(n, left, u, u);
			matrixMultiply(n, uInv, right, uInv);
		}
	}

	// write to outputs
	matrixFlatten(n, u, uFlat);
	matrixFlatten(n, b, bFlat);

	return;
}

void assignVals(int n, double left[][n], double right[][n], double phi, int i, int j){
	left[i-1][i-1] = cos(phi);
	left[i-1][i] = -1*sin(phi);
	left[i][i-1] = sin(phi);
	left[i][i] = cos(phi);

	right[i-1][i-1] = cos(phi);
	right[i-1][i] = sin(phi);
	right[i][i-1] = -1*sin(phi);
	right[i][i] = cos(phi);

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
// fixed for memory overlap issue btwn u, v, and w
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
	int i, j, sum = 0, n = 8;
	double a[n][n], u[n][n], b[n][n],
		aFlat[n*n], uFlat[n*n], bFlat[n*n],
		uInv[n][n], test[n][n];

	// init a
	for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a[i][j] = rand() % 100;//++sum;
		}
	}

	// flatten input, run upperhes, then expand outputs
	matrixFlatten(n, a, aFlat);
	upperhes_slow(n, aFlat, uFlat, bFlat);
	matrixExpand(n, uFlat, u);
	matrixExpand(n, bFlat, b);

	// check outputs.
	printf("a:\n");
	matrixPrint(n, a);
	printf("u:\n");
	matrixPrint(n, u);
	printf("b:\n");
	matrixPrint(n, b);
	printf("test code:\n\n");
	printf("u uT (should be identity)\n");
	matrixTranspose(n, u, uInv);
	matrixMultiply(n, u, uInv, test);
	matrixPrint(n, test);
	printf("uT b u (should be original a)\n");
	matrixMultiply(n, uInv, b, test);
	matrixMultiply(n, test, u, test);
	matrixPrint(n, test);
	printf("original a:\n");
	matrixPrint(n, a);

	printf("does it work?\n");
	for(i = 0; i < n; ++i){
		for(j = 0; j < n; j++){
			// approximate equals
			if(fabs(test[i][j] - a[i][j]) > EPS){
				printf("NO\n");
				printf("a[%d][%d]: %f, test[%d][%d]: %f\n", i, j, a[i][j], 
					i, j, test[i][j]);
				i = INT_MAX;
				break;
			}
		}
		if (i > n){
			break;
		}
	}
	if(i < INT_MAX)
		printf("YES\n");

	return 1;
}