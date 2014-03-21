#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

#define EPS (.0000001)

void upperhes(int n, double *a, double *u, double *b);
void fast_left_multiply(int n, double b[][n], int i, int j, double phi, double w[][n]);
void fast_right_multiply(int n, double b[][n], int i, int j, double phi, double w[][n]);
void upperhes_slow(int n, double *aFlat, double *uFlat, double *bFlat);


// MATRIX ROUTINES: user must allocate memory for retrun param w
void matrixCopy(int n, double u[][n], double w[][n]); // w := u
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

	// perform itrations
	for(j = 0; j < n; j++){
		for(i = n-1; i > j+1; i--){
			// find x and y, then calculate phi
			x = b[i-1][j];
			y = b[i][j];
			phi = atan(-1 * (y / x));

			// kill offenders of b.
			fast_left_multiply(n, b, i, j, phi, b);
			fast_right_multiply(n, b, i, j, phi, b);

			// kill corresponding "offenders" of u.
			fast_left_multiply(n, u, i, j, phi, u);
			fast_right_multiply(n, uInv, i, j, phi, uInv);
		}
	}

	// write to outputs
	matrixFlatten(n, u, uFlat);
	matrixFlatten(n, b, bFlat);

	return;
}

void fast_left_multiply(int n, double b[][n], int i, int j, double phi, double w[][n]){
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

void fast_right_multiply(int n, double b[][n], int i, int j, double phi, double w[][n]){
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

void upperhes_slow(int n, double *aFlat, double *uFlat, double *bFlat){
	int i, j;
	double a[n][n], u[n][n], b[n][n], left[n][n], right[n][n], 
		uInv[n][n], x = 0.0, y = 0.0, phi = 0.0;

	// expand a 
	matrixExpand(n, aFlat, a);

	// init b to a
	matrixCopy(n, a, b);

	// init u, uInv to Inxn
	identity(n, u);
	identity(n, uInv);

	// perform iterations
	for(j = 0; j < n; j++){
		for(i = n-1; i > j+1; i--){
			// init left and right to I each at start of each iter
			identity(n, left);
			identity(n, right);

			// find x and y, then calculate phi
			x = b[i-1][j];
			y = b[i][j];
			phi = atan(-1 * (y / x));

			// assign rotation values to left matrix
			left[i-1][i-1] = cos(phi);
			left[i][i-1] = sin(phi);
			left[i-1][i] = -1*sin(phi);
			left[i][i] = cos(phi);

			// assign rotation values to right matrix
			right[i-1][i-1] = cos(phi);
			right[i][i-1] = -1*sin(phi);
			right[i-1][i] = sin(phi);
			right[i][i] = cos(phi);

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

/**/
int main(){
	int i, j, sum = 0, n = 4;
	double a[n][n], u[n][n], b[n][n],
		aFlat[n*n], uFlat[n*n], bFlat[n*n],
		uInv[n][n], test[n][n];

	// init a
	for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a[i][j] = rand() % 100;
			// a[i][j] = ++sum;
		}
	}

	/**/
	// flatten input, run upperhes_slow, then expand outputs
	matrixFlatten(n, a, aFlat);
	upperhes_slow(n, aFlat, uFlat, bFlat);
	matrixExpand(n, uFlat, u);
	matrixExpand(n, bFlat, b);

	// check outputs.
	printf("a:\n");
	matrixPrint(n, a);
	printf("correct u:\n");
	matrixPrint(n, u);
	printf("correct b:\n");
	matrixPrint(n, b);

	/**/
	// flatten input, run upperhes, then expand outputs
	matrixFlatten(n, a, aFlat);
	upperhes(n, aFlat, uFlat, bFlat);
	matrixExpand(n, uFlat, u);
	matrixExpand(n, bFlat, b);

	// check outputs.
	printf("a:\n");
	matrixPrint(n, a);
	printf("u:\n");
	matrixPrint(n, u);
	printf("b:\n");
	matrixPrint(n, b);

	// test code
	/**/
	printf("test code:\n\n");
	printf("u uT (should be identity, confirms orthogonality)\n");
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
				printf("a[%d][%d]: %f =/= test[%d][%d]: %f\n", i, j, a[i][j], 
					i, j, test[i][j]);
				i = INT_MAX; // use i as inequality flag
				break;
			}
		}
		if (i > n){
			break;
		}
	}
	if(i < INT_MAX) // test flag.
		printf("YES\n");
	/**/

	return 1;
}
/**/