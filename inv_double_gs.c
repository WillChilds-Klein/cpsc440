#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>


void inv_double_gs(double* a, int n, double* u, double* b);

// VECTOR ROUTINES: user must allocate memory for return param w
void vectorSubtract(int n, double* u, double* v, double* w); 	// u - v = w
void vectorScale(int n, double* u, double c, double* w);		// cu = w
void dotProduct(int n, double* u, double* v, double* w);		// u * v = w
void vectorProject(int n, double* u, double* v, double* x, double* w);// ((u*v) / ||v||^2)x = w
void vectorNormalize(int n, double* u, double* v, double* w);	// u / ||v|| = w
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


void inv_double_gs(double* a, int n, double* u, double* b){
	double A[n][n], U[n][n], G[n][n], B[n][n], Ut[n][n], proj[n], gproj[n],
		aVecs[n][n], uVecs[n][n], gVecs[n][n], wVecs[n][n];
	int i, j;

	matrixExpand(n, a, A);

	// init G to Identity
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == j)
				G[i][j] = 1.0;
			else
				G[i][j] = 0.0;
		}
	}

	// init U to A
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			U[i][j] = A[i][j];
		}
	}

	// aVecs = aT b/c we want to work with the column vectors. 
	// re-transp both this and uVecs post-GS to get correct matrices.
	matrixTranspose(n, A, aVecs);
	matrixTranspose(n, U, uVecs);
	matrixTranspose(n, G, gVecs);

	for(i = 0; i < n; i++){
		// single GS
		for(j = i+1; j < n; j++){
			// store component to take off in temp
			vectorProject(n, aVecs[j], uVecs[i], uVecs[i], proj);
			vectorProject(n, aVecs[j], uVecs[i], gVecs[i], gproj);

			// subtract off appropriate component
			vectorSubtract(n, uVecs[j], proj, uVecs[j]);
			vectorSubtract(n, gVecs[j], gproj, gVecs[j]);
		}
		vectorNormalize(n, gVecs[i], uVecs[i], gVecs[i]);
		vectorNormalize(n, uVecs[i], uVecs[i], uVecs[i]);
	}

	// for(i = n-1; i >= 0; i--){
	// 	// double GS
	// 	for(j = i-1; j >= 0; j--){
	// 		// store component to take off in temp
	// 		vectorProject(n, uVecs[i], uVecs[j], proj);

	// 		// subtract off appropriate component
	// 		vectorSubtract(n, uVecs[j], proj, uVecs[j]);
	// 	}
	// 	vectorNormalize(n, aVecs[i], aVecs[i] uVecs[i]);
	// }

	// re-transp to get correct A, G, and U.
	matrixTranspose(n, aVecs, A);
	matrixTranspose(n, uVecs, U);
	matrixTranspose(n, gVecs, G);

	printf("A:\n");
	matrixPrint(n, A);
	printf("U:\n");
	matrixPrint(n, U);
	printf("G:\n");
	matrixPrint(n, G);

	// we need U's transpose
	matrixTranspose(n, U, Ut);

	// Now, find the inverse
	matrixMultiply(n, G, Ut, B);

	printf("A:\n");
	matrixPrint(n, A);
	printf("B:\n");
	matrixPrint(n, B);

	// check that AB is Identity.
	matrixMultiply(n, A, B, G);
	printf("AB: (should be Identity)\n");
	matrixPrint(n, G);

}

void vectorSubtract(int n, double* u, double* v, double* w){
	int i;

	for(i = 0; i < n; i++){
		w[i] = u[i] - v[i];
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

void vectorProject(int n, double* u, double* v, double* x, double* w){
	double product[1], magnitude[1], c;

	dotProduct(n, u, v, product);
	vectorMagnitude(n, v, magnitude);
	c = *product / pow(*magnitude, 2);
	vectorScale(n, x, c, w);

	return;
}

void vectorNormalize(int n, double* u, double* v, double* w){
	int i;
	double magnitude[1];

	vectorMagnitude(n, v, magnitude);

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

int main(){
	int i, j, sum = 1, n = 2;
	double a[n][n], u[n][n], b[n][n], x[n], y[1], 
		flatA[n*n], flatU[n*n], flatB[n*n];

	for(i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			a[i][j] = sum++;
		}
	}

	// Test subroutines.
	/** /
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
	vectorProject(n, a[1], a[0], x);
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
	/**/

	matrixFlatten(n, a, flatA);
	inv_double_gs(flatA, n, flatU, flatB);

}