#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gmath.h"
#include "matrix.h"
#include "ml6.h"


//vector functions
//normalize vetor, should modify the parameter
void normalize( double *vector ) {
}

//Return the dot porduct of a . b
double dot_product( double *a, double *b ) {
    int out = 0;
    for(int i = 0; i < 3; i++){
        out += a[i] * b[i];
    }
	return out;
}


//Calculate the surface normal for the triangle whose first
//point is located at index i in polygons
double * calculate_normal(struct matrix *polygons, int i) {
    double x0 = polygons->m[0][i];
    double y0 = polygons->m[1][i];
    double z0 = polygons->m[2][i];

    double x1 = polygons->m[0][i+1];
    double y1 = polygons->m[1][i+1];
    double z1 = polygons->m[2][i+1];

    double x2 = polygons->m[0][i+2];
    double y2 = polygons->m[1][i+2];
    double z2 = polygons->m[2][i+2];

    double A[] = {x0 - x1, y0 - y1, z0 - z1};
    double B[] = {x2 - x1, y2 - y1, z2 - z1};

    double * N = malloc(sizeof(double) * 3);
    N[0] = A[1]*B[2] - A[2]*B[1];
    N[1] = A[2]*B[0] - A[0]*B[2];
    N[2] = A[0]*B[1] - A[1]*B[0];

	return N;
}
