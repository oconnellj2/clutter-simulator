#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

// Calculates the center of a facet.
void center(double flist[][3], double *value) {
	value[0] = (flist[0][0] + flist[1][0] + flist[2][0]) / 3;
	value[1] = (flist[0][1] + flist[1][1] + flist[2][1]) / 3;
	value[2] = (flist[0][2] + flist[1][2] + flist[2][2]) / 3;
}

// Vector in a 3D Cartesian coordinate system.
void vec(double *l0, double *l1, double *value) {
	if (l1 == NULL) {
		value[0] = l0[0];
		value[1] = l0[1];
		value[2] = l0[2];
	} else {
		value[0] = l1[0] - l0[0];
		value[1] = l1[1] - l0[1];
		value[2] = l1[2] - l0[2];
	}
}

// Calculates the magnitude of a vector point(value).
double vlen(double *value) {
	return sqrt(pow(value[0], 2) + pow(value[1], 2) + pow(value[2], 2));
}

// Calculates the cross product of 2 vector points.
void xprd(double *vec1, double *vec2, double *value) {
	value[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	value[1] = -(vec1[0] * vec2[2] - vec1[2] * vec2[0]);
	value[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

// Calculates the normal value of a facet.
void norm(double flist[][3], double *value) {
	double v1[3];
	double v2[3];

	vec(flist[1], flist[0], v1);
	vec(flist[1], flist[2], v2);
	xprd(v1, v2, value);
}

/*********************************************************
 * A function called from python to calculate the twtt.
 *
 * Args:
 *      flist(3d array): Facet objects.
 *      xcvr(Array): Spacecrafts coordinate vector.
 *      twtt(Array): Empty.
 *      c(Float): Speed of light
 *      flen(int): Len of flist.
 *
 * Returns:
 *      twtt(Array): Two way travel time values to python.
 **********************************************************/
double *c_twtt(double flist[][3][3], double *xcvr, double *twtt, double c, int flen) {
	double cen[3];
	double vsc[3];
	for (int i = 0; i < flen; i++) {
		center(flist[i], cen);
		vec(cen, xcvr, vsc);
		double vsclen = vlen(vsc);
		twtt[i] = (vsclen)*2 / c;
	}
	return twtt;
}

/*********************************************************
 * A function called from python to calculate the pwr.
 *
 * Args:
 *      flist(3d array): Facet objects.
 *      xcvr(Array): Spacecrafts coordinate vector.
 *      pwr(Array): Empty.
 *      flen(int): Len of flist.
 *
 * Returns:
 *      pwr(Array): Power values to python.
 *********************************************************/
double *c_pwr(double flist[][3][3], double *xcvr, double *pwr, int flen) {
	double cen[3];
	double vsc[3];
	double fNorm[3];

	for (int i = 0; i < flen; i++) {
		center(flist[i], cen);
		vec(cen, xcvr, vsc);
		double vsclen = vlen(vsc);
		norm(flist[i], fNorm);
		double dprd = ((vsc[0] * fNorm[0]) + (vsc[1] * fNorm[1]) + (vsc[2] * fNorm[2]));
		double ct = dprd / (vsclen * vlen(fNorm));
		pwr[i] = fabs(pow((0.5 * vlen(fNorm) * ct), 2) / pow(vsclen, 4));
	}
	return pwr;
}