#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N_FEATURES 2
#define N_CLASSES 2
#define N_VECTORS 32
#define N_ROWS 2
#define N_COEFFICIENTS 1
#define N_INTERCEPTS 1
#define KERNEL_TYPE 'l'
#define KERNEL_GAMMA 0.001
#define KERNEL_COEF 0.0
#define KERNEL_DEGREE 3

double vectors[32][2] = {{15435.0, 158273.0}, {14454.0, 147972.0}, {37932.0, 137228.0}, {3562.0, 159910.0}, {30004.0, 150009.0}, {400.0, 152800.0}, {14454.0, 147972.0}, {2456.0, 171214.0}, {16011.0, 155246.0}, {7337.0, 156508.0}, {7337.0, 156508.0}, {4453.0, 154939.0}, {20209.0, 139233.0}, {30202.0, 130303.0}, {35213.0, 123239.0}, {32847.0, 120141.0}, {6735.0, 173337.0}, {22963.0, 96872.0}, {31379.0, 106762.0}, {11948.0, 149090.0}, {21132.0, 157042.0}, {25016.0, 124126.0}, {20614.0, 127182.0}, {2910.0, 174296.0}, {6943.0, 182012.0}, {16558.0, 149658.0}, {6316.0, 168658.0}, {25181.0, 125750.0}, {18510.0, 145149.0}, {12598.0, 196520.0}, {12598.0, 196520.0}, {10937.0, 150976.0}};
double coefficients[1][32] = {{-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.4376177493404631, -1.0, -1.0, -1.0, -0.22343836028456915, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.08915136143851529, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5744593845832183, 0.9974453636035223, 1.0}};
double intercepts[1] = {47202.775160972305};
int weights[2] = {16, 16};

int predict_impl(double features[]) {
    int i, j, k, d, l;

    double kernels[N_VECTORS];
    double kernel;
    switch (KERNEL_TYPE) {
        case 'l':
            // <x,x'>
            for (i = 0; i < N_VECTORS; i++) {
                kernel = 0.;
                for (j = 0; j < N_FEATURES; j++) {
                    kernel += vectors[i][j] * features[j];
                }
                kernels[i] = kernel;
            }
            break;
        case 'p':
            // (y<x,x'>+r)^d
            for (i = 0; i < N_VECTORS; i++) {
                kernel = 0.;
                for (j = 0; j < N_FEATURES; j++) {
                    kernel += vectors[i][j] * features[j];
                }
                kernels[i] = pow((KERNEL_GAMMA * kernel) + KERNEL_COEF, KERNEL_DEGREE);
            }
            break;
        case 'r':
            // exp(-y|x-x'|^2)
            for (i = 0; i < N_VECTORS; i++) {
                kernel = 0.;
                for (j = 0; j < N_FEATURES; j++) {
                    kernel += pow(vectors[i][j] - features[j], 2);
                }
                kernels[i] = exp(-KERNEL_GAMMA * kernel);
            }
            break;
        case 's':
            // tanh(y<x,x'>+r)
            for (i = 0; i < N_VECTORS; i++) {
                kernel = 0.;
                for (j = 0; j < N_FEATURES; j++) {
                    kernel += vectors[i][j] * features[j];
                }
                kernels[i] = tanh((KERNEL_GAMMA * kernel) + KERNEL_COEF);
            }
            break;
    }

    int starts[N_ROWS];
    int start;
    for (i = 0; i < N_ROWS; i++) {
        if (i != 0) {
            start = 0;
            for (j = 0; j < i; j++) {
                start += weights[j];
            }
            starts[i] = start;
        } else {
            starts[0] = 0;
        }
    }

    int ends[N_ROWS];
    for (i = 0; i < N_ROWS; i++) {
        ends[i] = weights[i] + starts[i];
    }

    if (N_CLASSES == 2) {

        for (i = 0; i < N_VECTORS; i++) {
            kernels[i] = -kernels[i];
        }

        double decision = 0.;
        for (k = starts[1]; k < ends[1]; k++) {
            decision += kernels[k] * coefficients[0][k];
        }
        for (k = starts[0]; k < ends[0]; k++) {
            decision += kernels[k] * coefficients[0][k];
        }
        decision += intercepts[0];

        if (decision > 0) {
            return 0;
        }
        return 1;

    }

    double decisions[N_INTERCEPTS];
    double tmp;
    for (i = 0, d = 0, l = N_ROWS; i < l; i++) {
        for (j = i + 1; j < l; j++) {
            tmp = 0.;
            for (k = starts[j]; k < ends[j]; k++) {
                tmp += kernels[k] * coefficients[i][k];
            }
            for (k = starts[i]; k < ends[i]; k++) {
                tmp += kernels[k] * coefficients[j - 1][k];
            }
            decisions[d] = tmp + intercepts[d];
            d = d + 1;
        }
    }

    int votes[N_INTERCEPTS];
    for (i = 0, d = 0, l = N_ROWS; i < l; i++) {
        for (j = i + 1; j < l; j++) {
            votes[d] = decisions[d] > 0 ? i : j;
            d = d + 1;
        }
    }

    int amounts[N_CLASSES];
    for (i = 0, l = N_CLASSES; i < l; i++) {
        amounts[i] = 0;
    }
    for (i = 0; i < N_INTERCEPTS; i++) {
        amounts[votes[i]] += 1;
    }

    int classVal = -1;
    int classIdx = -1;
    for (i = 0; i < N_CLASSES; i++) {
        if (amounts[i] > classVal) {
            classVal = amounts[i];
            classIdx= i;
        }
    }
    return classIdx;
}

int predict(int rows, int nnz) {
  double data[2] = { (double)rows, (double)nnz };

  if(nnz < 500) {
    return CL_NATIVE_IMPL;
  } else {
    return predict_impl(data);
  }
}
