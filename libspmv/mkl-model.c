#include <math.h>

#define N_FEATURES 2
#define N_CLASSES 2
#define N_VECTORS 50
#define N_ROWS 2
#define N_COEFFICIENTS 1
#define N_INTERCEPTS 1
#define KERNEL_TYPE 'l'
#define KERNEL_GAMMA 3
#define KERNEL_COEF 0.0
#define KERNEL_DEGREE 3

double vectors[50][2] = {{0.5478737760087475, 0.48572472450974125}, {0.5208301477534012, 0.4615685085801199}, {0.9181766255986351, 0.43451169160967423}, {0.8913861154296816, 0.4098467989300647}, {-0.05600356842963666, 0.4894181572331845}, {0.8216164188115967, 0.4664760488788056}, {0.6244468844608126, 0.413098827716405}, {0.641311017432123, 0.42796551392466065}, {0.42296992653557636, 0.4679946559958921}, {0.5208301477534012, 0.4615685085801199}, {-0.20911826965276956, 0.5139350224821257}, {0.5512220315780401, 0.4653328622806164}, {0.5629626546730778, 0.47879339245938085}, {0.24159031261102007, 0.48169945091801664}, {0.24159031261102007, 0.48169945091801664}, {0.49459661036617214, 0.4917383464479667}, {0.03100951330075902, 0.46734358792643965}, {0.03594059171396891, 0.47808287888641604}, {0.6489381657132236, 0.4626439266573576}, {0.6588588574103602, 0.43971815441786255}, {0.8243252272677997, 0.4159252033027964}, {0.8875445774575025, 0.39591885414864386}, {0.46252721333633623, 0.4417772131808212}, {0.46252721333633623, 0.4417772131808212}, {0.3933078178030242, 0.45324212526549773}, {-0.3113741101155798, 0.5624367040723388}, {0.2063323568886872, 0.5183584267163684}, {0.2281836582944631, 0.5116636806700052}, {1.23902007089804, 0.42318154782055367}, {0.7114730875891878, 0.3095088618824184}, {0.8400699106917294, 0.3444021964448098}, {0.9437896683092223, 0.4543328423256594}, {0.44241410954660115, 0.4642702974424343}, {0.6772515035378656, 0.4829220653710912}, {1.0248918465603793, 0.439798063442186}, {-0.20911826965276956, 0.5139350224821257}, {0.7345083466154286, 0.44811370130442335}, {0.5294578986952608, 0.4850823433148813}, {0.7467389362380008, 0.39849305176689087}, {0.6670306186670463, 0.407223229756064}, {0.9434936901586226, 0.45414963522548707}, {-1.6684941548960344, 0.5915373276359375}, {-0.1392635105778498, 0.5203388287618371}, {0.5767975238440887, 0.4656351901272112}, {0.7765505479832859, 0.44944848606887444}, {0.17987960693802682, 0.5085360849349662}, {0.7494463776076854, 0.4031588101857246}, {0.6226929358727764, 0.45465446637734525}, {0.9020572395001406, 0.448244634147585}, {0.4060029657182407, 0.4687824655552984}};
double coefficients[1][50] = {{-10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -301.42612474031506, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, -10000.0, 6228.526556095519, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 4072.8995686446924, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0, 10000.0}};
double intercepts[1] = {13.507306706186151};
int weights[2] = {25, 25};

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

int predict(int rows, int nnz)
{
  double log_rows = log(rows);
  double log_nnz = log(nnz);

  log_rows -= 8.314063520774488;
  log_rows /= 2.428167717486321;

  log_nnz -= 10.618865565675158;
  log_nnz /= 2.785962985081255;

  double data[2] = { log_rows, log_nnz };
  return predict_impl(data);
}
