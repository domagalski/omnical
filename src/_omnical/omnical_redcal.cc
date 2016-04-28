//Jeff 2012-07-19
#include <stdint.h>
#include <stdio.h>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <memory.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <numeric>
#include <algorithm>
#include <functional>
#include <numeric>
#include "include/omnical_redcal.h"
#include <algorithm>
#define uint unsigned int
using namespace std;
extern "C" {
    #include <cblas.h>
}
#include "include/complexCL.h"

const string FILENAME = "omnical_redcal.cc";
const float PI = atan2(0, -1);
const float UBLPRECISION = pow(10, -3);
const float MIN_NONE_ZERO = pow(10, -10);
const float MAX_NONE_INF = pow(10, 10);
const float MAX_POW_2 = pow(10, 10); //limiting the base of ^2 to be 10^10

void initcalmodule(calmemmodule* module, redundantinfo* info){
	int nant = info->nAntenna;
	int nbl = info->bl2d.size();
    int nubl = info->ublindex.size();
	int ncross = nbl;
	(module->amp1).resize(ncross);
	(module->amp2).resize(ncross);
	(module->amp3).resize(ncross);
	(module->pha1).resize(ncross);
	(module->pha2).resize(ncross);
	(module->pha3).resize(ncross);
	(module->x1).resize(nubl + nant);
	(module->x2).resize(nubl + nant);
	(module->x3).resize(nubl + nant);
	(module->x4).resize(nubl + nant);

	(module->g1).resize(nant);
	(module->g0) = (module->g1);
	(module->g2) = (module->g1);
	(module->g3) = (module->g1);

	(module->adata1).resize(nbl);
	(module->adata2) = (module->adata1);
	(module->adata3) = (module->adata1);

	(module->cdata1).resize(ncross);
	(module->cdata2) = (module->cdata1);
	(module->cdata3) = (module->cdata1);


	(module->ubl1).resize(nubl);
	(module->ubl0) = (module->ubl1);
	(module->ubl2) = (module->ubl1);
	(module->ubl3) = (module->ubl1);

	(module->ublgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ublgrp1)[i].resize(info->ublcount[i]);
	}
	(module->ublgrp2) = (module->ublgrp1);

	(module->ubl2dgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ubl2dgrp1)[i].resize(info->ublcount[i]);
	}
	(module->ubl2dgrp2) = (module->ubl2dgrp1);
	return;
}

// Overloading the function...
void initcalmodule(calmemmodule_clh* module, redundantinfo* info){
	int nant = info->nAntenna;
	int nbl = info->bl2d.size();
    int nubl = info->ublindex.size();
	int ncross = nbl;
	(module->amp1).resize(ncross);
	(module->amp2).resize(ncross);
	(module->amp3).resize(ncross);
	(module->pha1).resize(ncross);
	(module->pha2).resize(ncross);
	(module->pha3).resize(ncross);
	(module->x1).resize(nubl + nant);
	(module->x2).resize(nubl + nant);
	(module->x3).resize(nubl + nant);
	(module->x4).resize(nubl + nant);

	(module->g1).resize(nant);
	(module->g0) = (module->g1);
	(module->g2) = (module->g1);
	(module->g3) = (module->g1);

	(module->adata1).resize(nbl);
	(module->adata2) = (module->adata1);
	(module->adata3) = (module->adata1);

	(module->cdata1).resize(ncross);
	(module->cdata2) = (module->cdata1);
	(module->cdata3) = (module->cdata1);


	(module->ubl1).resize(nubl);
	(module->ubl0) = (module->ubl1);
	(module->ubl2) = (module->ubl1);
	(module->ubl3) = (module->ubl1);

	(module->ublgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ublgrp1)[i].resize(info->ublcount[i]);
	}
	(module->ublgrp2) = (module->ublgrp1);

	(module->ubl2dgrp1).resize(nubl);
	for (int i = 0; i < nubl; i++){
		(module->ubl2dgrp1)[i].resize(info->ublcount[i]);
	}
	(module->ubl2dgrp2) = (module->ubl2dgrp1);
	return;
}

float square(float x){
	return pow( max(min(x, MAX_POW_2), -MAX_POW_2), 2);
}

float amp(float x, float y){
	return sqrt( square(x)  + square(y) );
}

float amp(complex float cx){
    return sqrt(crealf(cx*conjf(cx)));
}

float amp(cfloat cx){
    return cmod(cx);
}

float phase(float re, float im){
	/*if (re == 0 and im == 0){
		return 0;
	}*/
	return atan2(im, re);
}

float phase(complex float cx){
    return cargf(cx);
}

float phase(cfloat cx){
    return carg(cx);
}

// BLAS version, maybe use arrays in stead of vectors.
float norm(vector<complex float> *vec){
    return cblas_scnrm2(vec -> size(), vec -> data(), 1);
}

// This will be made into a GPU version.
float norm(vector<cfloat> *vec){
    float tmp = 0;
    for (int i=0; i<vec->size(); i++)
        tmp += pow(cmod(vec->at(i)), 2);
    return pow(tmp, 0.5);
}

float phaseWrap(float x, float offset/*default -pi*/){
	while ( x <= offset ){
		x += 2 * PI;
	}
	while ( x > offset + 2 * PI ){
		x -= 2 * PI;
	}
	return x;
}

float median(vector<float> list){
	int l = list.size();
	if (l == 0) return 0;
	sort(list.begin(), list.end());
	int index = floor( l / 2 );
	if(l % 2 == 1) return list[index];
	else return (list[index] + list[index - 1]) / 2;
}

float medianAngle(vector<float> *list){
	string METHODNAME = "medianAngle";
	//cout << "#!#" << FILENAME << "#!#" << METHODNAME << " DBG ";
	vector<float> xList(list->size());
	vector<float> yList(list->size());
	for (unsigned int i = 0; i < list->size(); i++){
		//cout << list[i] << " ";
		xList[i] = cos(list->at(i));
		yList[i] = sin(list->at(i));
	}
	//cout << " median is " << atan2(median(yList), median(xList)) << endl;
	return atan2(median(yList), median(xList));
}

// XXX I want one where the input matrix is a 1D vector as to not waste
// time/space copying the input vector into a new array.
void vecmatmul(vector<vector<float> > * Afitting, vector<float> * v, vector<float> * ampfit){
    int i, j; // looping indices
	int nrows = Afitting -> size(), ncols = v -> size(); // Array dimensions
	float *in_v = v -> data(), *out_v = ampfit -> data(); // in/out vectors
	float *Amatrix = new float[nrows*ncols]; // matrix that's nice for cblas

	// Create a 1-D matrix that's nice for BLAS
	for (i=0; i<nrows; i++){
        for (j=0; j<ncols; j++)
            Amatrix[ncols*i + j] = (Afitting -> at(i))[j];
    }

    // Run the BLAS function
    cblas_sgemv(CblasRowMajor, CblasNoTrans,
                nrows, ncols, 1, Amatrix, ncols, in_v, 1, 0, out_v, 1);

    // Cleanup
    delete Amatrix;
	return;
}

// just going to overload this now and take out the first one later.
void vecmatmul(vector<float> * Afitting, vector<float> * v, vector<float> * ampfit){
	int nrows = ampfit -> size(), ncols = v -> size(); // Array dimensions
	float *in_v = v -> data(), *out_v = ampfit -> data(); // in/out vectors
	float *Amatrix = Afitting -> data();

    // Run the BLAS function
    cblas_sgemv(CblasRowMajor, CblasNoTrans,
                nrows, ncols, 1, Amatrix, ncols, in_v, 1, 0, out_v, 1);
	return;
}

complex float minimizecomplex(vector<complex float> *a, vector<complex float> *b){
    // A = B*c where A and B complex vecs, c complex number, solve for c
    // A            = B*c
    // A * {B}      = B * {B} * c
    // sum(A * {B}) = sum(B * {B}) * c
    // c            = sum(A * {B}) / sum(B * {B})
	complex float sum1 = 0;
	for (uint i =0; i < a->size(); i++){
	    sum1 += (a -> at(i)) * conjf(b -> at(i));
	}

	return sum1 / pow(norm(b), 2);
}

cfloat minimizecomplex(vector<cfloat> *a, vector<cfloat> *b){
    // A = B*c where A and B complex vecs, c complex number, solve for c
    // A            = B*c
    // A * {B}      = B * {B} * c
    // sum(A * {B}) = sum(B * {B}) * c
    // c            = sum(A * {B}) / sum(B * {B})
	cfloat sum1 = {0, 0}, sum2 = {0, 0};
	for (uint i =0; i < a->size(); i++){
	    sum1 = cadd(sum1, cmult(a -> at(i), conj(b -> at(i))));
	}
	sum2.x = pow(norm(b), 2);

	return cdiv(sum1, sum2);
}

void logcaladd(vector<cfloat> *data,
               vector<cfloat> *additivein,
               redundantinfo* info,
               vector<float>* calpar,
               vector<cfloat> *additiveout,
               int computeUBLFit,
               int compute_calpar,
               calmemmodule_clh* module){
    // if computeUBLFit is 1, compute the ubl estimates given data and calpars,
    // rather than read ubl estimates from input
	int nubl = info->ublindex.size();
    int ai, aj; // antenna indices

	////initialize data and g0 ubl0
	for (unsigned int b = 0; b < (module->cdata1).size(); b++){
		module->cdata1[b] = csub(data->at(b), additivein->at(b));
	}

	float amptmp;
	unsigned int cbl;
	for (int a = 0; a < info->nAntenna; a++){
		amptmp = pow(10, calpar->at(3 + a));
		module->g0[a].x = amptmp * cos(calpar->at(3 + info->nAntenna + a));
		module->g0[a].y = amptmp * sin(calpar->at(3 + info->nAntenna + a));
	}

	if (computeUBLFit != 1){
		for (int u = 0; u < nubl; u++){
			module->ubl0[u].x = calpar->at(3 + 2 * info->nAntenna + 2 * u);
			module->ubl0[u].y = calpar->at(3 + 2 * info->nAntenna + 2 * u + 1);
		}
	} else{
	    // if computeUBLFit is 1, compute the ubl estimates given data and
	    // calpars, rather than read ubl estimates from input
		for (int u = 0; u < nubl; u++){
			for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
				cbl = info->ublindex[u][i];
                ai = info->bl2d[cbl][0]; aj = info->bl2d[cbl][1];
				module->ubl2dgrp1[u][i] = module->cdata1[cbl];
				module->ubl2dgrp2[u][i] = cmult(conj(module->g0[ai]), module->g0[aj]);
			}

			module->ubl0[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
		}
	}

	int nant = info->nAntenna;
	int ncross = info->bl2d.size();
	////read in amp and args
	for (int b = 0; b < ncross; b++){
		ai = info->bl2d[b][0];
		aj = info->bl2d[b][1];
		cfloat cxdiff = csub(data->at(b), additivein->at(b));
		if (cxdiff.x == 0 && cxdiff.y == 0){//got 0, quit
			for(int i = 3; i < 3 + 2 * nant + 2 * nubl; i++){
				calpar->at(i) = 0;
			}
			calpar->at(1) = INFINITY;
			return;
		}

		module->amp1[b] = log10(amp(csub(data->at(b), additivein->at(b)))) - calpar->at(3+ai) - calpar->at(3+aj);
		module->pha1[b] = phaseWrap(phase(csub(data->at(b), additivein->at(b))) + calpar->at(3 + nant + ai) - calpar->at(3 + nant + aj));
	}

	////rewrap args//TODO: use module->ubl0
	for(int i = 0; i < nubl; i ++){
		for (uint j = 0; j < (module->ublgrp1)[i].size(); j ++){
			(module->ublgrp1)[i][j] = module->pha1[info->ublindex[i][j]];
		}
	}

	for (int i = 0; i < nubl; i++){
		module->ubl1[i].x = creal(module->ubl1[i]);
		module->ubl1[i].y = medianAngle(&((module->ublgrp1)[i]));
	}

	for (int b = 0; b < ncross; b++) {
		module->pha1[b] = phaseWrap(module->pha1[b], cimag((module->ubl1)[info->bltoubl[b]]) - PI);
	}

	fill(module->x3.begin(), module->x3.end(), 0);////At.y
	for (unsigned int i = 0; i < info->Atsparse.size(); i++){
		for (unsigned int j = 0; j < info->Atsparse[i].size(); j++){
			module->x3[i] += module->amp1[info->Atsparse[i][j]];
		}
	}

	fill(module->x4.begin(), module->x4.end(), 0);////Bt.y
	for (unsigned int i = 0; i < info->Btsparse.size(); i++){
		for (unsigned int j = 0; j < info->Btsparse[i].size(); j++){
			module->x4[i] += module->pha1[info->Btsparse[i][j][0]] * info->Btsparse[i][j][1];
		}
	}
	vecmatmul(&(info->AtAi), &(module->x3), &(module->x1));
	vecmatmul(&(info->BtBi), &(module->x4), &(module->x2));

	for(int b = 0; b < ncross; b++) {
		ai = info->bl2d[b][0];
		aj = info->bl2d[b][1];
		float amp = pow(10, module->x1[nant + info->bltoubl[b]] + module->x1[ai] + module->x1[aj] + calpar->at(3 + ai) + calpar->at(3 + aj));
		float phase =  module->x2[nant + info->bltoubl[b]] - module->x2[ai] + module->x2[aj] - calpar->at(3 + nant + ai) + calpar->at(3 + nant + aj);
		cfloat tmpvar = {amp*cos(phase), amp*sin(phase)};
		additiveout->at(b) = csub(data->at(b), tmpvar);
	}

	if(compute_calpar == 0){////compute additive term only
		calpar->at(1) = pow(norm(additiveout), 2);
		//cout << norm(additiveout) << endl;
		return;
	} else if(compute_calpar == 1){////compute full set of calpars
		for(int a = 0; a < nant; a++){
			calpar->at(3 + a) += module->x1[a];
			calpar->at(3 + nant + a) += module->x2[a];
		}
		for(int u = 0; u < nubl; u++){
			calpar->at(3 + 2 * nant + 2 * u) = pow(10, module->x1[nant + u]) * cos(module->x2[nant + u]);
			calpar->at(3 + 2 * nant + 2 * u + 1) = pow(10, module->x1[nant + u]) * sin(module->x2[nant + u]);
		}
		calpar->at(1) = pow(norm(additiveout), 2);
	}
	return;
}

void lincal(vector<cfloat> *data,
            vector<cfloat> *additivein,
            redundantinfo *info,
            vector<float> *calpar,
            vector<cfloat> *additiveout,
            int computeUBLFit,
            calmemmodule_clh* module,
            float convergethresh,
            int maxiter,
            float stepsize){
	int nubl = info->ublindex.size();
    int ai, aj; // antenna indices
	////initialize data and g0 ubl0
	for (unsigned int b = 0; b < (module->cdata1).size(); b++){
		module->cdata1[b] = csub(data->at(b), additivein->at(b));
	}
	float amptmp;
	unsigned int cbl;
	float stepsize2 = 1 - stepsize;
	for (int a = 0; a < info->nAntenna; a++){
		amptmp = pow(10, calpar->at(3 + a));
		module->g0[a].x = amptmp * cos(calpar->at(3 + info->nAntenna + a));
		module->g0[a].y = amptmp * sin(calpar->at(3 + info->nAntenna + a));
	}
	if (computeUBLFit != 1){
		for (int u = 0; u < nubl; u++){
			module->ubl0[u].x = calpar->at(3 + 2 * info->nAntenna + 2 * u);
			module->ubl0[u].y = calpar->at(3 + 2 * info->nAntenna + 2 * u + 1);
		}
	} else{
	    // if computeUBLFit is 1, compute the ubl estimates given data and
	    // calpars, rather than read ubl estimates from input
		for (int u = 0; u < nubl; u++){
			for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
				cbl = info->ublindex[u][i];
                ai = info->bl2d[cbl][0]; aj = info->bl2d[cbl][1];
				module->ubl2dgrp1[u][i] = module->cdata1[cbl];
				module->ubl2dgrp2[u][i] = cmult(conj(module->g0[ai]), module->g0[aj]);
			}

			module->ubl0[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
		}
	}

	float starting_chisq, chisq, chisq2;
	cfloat gain;
	int a1, a2; // antenna indices
	chisq = 0;
	for (unsigned int b = 0; b < (module->cdata2).size(); b++){
		a1 = info->bl2d[b][0];
		a2 = info->bl2d[b][1];
		gain = cmult(conj(module->g0[a1]), module->g0[a2]);
		//printf("gain: %g + %gi\n", gain.x, gain.y);
		//printf("module->ubl0: %g + %gi\n", module->ubl0[info->bltoubl[b]].x, module->ubl0[info->bltoubl[b]].y);
		module->cdata2[b] = csub(cmult(gain, module->ubl0[info->bltoubl[b]]), module->cdata1[b]);
	}
	chisq = pow(norm(&(module->cdata2)), 2);
	starting_chisq = chisq;

	////start iterations
	int iter = 0;
	float componentchange = 100;
	while(iter < maxiter and componentchange > convergethresh){
		iter++;

		for (unsigned int a3 = 0; a3 < module->g3.size(); a3++){
		    // g3 will be containing the final dg, g1, g2 will contain a and b
		    // as in the cost function LAMBDA = ||a + b*g||^2
			for (unsigned int a = 0; a < module->g3.size(); a++){
				cbl = info->bl1dmatrix[a3][a];
                // cbl is unsigned, so gauranteed not < 0
				if (cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){
				    //badbl or ubl has only 1 bl
					module->g1[a].x = 0; module->g1[a].y = 0;
					module->g2[a].x = 0; module->g2[a].y = 0;
				}else if(info->bl2d[cbl][1] == a3){
					module->g1[a] = module->cdata1[cbl];
					module->g2[a] = cmult(conj(module->g0[a]), module->ubl0[info->bltoubl[cbl]]);
				}else{
				    module->g1[a] = conj(module->cdata1[cbl]);
				    module->g2[a] = cmult(conj(module->g0[a]), conj(module->ubl0[info->bltoubl[cbl]]));
				}
			}
			module->g3[a3] = minimizecomplex(&(module->g1), &(module->g2));
		}

		////ubl M
		for (int u = 0; u < nubl; u++){
			for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
				cbl = info->ublindex[u][i];
                ai = info->bl2d[cbl][0]; aj = info->bl2d[cbl][1];
                module->ubl2dgrp1[u][i] = module->cdata1[cbl];
                module->ubl2dgrp2[u][i] = cmult(conj(module->g0[ai]), module->g0[aj]);
			}

			module->ubl3[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
		}


		// Update g and ubl, do not update single-bl bls since they are not
		// reversible. Will reverse this step later is chisq increased
		// float fraction;
		for (unsigned int a = 0; a < module->g3.size(); a++){
			module->g0[a] = cadd(cmult(module->g0[a], stepsize2), cmult(module->g3[a], stepsize));
		}
		for (unsigned int u = 0; u < module->ubl3.size(); u++){
			if ((info->ublcount)[u] > 1){
				module->ubl0[u] = cadd(cmult(module->ubl0[u], stepsize2),
				                       cmult(module->ubl3[u], stepsize));
			}
		}

		//compute chisq and decide convergence
		chisq2 = 0;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			if ((info->ublcount)[info->bltoubl[b]] > 1){
			    // automatically use 0 for single-bl ubls, their actual values
			    // are not updated yet
				a1 = info->bl2d[b][0];
				a2 = info->bl2d[b][1];
				gain = cmult(conj(module->g0[a1]), module->g0[a2]);
				module->cdata2[b] = csub(cmult(gain, module->ubl0[info->bltoubl[b]]),
				                         module->cdata1[b]);
			} else {
			    module->cdata2[b].x = 0;
			    module->cdata2[b].y = 0;
            }
		}
		chisq2 = pow(norm(&(module->cdata2)), 2);
		componentchange = (chisq - chisq2) / chisq;

        if (componentchange > 0){
            // if improved, keep g0 and ubl0 updates, and update single-bl ubls 
            // and chisq
			chisq = chisq2;
			for (unsigned int u = 0; u < module->ubl3.size(); u++){
			    // make sure there's no error on unique baselines with only 1 
			    // baseline
				for (unsigned int i = 0; i < module->ubl2dgrp1[u].size(); i++){
					cbl = info->ublindex[u][i];
                    ai = info->bl2d[cbl][0]; aj = info->bl2d[cbl][1];
                    module->ubl2dgrp1[u][i] = module->cdata1[cbl];
                    module->ubl2dgrp2[u][i] = cmult(conj(module->g0[ai]), module->g0[aj]);
				}
				module->ubl3[u] = minimizecomplex(&(module->ubl2dgrp1[u]), &(module->ubl2dgrp2[u]));
				module->ubl0[u] = module->ubl3[u];
			}
		} else {//reverse g0 and ubl0 changes
			iter--;
			for (unsigned int a = 0; a < module->g3.size(); a++){
				module->g0[a] = csub(module->g0[a], cmult(module->g3[a], stepsize));
				module->g0[a] = cdiv(module->g0[a], stepsize2);

			}
			for (unsigned int u = 0; u < module->ubl3.size(); u++){
				if ((info->ublcount)[u] > 1){
					module->ubl0[u] = csub(module->ubl0[u], cmult(module->ubl3[u], stepsize));
					module->ubl0[u] = cdiv(module->ubl0[u], stepsize2);
				}
			}
		}
	}


	////update calpar and additive term
	if(componentchange > 0 or iter > 1){
		for (unsigned int a = 0; a < module->g0.size(); a++){
			calpar->at(3 + a) = log10(amp(module->g0[a]));
			calpar->at(3 + info->nAntenna + a) = phase(module->g0[a]);
		}
		int tmp = 3 + 2 * info->nAntenna;
		for (unsigned int u = 0; u < module->ubl0.size(); u++){
			calpar->at(tmp + 2 * u) = creal(module->ubl0[u]);
			calpar->at(tmp + 2 * u + 1) = cimag(module->ubl0[u]);
		}

		calpar->at(0) += iter;
		calpar->at(2) = chisq;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			additiveout->at(b) = cneg(module->cdata2[b]);
		}

        // Create a chisq for each antenna. Right now, this is done at every
        // time and frequency, since that's how lincal is called, but that can
        // be summed later over all times and frequencies to get a chisq for
        // each antenna.
        float delta;
        int chisq_ant = 3 + 2*(info -> nAntenna + nubl);
        for (int b = 0; b < (module -> cdata2).size(); b++){
		    delta = pow(amp(module->cdata2[b]), 2);
            a1 = info->bl2d[b][0];
            a2 = info->bl2d[b][1];
            calpar -> at(chisq_ant + a1) += delta;
            calpar -> at(chisq_ant + a2) += delta;
        }
	
	}else{////if chisq didnt decrease, keep everything untouched
		calpar->at(0) += 0;
		calpar->at(2) = starting_chisq;
        // XXX do we need to put a dummy value in for chisq per ant?
	}

	return;
}

void gaincal(vector<complex float> *data,
             vector<complex float> *additivein,
             redundantinfo* info,
             vector<float>* calpar,
             vector<complex float> *additiveout,
             calmemmodule* module,
             float convergethresh,
             int maxiter,
             float stepsize){
    int nubl = info->ublindex.size();

	////initialize data and g0 ubl0
	for (unsigned int b = 0; b < (module->cdata1).size(); b++){
		module->cdata1[b] = data->at(b) - additivein->at(b);
	}
	float amptmp;
	unsigned int cbl;
	float stepsize2 = 1 - stepsize;
	for (int a = 0; a < info->nAntenna; a++){
		amptmp = pow(10, calpar->at(3 + a));
		module->g0[a]  = amptmp * cos(calpar->at(3 + info->nAntenna + a));
		module->g0[a] += amptmp*I*sin(calpar->at(3 + info->nAntenna + a));
	}

	for (int u = 0; u < nubl; u++){
		module->ubl0[u] = 1;
	}


	float starting_chisq, chisq, chisq2, delta;
	complex float gain, difftmp;
	int a1, a2;
	chisq = 0;
	for (unsigned int b = 0; b < (module->cdata2).size(); b++){
		a1 = info->bl2d[b][0];
		a2 = info->bl2d[b][1];
		gain = conjf(module->g0[a1]) * module->g0[a2];
		module->cdata2[b] = gain * module->ubl0[info->bltoubl[b]];
		difftmp = module->cdata2[b] - module->cdata1[b];
		delta = crealf(difftmp * conjf(difftmp));
		chisq += delta;
	}
	starting_chisq = chisq;

	////start iterations
	int iter = 0;
	float componentchange = 100;
	while(iter < maxiter and componentchange > convergethresh){
		iter++;

		for (unsigned int a3 = 0; a3 < module->g3.size(); a3++){
		    // g3 will be containing the final dg, g1, g2 will contain a and b
		    // as in the cost function LAMBDA = ||a + b*g||^2
			for (unsigned int a = 0; a < module->g3.size(); a++){
				cbl = info->bl1dmatrix[a3][a];
                // cbl is unsigned and so gauranteed >= 0
				if (cbl > module->cdata1.size() or info->ublcount[info->bltoubl[cbl]] < 2){//badbl or ubl has only 1 bl
					module->g1[a] = 0;
					module->g2[a] = 0;
				}else if(info->bl2d[cbl][1] == a3){
					module->g1[a] = module->cdata1[cbl];
					module->g2[a] = conjf(module->g0[a]) * module->ubl0[info->bltoubl[cbl]];
				}else{
				    module->g1[a] = conjf(module->cdata1[cbl]);
				    module->g2[a] = conjf(module->g0[a]) * conjf(module->ubl0[info->bltoubl[cbl]]);
				}
			}
			module->g3[a3] = minimizecomplex(&(module->g1), &(module->g2));
		}


		// Update g and ubl, do not update single-bl bls since they are not
		// reversible. Will reverse this step later is chisq increased
		// float fraction;
		for (unsigned int a = 0; a < module->g3.size(); a++){
			module->g0[a] = stepsize2 * module->g0[a] + stepsize * module->g3[a];
		}

		//compute chisq and decide convergence
		chisq2 = 0;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			if ((info->ublcount)[info->bltoubl[b]] > 1){
			    // automatically use 0 for single-bl ubls, their actaul values 
			    // are not updated yet
				a1 = info->bl2d[b][0];
				a2 = info->bl2d[b][1];
				gain = conjf(module->g0[a1]) * module->g0[a2];
				module->cdata2[b] = gain * module->ubl0[info->bltoubl[b]];
		        difftmp = module->cdata2[b] - module->cdata1[b];
		        chisq2 += crealf(difftmp * conjf(difftmp));
			}
		}
		componentchange = (chisq - chisq2) / chisq;

		if (componentchange > 0){
		    // if improved, keep g0 and ubl0 updates, and update single-bl ubls
		    // and chisq
			chisq = chisq2;
		} else {//reverse g0 and ubl0 changes
			iter--;
			for (unsigned int a = 0; a < module->g3.size(); a++){
				module->g0[a] = (module->g0[a] - stepsize * module->g3[a]) / stepsize2;
			}
		}
	}

	////update calpar and additive term
	if(componentchange > 0 or iter > 1){
		for (unsigned int a = 0; a < module->g0.size(); a++){
			calpar->at(3 + a) = log10(amp(module->g0[a]));
			calpar->at(3 + info->nAntenna + a) = phase(module->g0[a]);
		}

		calpar->at(0) += iter;
		calpar->at(2) = chisq;
		for (unsigned int b = 0; b < (module->cdata2).size(); b++){
			additiveout->at(b) = module->cdata1[b] - module->cdata2[b];
		}
	}else{////if chisq didnt decrease, keep everything untouched
		calpar->at(0) += 0;
		calpar->at(2) = starting_chisq;
	}
	return;
}

void removeDegen(vector<float> *calpar, redundantinfo *info, calmemmodule_clh* module){
    //forces the calibration parameters to have average 1 amp, and no shifting
    //the image in phase. Note: 1) If you have not absolute calibrated the data,
    //there's no point in running this, because this can only ensure that the
    //calpars don't screw up already absolute calibrated data. 2) the antloc and
    //ubl stored in redundant info must be computed from idealized antloc,
    //otherwise the error in antloc from perfect redundancy will roll into this
    //result, in an unknown fashion.
	////load data
    int nubl = info->ublindex.size();
	vector<float> pha1(info->nAntenna, 0);
	for (int a = 0 ; a < info->nAntenna; a ++){
		pha1[a] = calpar->at(3 + info->nAntenna + a);
	}
	for (int u = 0 ; u < nubl; u ++){
		module->ubl1[u].x = amp(calpar->at(3 + 2 * info->nAntenna + 2 * u), calpar->at(3 + 2 * info->nAntenna + 2 * u + 1));
		module->ubl1[u].y = phase(calpar->at(3 + 2 * info->nAntenna + 2 * u), calpar->at(3 + 2 * info->nAntenna + 2 * u + 1));
	}

	////compute amp delta
	float ampfactor = 0;
	// average |g|, divide ant calpar by this, multiply ublfit by square this
	for (int a = 0 ; a < info->nAntenna; a ++){
		ampfactor += pow(10, calpar->at(3 + a));
	}
	ampfactor = ampfactor / info->nAntenna;

	////compute phase delta
	vecmatmul(&(info->degenM), &(pha1), &(module->x1));
	// x1: add ant calpar and ubl fit by this

	////correct ant calpar
	for (int a = 0 ; a < info->nAntenna; a ++){
		calpar->at(3 + a) = calpar->at(3 + a) - log10(ampfactor);
		calpar->at(3 + info->nAntenna + a) = phaseWrap(calpar->at(3 + info->nAntenna + a) + module->x1[a]);
	}

	////correct ublfit
	for (int u = 0 ; u < nubl; u ++){
		module->ubl2[u].x = creal(module->ubl1[u]) * ampfactor * ampfactor;
		module->ubl2[u].y = cimag(module->ubl1[u]) + module->x1[info->nAntenna + u];
		calpar->at(3 + 2 * info->nAntenna + 2 * u) = creal(module->ubl2[u]) * cos(cimag(module->ubl2[u]));
		calpar->at(3 + 2 * info->nAntenna + 2 * u + 1) = creal(module->ubl2[u]) * sin(cimag(module->ubl2[u]));
	}
	return;
}
