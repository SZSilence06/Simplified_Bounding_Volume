#ifndef FASTLAP_H
#define FASTLAP_H

#include "mulStruct.h"
#include "mulGlobal.h"

int fastlap(int* plhsSize, int* prhsSize, int* pnumSing, double* px, int* pshape, int* pdtype, int* plhsType, int* prhsType, int* plhsIndex, int* prhsIndex,
            double* plhsVect, double* prhsVect, double* pxf, double* pxnrm, int* pnumLev, int* pnumMom, int* pmaxItr, double* ptol, int* pjob);

#endif
