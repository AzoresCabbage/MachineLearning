#pragma once
#include "stdafx.h"
#include "Matrix.h"

class HessenbergHouseholderReflection{
public:
	Matrix getIdentityMatrix(const int size);
	Matrix getHRatCol(const Matrix& M,const int col);
	Matrix getMatrixPatCol(const Matrix& M,const int col);
	Matrix getMatrixVatCol(const Matrix& M,const int col);
	Matrix getXatCol(const Matrix& M,const int col);
	Matrix getVatCol(const Matrix& x,const Matrix& w);
	Matrix getWatCol(const Matrix& M,const int col);
};