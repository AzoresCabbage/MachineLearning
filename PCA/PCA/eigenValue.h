#pragma once
#include "stdafx.h"
#include "Matrix.h"
#include "HessenbergHouseholderReflection.h"
#include <vector>

#define eps 1e-5

class eigenValue
{
private:
	HessenbergHouseholderReflection hessenberg;
	Matrix convertToHessenberg(const Matrix&);
	Matrix getRForAStep(const Matrix& M,const int i);
	bool isConvergence(const Matrix&);
public:
	eigenValue(){}
	std::vector <double> getEigenValues(const Matrix&);
};