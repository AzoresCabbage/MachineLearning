#include "stdafx.h"
#include "eigenValue.h"

bool eigenValue::isConvergence(const Matrix& M)
{
	for(int i=1;i<M.row;++i)
		if(fabs(M[i][i-1]) >= eps)
			return false;
	return true;
}

Matrix eigenValue::getRForAStep(const Matrix& R,const int i)
{
	Matrix partR = hessenberg.getIdentityMatrix(R.row);
	double r = R[i][i]*R[i][i] + R[i+1][i]*R[i+1][i];
	r = sqrt(r);
	partR[i][i] = partR[i+1][i+1] = R[i][i]/r;
	partR[i][i+1] = R[i+1][i]/r;
	partR[i+1][i] = -R[i+1][i]/r;
	return partR;
}

Matrix eigenValue::convertToHessenberg(const Matrix& A)
{
	Matrix R(A);
	for(int i=0,sz=A.row-2;i<sz;++i)
	{
		Matrix H = hessenberg.getHRatCol(R,i);
		if(H.row == 0) return R;
		R = H * R;
		R = R * H;
	}
	return R;
}

//QR decomposition to get all eigenvalues
std::vector <double> eigenValue::getEigenValues(const Matrix& M)
{
	Matrix hessenMat = convertToHessenberg(M);
	Matrix B(hessenMat);
	int row = hessenMat.row;
	int col = hessenMat.col;
	while(true)
	{
		Matrix Q = hessenberg.getIdentityMatrix(row);
		Matrix R(B);
		for(int i=0;i<row-1;++i)
		{
			Matrix partR = getRForAStep(R,i);
			Matrix tmp = partR * R;
			R = partR * R;
			Matrix transposeR = partR.transpose();
			Q = Q * transposeR;
		}
		B = R*Q;
		if(isConvergence(B))
		{
			std::vector <double> ret;
			for(int i=0;i<row;++i)
				ret.push_back(B[i][i]);
			return ret;
		}
	}
}