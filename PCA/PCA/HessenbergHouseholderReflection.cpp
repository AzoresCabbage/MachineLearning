#include "stdafx.h"
#include "HessenbergHouseholderReflection.h"


Matrix HessenbergHouseholderReflection::getIdentityMatrix(const int size)
{
	Matrix res(size,size);
	for(int i=0;i<size;++i)
		res[i][i] = 1;
	return res;
}

Matrix HessenbergHouseholderReflection::getHRatCol(const Matrix& M,const int col)
{
	//Matrix H = I - 2*P
	Matrix I = getIdentityMatrix(M.col);
	Matrix P = getMatrixPatCol(M,col);
	Matrix P2 = P*2;
	Matrix H = I - P2;
	return H;
}

Matrix HessenbergHouseholderReflection::getMatrixPatCol(const Matrix& M,const int col)
{
	//Matrix p = v*v^T / v^T*v
	Matrix v = getMatrixVatCol(M,col);//return a vector
	Matrix vt = v.transpose();
	Matrix vMulvt = v * vt;//is a value stored in mat[0][0]
	Matrix vtMulv = vt * v;
	if(vtMulv[0][0] == 0)
	{
		printf("Vector v multiply v^T = 0!");
		exit(-1);
	}
	Matrix P = vMulvt * (1.0/vtMulv[0][0]);
	return P;
}

Matrix HessenbergHouseholderReflection::getMatrixVatCol(const Matrix& M,const int col)
{
	//Matrix v = w - x
	Matrix x = getXatCol(M,col);
	Matrix w = getWatCol(x,col);
	Matrix v = getVatCol(x,w);
	return v;
}

Matrix HessenbergHouseholderReflection::getXatCol(const Matrix& M,const int col)
{
	Matrix ret(M.row,1);
	for(int i=0;i<M.row;++i)
		if(i > col)
			ret[i][0] = M[i][col];
	return ret;
}

Matrix HessenbergHouseholderReflection::getVatCol(const Matrix& x,const Matrix& w)
{
	int sign = 1;
	if(x[0][0] > 0) sign = -1;
	Matrix v(x.row,x.col);
	for(int i=0;i<v.row;++i)
		v[i][0] = x[i][0] + w[i][0];
	return v;
}

Matrix HessenbergHouseholderReflection::getWatCol(const Matrix& M,const int col)
{
	Matrix w(M.row,1);
	double tmp = 0;
	for(int i=0;i<M.row;++i)
		tmp += M[i][0]*M[i][0];
	tmp = sqrt(tmp);
	w[col+1][0] = tmp;
	return w;
}