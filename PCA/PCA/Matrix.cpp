#include "stdafx.h"
#include "Matrix.h"

Matrix::Matrix(const int n,const int m):row(n),col(m)
{
	mat = new double*[n];
	for(int i=0;i<n;++i)
	{
		mat[i] = new double[m];
		memset(mat[i],0,sizeof(double)*m);
	}
}

Matrix::Matrix(const int n,const int m,double** a):row(n),col(m)
{
	mat = new double*[n];
	for(int i=0;i<n;++i)
	{
		mat[i] = new double[m];
		memcpy(mat[i],a[i],sizeof(double)*m);
	}
}

Matrix::~Matrix()
{
	for(int i=0;i<row;++i)
		delete[] mat[i];
	delete[] mat;
}

Matrix::Matrix(const Matrix& a):row(a.row),col(a.col)
{
	mat = new double*[row];
	for(int i=0;i<row;++i)
	{
		mat[i] = new double[col];
		memcpy(mat[i],a.mat[i],sizeof(double)*col);
	}
}

Matrix& Matrix::operator=(const Matrix& a)
{
	double** tmp = mat;

	mat = new double*[a.row];
	for(int i=0;i<a.row;++i)
	{
		mat[i] = new double[a.col];
		memcpy(mat[i],a.mat[i],sizeof(double)*a.col);
	}
	for(int i=0;i<row;++i)
		delete[] tmp[i];
	delete[] tmp;

	row = a.row;
	col = a.col;
	return *this;
}

Matrix Matrix::transpose()
{
	Matrix ret(col,row);
	for(int i=0;i<col;++i)
		for(int j=0;j<row;++j)
			ret.mat[i][j] = mat[j][i];
	return ret;
}

void Matrix::print()
{
	for(int i=0;i<row;++i)
	{
		for(int j=0;j<col;++j)
			printf("%f ",mat[i][j]);
		puts("");
	}
	puts("");
}

Matrix Matrix::operator+(const Matrix& a)const
{
	if(row != a.row || col != a.col)
	{
		printf("Matrix(%d,%d) can't add Matrix(%d,%d)!\n",row,col,a.row,a.col);
		exit(-1);
	}
	Matrix ret(a);
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
			ret.mat[i][j] += mat[i][j];
	return ret;
}

Matrix Matrix::operator-(const Matrix& a)const
{
	if(row != a.row || col != a.col)
	{
		printf("Matrix(%d,%d) can't minus Matrix(%d,%d)!\n",row,col,a.row,a.col);
		exit(-1);
	}
	Matrix ret(*this);
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
			ret.mat[i][j] -= a.mat[i][j];
	return ret;
}

Matrix Matrix::operator*(const Matrix& a)const
{
	if(col != a.row)
	{
		printf("Matrix(%d,%d) can't multiply with Matrix(%d,%d)!\n",row,col,a.row,a.col);
		exit(-1);
	}
	Matrix ret(row,a.col);
	
	/*for(int i=0;i<row;i++)
		for(int k=0;k<col;k++){
			if(mat[i][k] == 0) continue;
			for(int j=0;j<a.col;j++)
				ret.mat[i][j] += mat[i][k] * a.mat[k][j];
		}
		*/
	for(int i=0;i<row;++i)
		for(int j=0;j<a.col;++j)
			for(int k=0;k<col;++k)
				ret.mat[i][j] += mat[i][k]*a.mat[k][j];
	return ret;
}

Matrix Matrix::operator*(const double a)const
{
	Matrix ret(*this);
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
			ret.mat[i][j] *= a;
	return ret;
}

bool Matrix::operator==(const Matrix& a)const
{
	if(row != a.row || col != a.col) return false;
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
			if(a.mat[i][j] != mat[i][j])
				return false;
	return true;
}

bool Matrix::operator!=(const Matrix& a)const
{
	return !operator==(a);
}

Matrix Matrix::operator-(const double a)const
{
	Matrix ret(*this);
	for(int i=0;i<row;++i)
		for(int j=0;j<col;++j)
			ret.mat[i][j] -= a;
	return ret;
}