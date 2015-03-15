#pragma once
#include "stdafx.h"

class Matrix{
public:
	int row,col;
	double** mat;
	Matrix(const int n,const int m);
	Matrix(const int n,const int m,double** matrix);
	Matrix(const Matrix& a);//¿½±´¹¹Ôìº¯Êý
	Matrix& operator = (const Matrix& a);
	
	~Matrix();

	Matrix transpose();
	void print();

	Matrix operator + (const Matrix& a)const;
	Matrix operator - (const Matrix& a)const;
	Matrix operator - (const double a)const;
	Matrix operator * (const Matrix& a)const;
	Matrix operator * (const double a)const;
	bool operator == (const Matrix& a)const;
	bool operator != (const Matrix& a)const;
};