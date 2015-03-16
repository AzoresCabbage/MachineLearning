// PCA.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "eigenValue.h"
#include "Matrix.h"
#include <algorithm>
#include <functional>
#include <fstream>
using namespace std;

int n,m,k;//n sample,every sample has m attributes,compress to k dim

//reference:<Introduction to Algorithms> ch28
//input:initial matrix A, a instance of matrix L,U,P
//output:the pi array
int* LUP_Decomposition(Matrix A,Matrix& L,Matrix& U,Matrix& P)
{
	int n = A.row;
	int *pi = new int[n];
	for(int i=0;i<n;++i)
		pi[i] = i;
	for(int k=0;k<n;++k)
	{
		double p = 0;
		int kk = 0;
		for(int i=k;i<n;++i)
		{
			if(fabs(A[i][k]) > p)
				p = fabs(A[i][k]) , kk=i;
		}
		if(p == 0)
		{
			printf("sigular matrix!\n");
			delete[] pi;
			exit(-1);
		}
		swap(pi[k],pi[kk]);
		for(int i=0;i<n;++i)
			swap(A[k][i],A[kk][i]);
		for(int i=k+1;i<n;++i)
		{
			A[i][k] /= A[k][k];
			for(int j=k+1;j<n;++j)
				A[i][j] -= A[i][k]*A[k][j];
		}
	}
	for(int i=0;i<n;++i)
		P[i][pi[i]] = 1;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
		{
			if(i > j)
				L[i][j] = A[i][j];
			else
				U[i][j] = A[i][j];
		}
		L[i][i] = 1;
	}
	return pi;
	//delete[] pi;
}

//reference:<Introduction to Algorithms> ch28
//input:a instance of matrix L,U. pi array, equation's answer vector b
//output:the solution vector
Matrix LUP_SolveEquationGroup(const Matrix& L,const Matrix& U,const int* const pi,const Matrix& b)
{
	int n = L.row;
	Matrix x(n,1);
	Matrix y(n,1);
	for(int i=0;i<n;++i)
	{
		y[i][0] = b[pi[i]][0];
		for(int j=0;j<i;++j)
		{
			y[i][0] -= L[i][j]*y[j][0];
		}
	}
	for(int i=n-1;i>=0;--i)
	{
		x[i][0] = y[i][0];
		for(int j=n-1;j>i;--j)
			x[i][0] -= U[i][j]*x[j][0];
		x[i][0] /= U[i][i];
	}
	return x;
}

//reference:<Introduction to Algorithms> ch28
//input:a matrix A
//output:the inverse matrix of A
Matrix getInvMatrix(Matrix A)
{
	int n = A.row;
	Matrix InvA(n,n);
	Matrix L(n,n);
	Matrix U(n,n);
	Matrix P(n,n);
	int *pi = LUP_Decomposition(A,L,U,P);
	for(int i=0;i<n;++i)
	{
		Matrix b(n,1);
		for(int j=0;j<n;++j)
			b[j][0] = i == j?1:0;
		Matrix x = LUP_SolveEquationGroup(L,U,pi,b);
		for(int j=0;j<n;++j)
			InvA[j][i] = x[j][0];
	}
	return InvA;
}

//normalization of a vector(Matrix) a's element
//input:vector(Matrix) a
//output:a vector(Matrix) after normalization
void normalization(Matrix& a)
{
	int n = a.row;
	double len = 0;
	for(int i=0;i<n;++i)
		len += a[i][0]*a[i][0];
	len = sqrt(len);
	for(int i=0;i<n;++i)
		a[i][0] /= len;
	return;
}

//get the error of two vector(Matrix) a and b
//input:vector(Matrix) a and b
//output:a double number of the error
double getErr(const Matrix& a,const Matrix& b)
{
	double res = 0;
	for(int i=0,n=a.row;i<n;++i)
		res += (fabs(a[i][0])-fabs(b[i][0]))
		*(fabs(a[i][0])-fabs(b[i][0]));
	return sqrt(res);
}

//inverse power method or named inverse iteration method
//input:init matrix A , an eigenValue of A
//output:a eigenVector responding to a
Matrix getEigenVector(Matrix A,const double lambda)
{
	int n = A.row;
	Matrix eigenVector(n,1);
	//create a random vector
	for(int i=0;i<n;++i)
		eigenVector[i][0] = rand()%100/10.0;

	Matrix I(n,n);
	for(int i=0;i<n;++i)
		I[i][i] = 1;

	A = A - I * lambda;
	
	Matrix L(n,n);
	Matrix U(n,n);
	Matrix P(n,n);
	int *pi = LUP_Decomposition(A,L,U,P);
	
	normalization(eigenVector);
	Matrix pre(n,1);
	double Err = 0;
	do
	{
		pre = eigenVector;
		eigenVector = LUP_SolveEquationGroup(L,U,pi,pre);
		normalization(eigenVector);
	}while((Err = getErr(pre,eigenVector)) >= eps);

	return eigenVector;
}

//initialize the init matrix A
void init(Matrix*& A)
{
	ifstream fin("input.txt");
	fin>>n>>m;
	A = new Matrix(n,m);
	for(int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			fin>>A->mat[i][j];
		}
	}
	fin.close();

	cout<<"init A = "<<endl;
	A->print();

	for(int i=0;i<m;++i)
	{
		double tmp = 0;
		for(int j=0;j<n;++j)
			tmp += A->mat[j][i];
		tmp /= n;
		for(int j=0;j<n;++j)
			A->mat[j][i] -= tmp;
		//A->print();
		/*
		tmp = 0;
		for(int j=0;j<n;++j)
			tmp += A->mat[j][i]*A->mat[j][i];
		tmp /= n;
		tmp = sqrt(tmp);
		for(int j=0;j<n;++j)
			A->mat[j][i] /= tmp;
		//A->print();
		*/
	}
	cout<<"after process A = "<<endl;
	A->print();
}

//benchmark
void test();
void test_LUP();
void testInv();

//main logic of PCA algorithm
Matrix PCA(const Matrix* const A)
{
	Matrix Cov(m,m);//covariance matrix
	double* avg = new double[m];//average value of each dimention
	vector < Matrix > col;
	for(int i=0;i<m;++i)
	{
		Matrix cur(n,1);
		double tmp = 0;
		for(int j=0;j<n;++j)
			tmp += A->mat[j][i] , cur[j][0] = A->mat[j][i];
		avg[i] = tmp/n;
		col.push_back(cur);
		//cout<<avg[i]<<" ";
	}
	//cout<<endl;

	for(int i=0;i<m;++i)
	{
		for(int j=0;j<m;++j)
		{
			Cov[i][j] = ((col[i] - avg[i]).transpose()*(col[j] - avg[j]))[0][0] / (n-1);

		}
	}

	//cout<<"Cov = "<<endl;
	//Cov.print();

	delete[] avg;
	eigenValue eigenval;
	vector <double> val = eigenval.getEigenValues(Cov);

	//find top k eigenvalue and get the eigenvector responding to them
	sort(val.begin(),val.end(),greater<double>());
	Matrix eigenMatrix(m,k);
	for(int i=0;i<k;++i)
	{
		Matrix x = getEigenVector(Cov,val[i]);
		x.print();
		for(int j=0;j<m;++j)
		{
			eigenMatrix[j][i] = x[j][0];
		}
	}

	//calc the final matrix after PCA's process 
	Matrix Final(n,k);
	Final = *(A) * eigenMatrix;
	//A->print();
	//eigenMatrix.print();
	cout<<"Fianl = "<<endl;
	Final.print();
	return Final;
}

int _tmain(int argc, _TCHAR* argv[])
{
	//testInv();
	//test_LUP();
	//test();
	Matrix* A = NULL;
	init(A);
	printf("please input how many dimention will you compress from %d dimention:\n",m);
	scanf_s("%d",&k);
	if(k > m)
	{
		printf("PCA is used to reduce dimention!");
		delete A;
		return 0;
	}
	Matrix Fianl = PCA(A);
	delete A;
	return 0;
}

void test()
{
	eigenValue val;
	
	//double data[3][3] = {
	//	5, -3, 2,
	//	6, -4, 4,
	//	4, -4, 5
	//};
	
	double data[3][3] = {
		0,11,-5,
        -2,17,-7,
        -4,26,-10
	};
	Matrix mat(3,3);
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j)
			mat[i][j] = data[i][j];
	std::vector <double> res = val.getEigenValues(mat);
	for(auto i = res.begin();i!=res.end();i++)
	{
		cout<<"eigenValue = "<<*i<<endl;
		Matrix x = getEigenVector(mat,*i);
		cout<<"eigenVector = "<<endl;
		x.print();
		cout<<endl;
	}
}

void test_LUP()
{
	const int maxn = 4;
	FILE* fin;
	freopen_s(&fin,"output.txt","w",stdout);
	double data[maxn][maxn] = {
		2, 0, 2, 0.6,
		3, 3, 4, -2,
		5, 5, 4, 2,
		-1, -2, 3.4, -1
	};
	Matrix mat(maxn,maxn);
	for(int i=0;i<maxn;++i)
		for(int j=0;j<maxn;++j)
			mat[i][j] = data[i][j];
	Matrix L(maxn,maxn);
	Matrix U(maxn,maxn);
	Matrix P(maxn,maxn);
	LUP_Decomposition(mat,L,U,P);
	cout<<"mat = "<<endl;
	mat.print();
	cout<<endl;
	cout<<"L = "<<endl;
	L.print();
	cout<<endl;
	cout<<"U = "<<endl;
	U.print();
	cout<<endl;
	cout<<"P = "<<endl;
	P.print();
	cout<<endl;
	cout<<"L*U = "<<endl;
	(L*U).print();
	cout<<endl;
	cout<<"P*A = "<<endl;
	(P*mat).print();
	cout<<endl;
	fclose(fin);
}

void testInv()
{
	const int maxn = 3;
	double data[maxn][maxn] = {
		5, -3, 2,
		6, -4, 4,
		4, -4, 5
	};
	Matrix mat(maxn,maxn);
	for(int i=0;i<maxn;++i)
		for(int j=0;j<maxn;++j)
			mat[i][j] = data[i][j];
	Matrix L(maxn,maxn);
	Matrix U(maxn,maxn);
	Matrix P(maxn,maxn);
	int *pi = LUP_Decomposition(mat,L,U,P);
	for(int i=0;i<maxn;++i)
	{
		Matrix b(maxn,1);
		for(int j=0;j<maxn;++j)
			b[j][0] = i == j?1:0;
		Matrix x = LUP_SolveEquationGroup(L,U,pi,b);
		x.print();
		cout<<endl;
	}
}