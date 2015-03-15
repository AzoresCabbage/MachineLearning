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

int* LUP_Decomposition(Matrix a,Matrix& L,Matrix& U,Matrix& P)
{
	int n = a.row;
	int *pi = new int[n];
	for(int i=0;i<n;++i)
		pi[i] = i;
	for(int k=0;k<n;++k)
	{
		double p = 0;
		int kk = 0;
		for(int i=k;i<n;++i)
		{
			if(fabs(a.mat[i][k]) > p)
				p = fabs(a.mat[i][k]) , kk=i;
		}
		if(p == 0)
		{
			printf("sigular matrix!\n");
			delete[] pi;
			exit(-1);
		}
		swap(pi[k],pi[kk]);
		for(int i=0;i<n;++i)
			swap(a.mat[k][i],a.mat[kk][i]);
		for(int i=k+1;i<n;++i)
		{
			a.mat[i][k] /= a.mat[k][k];
			for(int j=k+1;j<n;++j)
				a.mat[i][j] -= a.mat[i][k]*a.mat[k][j];
		}
	}
	for(int i=0;i<n;++i)
		P.mat[i][pi[i]] = 1;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
		{
			if(i > j)
				L.mat[i][j] = a.mat[i][j];
			else
				U.mat[i][j] = a.mat[i][j];
		}
		L.mat[i][i] = 1;
	}
	return pi;
	//delete[] pi;
}

Matrix LUP_SolveEquationGroup(Matrix& L,Matrix& U,int* pi,Matrix& b)
{
	int n = L.row;
	Matrix x(n,1);
	Matrix y(n,1);
	for(int i=0;i<n;++i)
	{
		y.mat[i][0] = b.mat[pi[i]][0];
		for(int j=0;j<i;++j)
		{
			y.mat[i][0] -= L.mat[i][j]*y.mat[j][0];
		}
	}
	for(int i=n-1;i>=0;--i)
	{
		x.mat[i][0] = y.mat[i][0];
		for(int j=n-1;j>i;--j)
			x.mat[i][0] -= U.mat[i][j]*x.mat[j][0];
		x.mat[i][0] /= U.mat[i][i];
	}
	return x;
}

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
			b.mat[j][0] = i == j?1:0;
		Matrix x = LUP_SolveEquationGroup(L,U,pi,b);
		for(int j=0;j<n;++j)
			InvA.mat[j][i] = x.mat[j][0];
	}
	return InvA;
}

void constraint(Matrix& a,int n)
{
	double len = 0;
	for(int i=0;i<n;++i)
		len += a.mat[i][0]*a.mat[i][0];
	len = sqrt(len);
	for(int i=0;i<n;++i)
		a.mat[i][0] /= len;
	return;
	double Ck = -DBL_MAX_10_EXP;
	for(int i=0;i<n;++i)
		Ck = max(Ck,a.mat[i][0]);
	for(int i=0;i<n;++i)
		a.mat[i][0] /= Ck;
}

double getErr(const Matrix& a,const Matrix& b)
{
	double res = 0;
	for(int i=0,n=a.row;i<n;++i)
		res += (fabs(a.mat[i][0])-fabs(b.mat[i][0]))
		*(fabs(a.mat[i][0])-fabs(b.mat[i][0]));
	return sqrt(res);
}

Matrix getEigenVector(Matrix A,double lambda)
{
	int n = A.row;
	Matrix eigenVector(n,1);
	for(int i=0;i<n;++i)
		eigenVector.mat[i][0] = rand()%100/10.0;

	Matrix I(n,n);
	for(int i=0;i<n;++i)
		I.mat[i][i] = 1;

	A = A - I * lambda;
	
	Matrix L(n,n);
	Matrix U(n,n);
	Matrix P(n,n);
	int *pi = LUP_Decomposition(A,L,U,P);
	
	constraint(eigenVector,n);
	Matrix pre(n,1);
	double Err = 0;
	do
	{
		pre = eigenVector;
		eigenVector = LUP_SolveEquationGroup(L,U,pi,pre);
		constraint(eigenVector,n);
	}while((Err = getErr(pre,eigenVector)) >= eps);

	return eigenVector;
}

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

/*
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
			mat.mat[i][j] = data[i][j];
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
*/

/*
void test_LUP()
{
#define maxn 4
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
			mat.mat[i][j] = data[i][j];
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
*/

/*
void testInv()
{
#define maxn 3
	double data[maxn][maxn] = {
		5, -3, 2,
		6, -4, 4,
		4, -4, 5
	};
	Matrix mat(maxn,maxn);
	for(int i=0;i<maxn;++i)
		for(int j=0;j<maxn;++j)
			mat.mat[i][j] = data[i][j];
	Matrix L(maxn,maxn);
	Matrix U(maxn,maxn);
	Matrix P(maxn,maxn);
	int *pi = LUP_Decomposition(mat,L,U,P);
	for(int i=0;i<maxn;++i)
	{
		Matrix b(maxn,1);
		for(int j=0;j<maxn;++j)
			b.mat[j][0] = i == j?1:0;
		Matrix x = LUP_SolveEquationGroup(L,U,pi,b);
		x.print();
		cout<<endl;
	}
}
*/

Matrix PCA(Matrix* A)
{
	Matrix Cov(m,m);
	double* avg = new double[m];
	vector < Matrix > col;
	for(int i=0;i<m;++i)
	{
		Matrix cur(n,1);
		double tmp = 0;
		for(int j=0;j<n;++j)
			tmp += A->mat[j][i] , cur.mat[j][0] = A->mat[j][i];
		avg[i] = tmp/n;
		col.push_back(cur);
		//cout<<avg[i]<<" ";
	}
	//cout<<endl;

	for(int i=0;i<m;++i)
	{

		for(int j=0;j<m;++j)
		{
			Cov.mat[i][j] = ((col[i] - avg[i]).transpose()*(col[j] - avg[j])).mat[0][0] / (n-1);

		}
	}

	//cout<<"Cov = "<<endl;
	//Cov.print();

	delete[] avg;
	eigenValue eigenval;
	vector <double> val = eigenval.getEigenValues(Cov);
	sort(val.begin(),val.end(),greater<double>());
	Matrix eigenMatrix(m,k);
	for(int i=0;i<k;++i)
	{
		Matrix x = getEigenVector(Cov,val[i]);
		x.print();
		for(int j=0;j<m;++j)
		{
			eigenMatrix.mat[j][i] = x.mat[j][0];
		}
	}
	Matrix Final(n,k);
	Final = *(A) * eigenMatrix;
	//A->print();
	//eigenMatrix.print();
	//Final.print();
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