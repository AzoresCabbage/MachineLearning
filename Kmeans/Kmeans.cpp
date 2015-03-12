/*
 * Name:Kmeans.cpp
 * Author:wyj
 * Date:2015.03.12
 * Version:1.0
 * Description:implemention of Kmeans algorithm
 */

#include <bits/stdc++.h>
using namespace std;

#define eps 1e-9
#define TYPE double
#define INF 1000000000

int n,dim,k;//n:point num ; dim:point dimention ; k:find k center

vector < vector <TYPE> > input;//sample set input
vector < int > c;//center of point i
vector < vector <TYPE> > u;//position of center i

template <class T>
T SQR(T x){return x*x;}

TYPE sqrSum(vector <TYPE> x,vector <TYPE> y)
{
    double res = 0;
    for(int j=0;j<dim;++j)
        res += SQR(x[j]-y[j]);
    return res;
}

TYPE eval()
{
    TYPE res = 0;
    for(int i=0;i<n;++i)
        res += sqrSum(input[i],u[c[i]]);
    return res;
}

int main()
{
    cout<<"please input how many center will you find in data set?"<<endl;
    cin>>k;

    freopen("input.txt","r",stdin);
    cin>>n>>dim;
    c.resize(n);
    input.resize(n);
    for(int i=0;i<n;++i)
    {
        input[i].resize(dim);
        for(int j=0;j<dim;++j)
            cin>>input[i][j];
    }
    fclose(stdin);

    srand((unsigned int)time(NULL));
    u.resize(k);
    for(int i=0;i<k;++i)
    {
        u[i].resize(dim);
        int idx = rand()%n;
        u[i] = input[idx];
    }
    
    TYPE pre = 0,curErr = 0;
    do
    {
        for(int i=0;i<n;++i)
        {
            TYPE mn = INF;
            int pos = c[i];
            for(int j=0;j<k;++j)
            {
                TYPE tmp = sqrSum(input[i],u[j]);
                if(tmp < mn)
                {
                    mn = tmp;
                    pos = j;
                }
            }
            c[i] = pos;
        }

        for(int j=0;j<k;++j)
        {
            vector <TYPE> tmp(dim,0);
            int tot = 0;
            for(int i=0;i<n;++i)
            {
                if(c[i] != j) continue;
                tot++;
                for(int d=0;d<dim;++d)
                    tmp[d] += input[i][d];
            }
            if(tot == 0) continue;
            for(int d=0;d<dim;++d)
                u[j][d] = 1.0*tmp[d]/tot;
        }
        
        TYPE curVal = eval();
        curErr = fabs(curVal-pre);
        pre = curVal;
    }while(curErr > eps);

    for(int i=0;i<k;++i)
    {
        cout<<"center "<<i<<":"<<endl;
        for(int j=0;j<dim;++j)
            cout<<u[i][j]<<" ";
        cout<<endl;
        cout<<endl;
    }

    return 0;
}
