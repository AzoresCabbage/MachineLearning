/*
 * Name: SMO Implement
 * Author£ºwyj
 * Date£º2015.03.10
 * Version£º1.0
 *
 * Reference£º"Sequential minimal Optimization A Fast Algorithm for Training Support Vector Machines"  pseudo code
 * 
 * none liner SVM's u:u[i]=sigma(j=1->N) y[j]*alpha[j]*K(x[j],x[i])-b[i]
 *
 * Training Input:Data set
 * Training Output:None
 *
 * Runing Input:One Sample Input
 * Runing Output:the Result
 */

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <fstream>
#include <time.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

#define scanf_s scanf
#define tol 0.001
#define eps 0.001
#define twoSigmaSqr 2
#define LINER 0
#define NONLINER 1

const double C = 1;

vector < vector <double> > input;//training set input
vector <int> target;//desired output
vector <double> alpha;//Lagrange multiplier
vector <double> w;
vector <double> error_cache;
vector < vector<double> > precomputeDot;//prepare compute the dot between samples
int n;//training set size
int dim;//dimention
double b;//threshold

template <class T>
T SQR(const T& a) {return a*a;}

void see(const vector <double>& v)
{
    for(int i=0,j=v.size();i<j;++i)
        cout<<v[i]<<" ";
    cout<<b<<endl;
}

//input:vector v1,v2
//target:v1*v2^T
//dot computation
double dot(const vector <double>& v1,const vector <double>& v2)
{
    double res = 0;
    for(int i=0;i<dim;++i)
        res += v1[i]*v2[i];
    return res;
}

//input:index of vector x1,x2 and flag kernelType to judge which method to use
//target:kernel function value
//measures the similarity or distance between x1 && x2
//Liner kernel && Radial Basis Function kernel(rbf kernel)
double kernel(int i1,int i2,int kernelType)
{
    //liner kernel
    if(kernelType == LINER)
        return precomputeDot[i1][i2];
    //rbf kernel
    else{
        double res = SQR(precomputeDot[i1][i1]) + SQR(precomputeDot[i2][i2]) - 2*precomputeDot[i1][i2];
        return exp(- res / twoSigmaSqr);
    }
}

//input:vector x,y and flag kernelType to judge which method to use
//target:kernel function value
//Liner kernel && Radial Basis Function kernel(rbf kernel)
double kernel(const vector <double>& x,const vector <double>& y,int kernelType)
{
    double res = 0;
    if(kernelType == LINER)
    {
        res = dot(x,y);
    }
    else
    {
        double xx = dot(x,x);
        double yy = dot(y,y);
        double xy = dot(x,y);
        res = exp(- (SQR(xx) + SQR(yy) - 2*xy) / twoSigmaSqr);
    }
    return res;
}

//u = WX - b
double liner_learned_func(int idx)
{
    double res = dot(w,input[idx]);
    return res - b;
}

//u = sigma(alpha[i]*target[i]*Kij) - b
double nonliner_learned_func(int idx)
{
    double res = 0;
    for(int i=0;i<n;++i)
    {
        res += alpha[i]*target[i]*kernel(idx,i,NONLINER);
    }
    return res - b;
}

//input:the index of sample and flag kernelType to judge which method to use
double learned_func(int idx,int kernelType)
{
    if(kernelType == LINER)
        return liner_learned_func(idx);
    else
        return nonliner_learned_func(idx);
}

//input:index of two candidate Lagrange multipliers and flag kernelType
//target:if this pair is ok return true,else return false
bool takeStep(int i1,int i2,int kernelType)
{
    if(i1 == i2) return false;
    double alpha1 = alpha[i1];//old alpha
    double alpha2 = alpha[i2];
    double y1 = target[i1];
    double y2 = target[i2];
    double E1 = 0;
    double E2 = 0;
    double s = y1*y2;
    double L,H;
    double a1,a2;//new alpha

    if(0 < alpha1 && alpha1 < C)
        E1 = error_cache[i1];
    else
        E1 = learned_func(i1,kernelType) - y1;

    if(0 < alpha2 && alpha2 < C)
        E2 = error_cache[i2];
    else
        E2 = learned_func(i2,kernelType) - y2;

    if(s < 0)
    {
        L = max(0.0,alpha2-alpha1);
        H = min(C,C+alpha2-alpha1);
    }
    else
    {
        L = max(0.0,alpha1+alpha2-C);
        H = min(C,alpha1+alpha2);
    }
    if(L == H)  return 0;
    double k11 = kernel(i1,i1,kernelType);
    double k12 = kernel(i1,i2,kernelType);
    double k22 = kernel(i2,i2,kernelType);
    double eta = k11+k22-2*k12;
    if(eta > 0)
    {
        a2 = alpha2 + y2*(E1-E2)/eta;
        if(a2 < L) a2 = L;
        else if(a2 > H) a2 = H;
    }
    else
    {
        double f1 = y1*(E1+b) - alpha1*k11 - s*alpha2*k12;
        double f2 = y2*(E2+b) - alpha2*k22 - s*alpha1*k12;
        double L1 = alpha1 + s*(alpha2-L);
        double H1 = alpha1 + s*(alpha2-H);
        double Lobj = L1*f1 + L*f2 + 0.5*L1*L1*k11 + 0.5*L*L*k22 + s*L*L1*k12;
        double Hobj = H1*f1 + H*f2 + 0.5*H1*H1*k11 + 0.5*H*H*k22 + s*H*H1*k12;

        if(Lobj < Hobj-eps)
            a2 = L;
        else if(Lobj > Hobj+eps)
            a2 = H;
        else
            a2 = alpha2;
    }
    //incresment nearly to be 0
    if(fabs(a2-alpha2) < eps*(a2+alpha2+eps))
    {
        //puts("incresment nearly to be 0");
        return 0;
    }
    a1 = alpha1+s*(alpha2-a2);

    //copy from smo.java,but why? paper say a1,not a1^cliped
    /*if(a1 < 0)
    {
        a2 += s*a1;
        a1 = 0;
    }
    else if(a1 > C)
    {
        a2 += s*(a1-C);
        a1 = C;
    }*/

    double preb = 0;
    double b1 = E1 + y1*(a1-alpha1)*k11 + y2*(a2-alpha2)*k12 + b;
    double b2 = E2 + y1*(a1-alpha1)*k12 + y2*(a2-alpha2)*k22 + b;
    if(0 < a1 && a1 < C)
        preb = b , b = b1;
    else if(0 < a2 && a2 < C)
        preb = b , b = b2;
    else
        preb = b , b = 0.5*(b1+b2);
    
    double coef1 = y1*(a1-alpha1);
    double coef2 = y2*(a2-alpha2);

    
    if(kernelType == LINER)
        for(int i=0;i<dim;++i)
            w[i] += coef1*input[i1][i] + coef2*input[i2][i];

    for(int i=0;i<n;++i)
        if(0 < alpha[i] && alpha[i] < C)
            error_cache[i] += coef1*kernel(i1,i,kernelType)
                + coef2*kernel(i2,i,kernelType) + preb - b;
    error_cache[i1] = error_cache[i2] = 0;

    alpha[i1] = a1;
    alpha[i2] = a2;

    return 1;
}

//fix one Langrange multiplier
//input:index of the second Lagrange multiplier and flag kernelType
//target:if i2 can find another Lagrange multiplier to optimize the result,then return true,else return false
int examineExample(int i2,int kernelType)
{
    double y2 = target[i2];
    double alpha2 = alpha[i2];
    double E2 = 0;

    if(0 < alpha2 && alpha2 < C)
        E2 = error_cache[i2];
    else
        E2 = learned_func(i2,kernelType) - y2;

    double r2 = E2*y2;//y2*y2 = 1
    if((r2 < -tol && alpha2 < C) || (r2 > tol && alpha2 > 0))
    {
        int i1 = -1;
        double mx = 0;

        //second heuristic
        for(int i=0;i<n;++i)
        {
            double tmp = fabs(E2 - error_cache[i]);
            if( 0 < alpha[i] && alpha[i] < C && mx < tmp )
                mx = tmp,i1 = i;
        }
        if(i1>=0)
            if(takeStep(i1,i2,kernelType))
                return 1;

        //unusual circumstances
        //use random postition to avoid bias towards beginning or ending of dataset
        int st = rand() % n;
        for(int i=st;i<st+n;++i)
        {
            i1 = i;
            if(i1 >= n) i1 -= n;
            if(alpha[i1] == 0 || alpha[i1] == C) continue;
            if(takeStep(i1,i2,kernelType)) return 1;
        }
        st = rand() % n;
        for(int i=st;i<st+n;++i)
        {
            i1 = i;
            if(i1 >= n) i1 -= n;
            if(takeStep(i1,i2,kernelType)) return 1;
        }
    }

    return 0;
}

//main logic
void train(int kernelType)
{
    srand((unsigned int)time(NULL));
    //logic from the paper in section "Heuristics for Choosing Which Multipliers to Optimize"
    int numChanged = 0;
    int examineAll = 1;
    int iteration = 0;
    int maxIter = 100000;
    while(numChanged > 0 || examineAll)
    {
        numChanged = 0;
        if(examineAll)
        {
            for(int i=0;i<n;++i)
            {
                numChanged += examineExample(i,kernelType);
            }
        }
        else
        {
            for(int i=0;i<n;++i)
            {
                if(0 < alpha[i] && alpha[i] < C)
                    numChanged += examineExample(i,kernelType);
            }
        }
        if(examineAll == 1)
            examineAll = 0;
        else if(numChanged == 0)
            examineAll = 1;
        printf("Iteration = %d\n",iteration++);
        if(maxIter < iteration) break;
    }
    if(kernelType == LINER)
        see(w);
    else
        see(alpha);
}

//read training set data from file and initialize variables
bool init(string fin,string fout)
{
	ifstream in(fin.c_str());
    double Dvalue;
    in>>n>>dim;
    for(int i=1;i<=n;++i)
    {
        vector <double> tmp;
        for(int j=1;j<=dim;++j)
        {
			in>>Dvalue;
            tmp.push_back(Dvalue);
        }
        input.push_back(tmp);
    }
	in.close();

    cout<<"read in.txt done!"<<endl;

	in.open(fout.c_str());
    int targetRow;
	in>>targetRow;
    
    if(n != targetRow)
        return false;
    
    int Ivalue;
    for(int i=1;i<=n;++i)
    {
		in>>Ivalue;
        target.push_back(Ivalue);
    }
	in.close();

    cout<<"read out.txt done!"<<endl;
    
    precomputeDot.resize(n);
    for(int i=0;i<n;++i)
        precomputeDot[i].resize(n);

    w.resize(dim);
    error_cache.resize(n);
    alpha.resize(n,0);
    
    cout<<"training set size = "<<n<<" dimention = "<<dim<<endl;
	for(int i=0;i<n;++i)
        error_cache[i] = 0-target[i];
    for(int i=0;i<n;++i)
        for(int j=i;j<n;++j)
            precomputeDot[i][j] = precomputeDot[j][i] = dot(input[i],input[j]);
    return true;
}

//input:predict vector x and flag indicate which method to use
//target:return which class belong to for x
int predict(const vector <double>& x,int kernelType)
{
    double res = 0;
    if(kernelType == LINER)
        res = dot(w,x)-b;
    else
    {
        for(int i=0;i<n;++i)
        {
            if(alpha[i] == 0) continue;
            res += target[i]*alpha[i]*kernel(x,input[i],kernelType);
            //cout<<target[i]<<" "<<alpha[i]<<" "<<kernel(x,input[i],kernelType)<<endl;
        }
        res -= b;
    }
    cout<<"res = "<<res<<endl;
    if(fabs(1-res) > fabs(-1-res)) return -1;
    else return 1;
}

void readModel(char* filename,int &kernelType)
{
    ifstream fin(filename);
    fin>>kernelType;
    fin>>n>>dim;
    if(kernelType == LINER)
    {
        for(int i=0;i<n;++i)
            fin >> w[i];
        fin >> b;
    }
    else
    {
        target.resize(n);
        alpha.resize(n);
        input.resize(n);
        for(int i=0;i<n;++i)
            fin >> target[i];
        for(int i=0;i<n;++i)
            fin >> alpha[i];
        for(int i=0;i<n;++i)
        {
            input[i].resize(dim);
            for(int j=0;j<dim;++j)
                fin>>input[i][j];
        }
        fin >> b;
    }
    fin.close();
}

void writeModelToFile(char* filename,int kernelType)
{
    ofstream fout(filename);
    fout<<kernelType<<endl;
    fout<<n<<" "<<dim<<endl;
    if(kernelType == LINER)
    {
        //u = WX - b
        for(int i=0;i<n;++i)
            fout<<w[i]<<" ";
        fout<<b<<endl;
    }
    else
    {
        //u = sigma(y[i]*alpha[i]*kernel(X,X[i])) - b
        for(int i=0;i<n;++i)
            fout<<target[i]<<" ";
        fout<<endl;
        for(int i=0;i<n;++i)
            fout<<alpha[i]<<" ";
        fout<<endl;
        for(int i=0;i<n;++i)
        {
            for(int j=0;j<dim;++j)
                fout<<input[i][j]<<" ";
            fout<<endl;
        }
        fout<<b<<endl;
    }
    fout.close();
}

int main()
{
    int opt;
    printf("1.build new model \n2.classify from old model\n");
    while(true)
    {
        scanf("%d",&opt);
        if(opt == 1 || opt == 2)
            break;
        printf("1.build new model \n2.classify from old model\n");
    }
    if(opt == 1)
    {
        printf("1.liner kernel\n2.nonliner kernel\n");
        while(true)
        {
            scanf("%d",&opt);
            if(opt == 1 || opt == 2)
                break;
            printf("1.liner kernel\n2.nonliner kernel\n");
        }
        cout<<"init..."<<endl;
        if(!init("in.txt","out.txt"))
        {
            printf("Training set's input file and target file not match!\n");
            return 1;
        }
        cout<<"init done!"<<endl;
        cout<<"trainning..."<<endl;
        train(opt==1?LINER:NONLINER);
        cout<<"trainning done!"<<endl;
        char path[500];
        printf("please input model save path:\n");
        scanf("%s",path);
        writeModelToFile(path,opt==1?LINER:NONLINER);
    }
    else
    {
        char path[500];
        printf("please input model path:\n");
        scanf("%s",path);
        int kernelType = LINER;
        readModel(path,kernelType);
        printf("This is a %s model,dimention = %d\n",kernelType == LINER?"liner":"nonliner",dim);
        double val;
        vector <double> x(dim);
        while(true)
        {
            printf("please input predict vector:\n");
            for(int i=0;i<dim;++i)
                scanf_s("%lf",&val) , x[i] = val;
            int res = predict(x,kernelType);
            printf("the prediction is %d\n",res);
        }
    }
    return 0;
}
