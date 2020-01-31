#include <iostream>
#include <fstream>

#include <cstdlib>
#include <time.h>
#include <math.h>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric;

typedef ublas::vector<float> vec;
typedef ublas::matrix<float> mat;
typedef ublas::zero_matrix<float> zeros;

mat Wprime(mat r);
mat Phi(mat r);

// Generate a matrix with random numbers (size N1 x N2)
mat random_matrix(int N1, int N2);
// Repeat the A matrix mxn times
mat repmat(mat A, int m, int n);
// Return a matrix with each value of x raised to the power of y
mat pow_matrix(mat x, float y);
// Generate a diagonal matrix
mat diagonal_matrix(vec v);
// Return a matrix with each value of x rounded downward
mat floor(mat a);
// Return the sum of items in a
mat sum_matrix(mat a);
// Save results in a .txt file
void write_file(int n, mat x1, mat x2, mat v1, mat v2);

mat Wprime(mat r)
{
    double a = 1.001;
    double b = 3.;

    int N1 = r.size1();
    int N2 = r.size2();

    mat Wp(N1,N2);

    Wp = pow_matrix(r,a-1) - pow_matrix(r,b-1);

    return Wp;
}

mat Phi(mat r)
{
    int i,j;

    int N1 = r.size1();
    int N2 = r.size2();

    mat P(N1,N2);

    for(i=0;i<N1;i++)
        for(j=0;j<N2;j++)
            P(i,j) = 1./sqrt(1+pow(r(i,j),2));
     
    return P;
}

mat random_matrix(int N1, int N2)
{
    typedef boost::mt19937 RNGType;
    RNGType rng(clock());
    boost::uniform_real<> random_float(0,1);
    boost::variate_generator< RNGType, boost::uniform_real<> > rand(rng, random_float);

    int i,j;
    mat rand_m(N1,N2);

    for(i=0;i<N1;i++)
        for(j=0;j<N2;j++)
            rand_m(i,j) = rand();

    return rand_m;
}

mat repmat(mat A, int m, int n)
{
    int i,j;

    int N1 = A.size1();
    int N2 = A.size2();

    mat A2(m*N1,n*N2);

    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
            ublas::project(A2,ublas::range(i*N1,(i+1)*N1),ublas::range(j*N2,(j+1)*N2)) = A;

    return A2;
}

mat pow_matrix(mat x, float y)
{
    int i,j;
    mat result(x.size1(), x.size2());

    for(i=0;i<int(x.size1());i++)
        for(j=0;j<int(x.size2());j++)
            result(i,j) = pow(x(i,j),y);

    return result;
}

mat diagonal_matrix(vec v)
{
    int i;
    int N = v.size();
    mat D = zeros(N,N);

    for(i=0;i<N;i++)
        D(i,i) = v(i);

    return D;
}

mat floor(mat A)
{
    int i,j;
    mat result(A.size1(), A.size2());

    for(i=0;i<int(A.size1());i++)
        for(j=0;j<int(A.size2());j++)
            result(i,j) = int(A(i,j));

    return result;
}

mat sum_matrix(mat A)
{
    int i,j;
    int N1 = A.size1();
    int N2 = A.size2();
    mat s = zeros(N1,1);

    for(i=0;i<N1;i++)
        for(j=0;j<N2;j++)
            s(i,0) += A(i,j);

    return s;
}

void write_file(int n, mat x1, mat x2, mat v1, mat v2)
{
    int i;
    string const filename("data/" + to_string(n) + ".txt");
    ofstream flux(filename.c_str());

    if(flux)    
        for(i=0;i<x1.size1();i++)
            flux << x1(i,0) << " ; " << x2(i,0) << " ; " 
                 << v1(i,0) << " ; " << v2(i,0) << endl;
    else
        cout << "ERROR: CAN'T OPEN'T THE FILE!" << endl;
}
