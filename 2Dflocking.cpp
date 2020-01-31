#include "2Dflocking.hpp"

using namespace std;
using namespace boost::numeric;

typedef ublas::vector<float> vec;
typedef ublas::matrix<float> mat;
typedef ublas::zero_matrix<float> zeros;

int main()
{
    // Parameter of the code
    int N = 1000;            // Number of particles
    float dt = 0.1;          // Time step
    float T = 10;            // Final time
    int nT = int(T/dt);      // Number of time step
    int sav = 1;             // Saving steps
    int a = 5;               // Size of the box

    mat alpha(N,1,1.);
    float beta = 0.05;

    clock_t t_start,t_stop;
    float compute_time;

    int k;

    cout << "Number of particles = " << N << endl;
    cout << "Final time = " << T << endl;
    cout << "Time step = {0:" << dt << "}" << endl;

    // Initial conditions
    mat x1(N,1);
    x1 = random_matrix(N,1) + mat(N,1,a);
    mat x2(N,1);
    x2 = random_matrix(N,1) + mat(N,1,a);

    mat v1(N,1);
    v1 = mat(N,1,-1) + 2*random_matrix(N,1);
    mat v2(N,1);
    v2 = mat(N,1,-1) + 2*random_matrix(N,1);

    mat dx1(N,1);
    dx1 = zeros(N,1);
    mat dx2(N,1);
    dx2 = zeros(N,1);

    mat dv1(N,1);
    dv1 = zeros(N,1);
    mat dv2(N,1);
    dv2 = zeros(N,1);

    mat x1Diff(N,N);
    mat x2Diff(N,N);
    mat normxDiff(N,N);

    mat v1Diff(N,N);
    mat v2Diff(N,N);

    mat WpX(N,N);

    // Identity matrix
    mat Id(N,N);
    Id = diagonal_matrix(vec(N,1));

    // Evolution: Forward Euler Scheme
    int cnt = 0;
    int cnt_sav = 0;

    t_start = clock();

    for(k=0;k<nT;k++)
    {
        // Compute the evolution of the position
        dx1 = v1;
        dx2 = v2;

        // and velocity
        x1Diff = repmat(x1,1,N) - repmat(ublas::trans(x1),N,1);
        x2Diff = repmat(x2,1,N) - repmat(ublas::trans(x2),N,1);
        normxDiff = pow_matrix((pow_matrix(x1Diff,2) + pow_matrix(x2Diff,2)),0.5);

        v1Diff = repmat(v1,1,N) - repmat(ublas::trans(v1),N,1);
        v2Diff = repmat(v2,1,N) - repmat(ublas::trans(v2),N,1);

        dv1 = element_prod(alpha - beta*(pow_matrix(v1,2)),v1);
        dv2 = element_prod(alpha - beta*(pow_matrix(v2,2)),v2);

        WpX = Wprime(normxDiff);

        normxDiff += Id;
        
        // dv1 = 1./N * sum(x1Diff*WpX/normxDiff)
        dv1 += 1./N * sum_matrix(element_div(element_prod(x1Diff,WpX),normxDiff));
        dv2 += 1./N * sum_matrix(element_div(element_prod(x2Diff,WpX),normxDiff));

        x1 += dt*dx1;
        x2 += dt*dx2;

        x1 -= 2*a*floor(x1/(2*a));
        x2 -= 2*a*floor(x2/(2*a));

        v1 += dt*dv1;
        v2 += dt*dv2;

        // Every sav iterations
        if(cnt%sav == 0)
        {
	    // Save data
            write_file(cnt,x1,x2,v1,v2);
            cnt++;
        }
    }

    t_stop = clock();
    compute_time = ((float)t_stop - t_start) / CLOCKS_PER_SEC;

    cout << "Compute time: " << compute_time << "sec" << endl;

    return 0;
}
