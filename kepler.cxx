#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;
void symplectic(double* p, double* q, const double& dt, double& H);

int main(){
const double e = 0.6;
const double pi = acos(-1.0);
const int dim = 2;
const double tEnd = 20.0*pi;
const double dt = 0.05;
double H;
double q[dim];
double p[dim];

double t = 0.0;
q[0] = 1.0-e;
q[1] = 0.0;
p[0] = 0.0;
p[1] = sqrt((1.0+e)/(1.0-e));

H = 0.5*(p[0]*p[0] + p[1]*p[1]) - 1.0/sqrt(q[0]*q[0] + q[1]*q[1]);

ofstream out("data");
out << t << "\t" << q[0] << "\t" << q[1] << "\t" << H << endl;

while (t<tEnd){
	symplectic(p, q, dt, H);
	out << t << "\t" << q[0] << "\t" << q[1] << "\t" << H << endl;
	t += dt;
}
out.close();
return 0;
}

void symplectic(double* p, double* q, const double& dt, double& H){
	double x = sqrt(q[0]*q[0] + q[1]*q[1]);
	for(int i=0; i<2; i++){
		p[i] -= dt*q[i]/pow(x,3);
		q[i] += dt*p[i];
	}
	H = 0.5*(p[0]*p[0] + p[1]*p[1]) - 1.0/sqrt(q[0]*q[0] + q[1]*q[1]);

}
