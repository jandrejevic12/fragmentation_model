#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

//#include "omp.h"

const double b=1.8,lb=log(b);				// base to use for logarithmically sized grid
const double dx=0.0125,ei=-50+dx,ef=dx+1e-5;	// range of exponents to use, and increment size (grid spacing)
const int N=int((ef-ei)/dx);				// number of grid lengths along each dimension (grid points = N+1)
const long int T=80100;						// number of timesteps
const long int out=100;						// how often to output to file
const double dt=0.001;						// step size
const double lam=0.5;						// exponent for overall rate of fragmentation
const double a=1.;							// shape parameter for conditional probability of breakup

// Function to return the wall clock time.
//double wtime() {return omp_get_wtime();}

// Tell the compiler about the existence of the required LAPACK functions.
extern "C" {
    void dtrtrs_(char *uplo,char *trans,char *diag,int *n,int *nrhs,
    			 double *a,int *lda,double *b,int *ldb,int *info);
}

// overall rate of fragmentation
double r(double x,double lam) {
	return pow(x,lam);
}

// conditional probability that x results from breakup of y
double f(double x,double y,double a) {
	return 0.5/y*(a+2)*pow(x/y,0.5*a-1);
}

// average facet area (s)
double xavg(double *c) {
	double s=0.5*((*c)*pow(b,ei)+c[N]*pow(b,ei+N*dx));
	for(int i=1;i<N;i++) s+=(c[i])*pow(b,(ei+i*dx));
	return 1./(s*dx*lb);
}

// total area (should be conserved and equal to 1)
double xtot(double *c) {
	double s=0.5*((*c)*pow(b,2*ei)+c[N]*pow(b,2*(ei+N*dx)));
	for(int i=1;i<N;i++) s+=(c[i])*pow(b,2*(ei+i*dx));
	return s*dx*lb;
}

// print the linear system of equations
void print_matrix(double *A) {
	for(int i=0;i<N+1;i++) {
		for(int j=0;j<N+1;j++) {
			printf("%g ",A[j*(N+1)+i]);
		}
		printf("\n");
	}
}

// assemble matrix for linear system to solve, in column-major order for LAPACK
void assemble_matrix(double *A) {
	int i,j;
	double x,y;
	for(i=0;i<N;i++) {
		x=pow(b,ei+i*dx);
		for(j=i+1;j<N;j++) {
			y=pow(b,ei+j*dx);
			A[j*(N+1)+i]=-dx*dt*r(y,lam)*f(x,y,a)*y*lb;
		}
		y=pow(b,ei+N*dx);
		A[N*(N+1)+i]=-0.5*dx*dt*r(x,lam)*f(x,y,a)*y*lb;
	}
	for(i=0;i<N;i++) {
		x=pow(b,ei+i*dx);
		A[i*(N+1)+i]=1+dt*r(x,lam)*(1-0.5*dx*f(x,x,a)*x*lb);
	}
	x=pow(b,ei+N*dx);
	// last grid point at x>1 has ommitted depletion term
	A[(N+1)*(N+1)-1]=1;//+dt*r(x,lam);
}

// Solve the fragmentation rate equation.
int main () {
	// settings for LAPACK linear solve
	char uplo='U',trans='N',diag='N';
	int n=N+1,nrhs=1,info;

	// create matrix and vectors, and zero out entries
	double *A=new double[(N+1)*(N+3)],*c1=A+(N+1)*(N+1),*c2=c1+(N+1);
	for(double *Ap=A;Ap<A+(N+1)*(N+3);Ap++) (*Ap)=0.;

	// output file
	FILE *fp;
	fp=fopen("frag_data.txt","w");
	if(fp==NULL) {
    	fprintf(stderr,"Can't open file.\n");
    	exit(1);
	}
	// output discrete x to file, and lam and a
	for(int i=0;i<N+1;i++) fprintf(fp,"%g ",pow(b,ei+i*dx));
	fprintf(fp,"%g %g\n",lam,a);

	// assemble matrix
	assemble_matrix(A);

	// print matrix
	//print_matrix(A);
    
    // specify initial condition, in both solution vectors
    //c1[N]=c2[N]=2./dx/lb*pow(b,-2*(ei+N*dx));
    c1[N-1]=c2[N-1]=1./dx/lb*pow(b,-2*(ei+(N-1)*dx));
    //c1[1810]=c2[1810]=1./dx/lb*pow(b,-2*(ei+1810*dx))*99./100.;
    //c1[1591]=c2[1591]=1./dx/lb*pow(b,-2*(ei+1591*dx));

    // variables for average facet area and total area
    double s,m;
    s=xavg(c1); m=xtot(c1);
    printf("%g %g\n",s,m);
    fprintf(fp,"0 ");
    for(int i=0;i<N+1;i++) fprintf(fp,"%g ",c1[i]);
    fprintf(fp,"%g\n",s);
	
	// carry out the first step; first order backward Euler, solution computed in-place
	dtrtrs_(&uplo,&trans,&diag,&n,&nrhs,A,&n,c1,&n,&info);
	if(1%out==0) {
		s=xavg(c1); m=xtot(c1);
		printf("%g %g\n",s,m);
		fprintf(fp,"%g ",dt);
		for(int i=0;i<N+1;i++) fprintf(fp,"%g ",c1[i]);
		fprintf(fp,"%g\n",s);
    }
    // solve remaining steps with second order accuracy, updating matrix
    for(int i=0;i<N+1;i++) A[i*(N+1)+i]+=0.5;
    double *tmp;
    for(int t=2;t<T;t++) {
    	for(double *cp1=c1,*cp2=c2;cp1<c1+(N+1);cp1++,cp2++) *cp2=2*(*cp1)-0.5*(*cp2);
    	dtrtrs_(&uplo,&trans,&diag,&n,&nrhs,A,&n,c2,&n,&info);
    	tmp=c2;c2=c1;c1=tmp;
    	if(t%out==0) {
    		s=xavg(c1); m=xtot(c1);
    		printf("%g %g\n",s,m);
    		fprintf(fp,"%g ",t*dt);
    		for(int i=0;i<N+1;i++) fprintf(fp,"%g ",c1[i]);
    		fprintf(fp,"%g\n",s);
    	}
    }
	fclose(fp);
	delete [] A;

}
