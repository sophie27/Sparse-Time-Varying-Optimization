/* 7 JAN 2020 */
/** Code written by Sophie M. Fosson **/
/** Compressed System Identification via online iterative soft thresholding, Douglas-Rachford splitting (ADMM), Distributed iterative soft thresholding **/
// Time varying parameters evoving as STEPWISE CONSTANT functions

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
using namespace std;
using namespace Eigen;
#define T 1000
#define LAMBDA 1e-2
//#define TAU 8e-3
#define MU 1e-6
#define pi 3.1415926535897932384626433832795028841971


double q, SNR, MAE_couple, MAE, TOT_MAE, thresh, thresh_dist,  np, admitted_run_time, auxdista, TAU;
double alpha0, alpha1, alpha, gap, gap2, fopt, fist, fadmm, fdista, dinreg_ist, dinreg_admm, dinreg_dista, dinreg_ist_theo, dinreg_admm_theo, dinreg_dista_theo;
int N,M,K,numNodes,v,h,i,j,m,n, k,t, run, num_neighbors, measuresNode;

// required functions
double Gaussian_Noise(double std) {
    return std*(sqrt(-2*log( ((double) rand() / (RAND_MAX))  )))*cos(2*pi* ((double) rand() / (RAND_MAX))  );
}

double soft_thresholding(double v,  double thresh) {
	double out;
	if ( abs(v) < thresh )
		out=0;
	else if (v>0)
		out= v-thresh;
	else if (v<0)
		out= v+thresh;
	return out;
}


int main (int argc, char *argv[]) {
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << "  <N> <K> <M> <run> <noise power>" << endl;
        cerr << "Ex: ./dynreg 20 2 15 1 0.01" << endl;
		return EXIT_FAILURE;
	}
    N = atoi (argv[1]);
	K = atoi (argv[2]);
	M = atoi (argv[3]);
	run=atoi(argv[4]);
	np=atof(argv[5]);
	admitted_run_time=M/2; //6 ms
	q=5e-1;	
	n = N/2;
	numNodes = M/3; // number of nodes for DISTA
	srand(10*run*M);
	VectorXd y,  xtrue, xist, xadmm,  z, dual, noise, auxist, auxadmm, aux, u, xmdista[numNodes], xdista[numNodes], final_mean, zc, xadmmc, dualc, xprevious, auxs, previousglobal, TAUs;
	MatrixXd A, A_T,B, As[numNodes], As_T[numNodes]; 
    A.resize(M,N);
	A_T.resize(N,M);
	B.resize(N,N);
 	xtrue.resize(N);
	xist.resize(N);
	xadmm.resize(N);
	final_mean.resize(N);
	z.resize(N);
	dual.resize(N);
 	u.resize(T);
	y.resize(T);
	noise.resize(T);
	aux.resize(N);
	auxs.resize(M);
	xadmmc.resize(N);
	zc.resize(N);
	dualc.resize(N);
	xprevious.resize(N);
 	TAUs.resize(numNodes);
    measuresNode = M/numNodes;

    
	
	// initialize

    // control input
	for (t = 0; t < M; t++)    {
		u(t)=Gaussian_Noise(1);
        noise(t)= Gaussian_Noise(np);
	}
	for (t = M; t < T; t++)    {
        	//u(t)=(2*(rand() % 2)-1); //a binary sequence
        	//u(t)= (rand() % 2); //a binary sequence
		u(t)=u(t % M);
        noise(t)=Gaussian_Noise(np);
	}

	xtrue.setZero();
	y.setZero();

    xist.setZero();
    xadmm.setZero();
    z.setZero();
    dual.setZero();

    for (j = 0; j < numNodes; j++) {
 		As[j].resize(measuresNode, N);
 		xdista[j].resize (N);
 		xmdista[j].resize (N);
 	}
    
    	
    // parameters to identify
	for (t = 1; t < T; t++)    {
		// a
		if (t < 0.3*T)
			xtrue(0)= -0.9;
		else if (t < 0.5*T)
			xtrue(0)= -0.9;
		else if (t < 0.7*T)
			xtrue(0)= 0.9;
		else
			xtrue(0)= 0.9;
		/*****/

		// b
		if (t < 0.2*T)
			xtrue(n)= 0.7;
		else if (t < 0.4*T)
			xtrue(n)= -0.8;
		else if (t < 0.7*T)
			xtrue(n)= 0.8;
		else
			xtrue(n)= -0.7;

        y(t) = xtrue(0)*y(t-1) + xtrue(n)*u(t-1)+ noise(t);
    }

    
    auto t1 = Clock::now();
	auto t2 = Clock::now();

	// MAIN CYCLE
	for (t = 1; t < T-M; t += M)   {

		// a
		if (t < 0.3*T)
			xtrue(0)= -0.9;
		else if (t < 0.5*T)
			xtrue(0)= -0.9;
		else if (t < 0.7*T)
			xtrue(0)= 0.9;
		else
			xtrue(0)= 0.9;

		/*****/
		// b
		if (t < 0.2*T)
			xtrue(n)=  0.7;
		else if (t < 0.4*T)
			xtrue(n)= -0.8;
		else if (t < 0.7*T)
			xtrue(n)= 0.8;
		else
			xtrue(n)= -0.7;


        // fill A (A is Toeplitz)
		for (i = 0; i < M; i++)  {
		    for (j = 0; j < n; j++)  {
                if (i+t-j-1 <0) {
                    A(i,j)=0;
                    A(i,j+n)=0;
                }
                else {
                	A(i,j)=y(i+t-j-1);
                	A(i,j+n)=u(i+t-j-1);
            	}
        	}
        }
        A=1./sqrt(M)*A;
        A_T=A.transpose();
		TAU=2./(A.squaredNorm());
		thresh =  LAMBDA*TAU;
        thresh_dist =  LAMBDA*TAU/2;
        for (j = 0; j < numNodes; j++) {
 			for (i = 0; i < measuresNode; i++) {
				for (h=0; h< N; h++) {
					As[j](i,h)= sqrt(M)/sqrt(measuresNode)*A(j*measuresNode+i, h) ;
				}
			}
			As_T[j]=As[j].transpose();
			TAUs(j)=2./(As[j].squaredNorm());
 		}
		        
        /** Compute GLOBAL MINIMUM to compate**/
        B=A_T*A+(1+MU)* MatrixXd::Identity(N, N);
		B=B.inverse();
		for (k = 0; k < 1e6; k++)  {   
        	xadmmc=B*(A_T*1./sqrt(M)*y.segment(t,M) + zc - dualc);			
            
            for (j = 0; j < N; j++) {
				zc(j)=soft_thresholding(xadmmc(j)+dualc(j), LAMBDA);
			}
			dualc=dualc+xadmmc-zc;
			aux = xadmmc - xprevious;
			xprevious = xadmmc;
			if (aux.squaredNorm() < 1e-12) {
				//t2 = Clock::now();
				break;
			}
		}
		auxs= y.segment(t,M) - sqrt(M)*A*xadmmc;
		fopt = 0.5*auxs.squaredNorm()+ MU*0.5*xadmmc.squaredNorm()+LAMBDA*xadmmc.lpNorm<1>(); 
		if ( t>1 ) {
			aux = previousglobal - xadmmc;
			gap = gap + aux.lpNorm<2>();
			gap2 = gap2 + aux.squaredNorm();
		}
		previousglobal=xadmmc;
        
        if ( t > 1 ) {
 			auxs= y.segment(t,M) - sqrt(M)*A*xist;
			fist = 0.5*auxs.squaredNorm()+ MU*0.5*xist.squaredNorm()+LAMBDA*xist.lpNorm<1>();
			dinreg_ist = dinreg_ist + fist - fopt;
			dinreg_ist_theo = gap+gap2;

			auxs= y.segment(t,M) - sqrt(M)*A*xadmm;
			fadmm = 0.5*auxs.squaredNorm()+ MU*0.5*xadmm.squaredNorm()+LAMBDA*xadmm.lpNorm<1>(); 
			dinreg_admm = dinreg_admm + fadmm - fopt;

			fdista=0;
			for (v = 1; v < numNodes-1; v++) {
				aux = xmdista[v]-xdista[v];
				fdista = fdista + aux.squaredNorm()/3;
				aux = xmdista[v-1]-xdista[v];
				fdista = fdista + aux.squaredNorm()/3;
				aux = xmdista[v+1]-xdista[v];
				fdista = fdista + aux.squaredNorm()/3;

				auxs= y.segment(t+v*measuresNode,measuresNode) - sqrt(measuresNode)*As[v]*xdista[v];
				fdista =  fdista + 0.5*auxs.squaredNorm()+ (MU*0.5*xdista[v].squaredNorm()+LAMBDA*xdista[v].lpNorm<1>())/numNodes; 
			}
			aux = xmdista[0]-xdista[0];
			fdista = fdista + aux.squaredNorm()/3;
			aux = xmdista[numNodes-1]-xdista[0];
			fdista = fdista + aux.squaredNorm()/3;
			aux = xmdista[1]-xdista[0];
			fdista = fdista + aux.squaredNorm()/3;
			auxs= y.segment(t+0*measuresNode,measuresNode) -sqrt(measuresNode)*As[0]*xdista[0];
			fdista =  fdista + 0.5*auxs.squaredNorm()+ (MU*0.5*xdista[0].squaredNorm()+LAMBDA*xdista[0].lpNorm<1>())/numNodes; 

			aux = xmdista[numNodes-1]-xdista[numNodes-1];
			fdista = fdista + aux.squaredNorm()/3;
			aux = xmdista[numNodes-2]-xdista[numNodes-1];
			fdista = fdista + aux.squaredNorm()/3;
			aux = xmdista[0]-xdista[numNodes-1];
			fdista = fdista + aux.squaredNorm()/3;
			auxs= y.segment(t+(numNodes-1)*measuresNode,measuresNode) -sqrt(measuresNode)*As[numNodes-1]*xdista[numNodes-1];
			fdista =  fdista + 0.5*auxs.squaredNorm()+ (MU*0.5*xdista[numNodes-1].squaredNorm()+LAMBDA*xdista[numNodes-1].lpNorm<1>())/numNodes; 

			dinreg_dista = dinreg_dista + fdista - fopt;

 		}

         /** IST**/
        t1 = Clock::now();
		for (k = 0; k < 1e6; k++)  {
        	xist=xist+TAU*A_T*(1./sqrt(M)*y.segment(t,M)-A*xist)-TAU*MU*xist;      
            		
			for (j = 0; j < N; j++) {
				xist(j)=soft_thresholding(xist(j), thresh);
			}
			
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /1e6 > admitted_run_time) {
				break;
			}
    	}
    	 
 		/** DOUGLAS-RACHFORD (equivalent to ADMM) **/
 		// N.B. The notation "dual" and "admm" comes from previous ADMM-based version; however, this is DR (to be modified) 
 		t1 = Clock::now();
 		B=A_T*A+(1+MU)* MatrixXd::Identity(N, N);
		B=B.inverse();
        
        for (k = 0; k < 1e6; k++)  {   
        	xadmm=B*(A_T*1./sqrt(M)*y.segment(t,M)+dual);			
            
            for (j = 0; j < N; j++) {
				z(j)=soft_thresholding(2*xadmm(j)-dual(j), LAMBDA);
			}
			dual=dual+2*(-xadmm+z);
			
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /1e6 > admitted_run_time) {
				break;
			} 		
		}
		
		/** DISTA**/
        t1 = Clock::now();
		for (k = 0; k < 1e6; k++)  { 
			// FIRST ROUND: COMPUTE THE LOCAL MEANS
			for (v = 1; v < numNodes -1; v++) {
				xmdista[v]=(xdista[v]+ xdista[v+1]+xdista[v-1])/3;
			}
			xmdista[0]=(xdista[0]+xdista[1]+xdista[numNodes-1])/3;
			xmdista[numNodes-1]=(xdista[numNodes-1]+xdista[numNodes-2]+xdista[0])/3;
			
			for (v = 1; v < numNodes-1; v++) {
				// MEAN OF MEANS
				aux=(xmdista[v]+ xmdista[v+1]+xmdista[v-1])/3;
				xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(1./sqrt(measuresNode)*y.segment(t+v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
		
				for (j = 0; j < N; j++) {
					xdista[v](j)=soft_thresholding(xdista[v](j), LAMBDA*TAUs(v)/(2));
				}
			}
			aux=(xmdista[0]+xmdista[1]+xmdista[numNodes-1])/3;
			xdista[0]=(1-q)/(1)*aux+q/1*(xdista[0]+TAUs(0)*As_T[0]*(1./sqrt(measuresNode)*y.segment(t+0*measuresNode,measuresNode)-As[0]*xdista[0])-TAUs[0]*MU*xdista[0]);
			
			for (j = 0; j < N; j++) {
				xdista[0](j)=soft_thresholding(xdista[0](j), LAMBDA*TAUs(0)/(2));
			}


			aux=(xmdista[numNodes-1]+xmdista[numNodes-2]+xmdista[0])/3;
			//xdista[numNodes-1]=(1-q)/(1+MU*TAU)*aux+q/(1+MU*TAUs(numNodes-1))*(xdista[numNodes-1]+TAUs(numNodes-1)*As_T[numNodes-1]*(1./sqrt(measuresNode)*y.segment(t+(numNodes-1)*measuresNode,measuresNode)-As[numNodes-1]*xdista[numNodes-1]));
			xdista[numNodes-1]=(1-q)/(1)*aux+q/(1)*(xdista[numNodes-1]+TAUs(numNodes-1)*As_T[numNodes-1]*(1./sqrt(measuresNode)*y.segment(t+(numNodes-1)*measuresNode,measuresNode)-As[numNodes-1]*xdista[numNodes-1])-TAUs[numNodes-1]*MU*xdista[numNodes-1]);
			
			for (j = 0; j < N; j++) {
				xdista[numNodes-1](j)=soft_thresholding(xdista[numNodes-1](j), LAMBDA*TAUs(numNodes-1)/(2));
			}

			
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /(1e6*numNodes) > admitted_run_time) {
				//cout << k <<" ehi dista"<< endl;
				break;
			}
    	}



		/* output */
		// At each time step t, I print the result
		// I do not save with the script the variables over all time 

		auxist= xist-xtrue;
		auxadmm= xadmm-xtrue;
		auxdista=0;
		final_mean.setZero();
		for (v = 0; v < numNodes; v++) {
			aux=xdista[v]-xtrue;
			auxdista=auxdista+aux.squaredNorm();
			final_mean=final_mean+xdista[v];
		}
		auxdista=auxdista/numNodes;
		final_mean=final_mean/numNodes;
		
		for (i = t; i < t+M; i++)    { 
			if (i+M+admitted_run_time < 0.3*T)
					xtrue(0)= -0.9;
			else if (i+M+admitted_run_time < 0.5*T)
				xtrue(0)= -0.9;
			else if (i+M+admitted_run_time < 0.7*T)
				xtrue(0)= 0.9;
			else
				xtrue(0)= 0.9;

			// b
			if (i+M+admitted_run_time < 0.2*T)
				xtrue(n)=  0.7;
			else if (i+M+admitted_run_time < 0.4*T)
				xtrue(n)= -0.8;
			else if (i+M+admitted_run_time < 0.7*T)
				xtrue(n)= 0.8;
			else
				xtrue(n)= -0.7;
	    		
	    	cout << i+M+admitted_run_time << " "<<xtrue(0) << " " << xist(0) << " " << xtrue(n) << " " << xist(n) << "  ";
	    	// I write the estimation of the null parameters
        	for (j = 1; j < n; j++)    {    
            	cout << xist(j) << "  ";
        	}
        	for (j = n+1; j < N; j++)    {    
        		cout << xist(j) << "  ";
        	}
        	//... and to conclude the MSE
        	cout << (auxist.squaredNorm()) / N <<  "     ";

               
	    	cout << i +M+admitted_run_time << " "<<xtrue(0) << " " << xadmm(0) << " " << xtrue(n) << " " << xadmm(n) << "  ";
	    	// I write the estimation of the null parameters
        	for (j = 1; j < n; j++)    {    
        		cout << xadmm(j) << "  ";
        	}
        	for (j = n+1; j < N; j++)    {    
        		cout << xadmm(j) << "  ";
        	}
        	//... and to conclude the MSE
        	cout << (auxadmm.squaredNorm()) / N <<  "     ";


			cout << i +M+admitted_run_time << " "<<xtrue(0) << " " << final_mean(0) << " " << xtrue(n) << " " << final_mean(n) << "  ";
	    	// I write the estimation of the null parameters
        	for (j = 1; j < n; j++)    {    
            	cout << final_mean(j) << "  ";
        	}
        	for (j = n+1; j < N; j++)    {    
            	cout << final_mean(j) << "  ";
        	}
        	//... and to conclude the MSE
         	cout << auxdista / N <<  "     ";
         	cout << dinreg_ist_theo << "  " << dinreg_ist << "  " << dinreg_admm << "  " << dinreg_dista << "  "<<endl;
    	} // end i
	} //end time t
		
	return EXIT_SUCCESS;
}
