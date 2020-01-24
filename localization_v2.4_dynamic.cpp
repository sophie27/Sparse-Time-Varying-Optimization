/* 7 JAN 2020 */
/** Code written by Sophie M. Fosson **/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <Eigen/Dense>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;
using namespace std;
using namespace Eigen;

#define MAX_ITERATIONS 1e5
#define T 105	
#define LAMBDA 1e-4

#define MU 1e-6
#define TOL 1e-5
//#define ERROR_TOL 1e-10
#define pi 3.1415926535897932384626433832795028841971
#define BETA 1 /** **/

double q, SNR, sigma, SRER, aux, TAU;
double W, np, nonserve, admitted_run_time, auxdistsa, shu;
int H,L,N,M,h,i,j, k,t, t_in, run;
int row, col, new_row, new_col, est_ist_row, est_ist_col, new_est_ist_row, new_est_ist_col, est_admm_row, est_admm_col, new_est_admm_row, new_est_admm_col;
int numNodes, num_neighbors, measuresNode, est_dist_row, est_dist_col, new_est_dist_row, new_est_dist_col,v;

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
//**//


int main (int argc, char *argv[]) {
    if (argc < 5) {
        cerr << "Usage:  ./locrun 1 144 1976 0.0022 > loc.dat" << endl;
        return EXIT_FAILURE;
    }
    k = atoi (argv[1]); // number of targets to detect; set to 1 up to nowm
    M = atoi (argv[2]); // number of measurements (ex: 144)
    run = atoi(argv[3]); // number of the run (any number is ok)
    np = atof(argv[4]); // noise power
   
    H = 25; //height of the grid (# cells)
    L = 25; //length of the grid (# cells)
    W = 100; //width of a square cell (cm)
    N = H*L; //number of  cells = 625

    measuresNode = 4; // measurements for each node
    numNodes = M/measuresNode; // 144/4
    cerr << numNodes << " " << measuresNode << endl;
    q=0.5;
    srand(10*run*M);
    VectorXd y(M),  s(N), xOriginal(N),  xist_old(N), xist(N), xadmm_old(N), xadmm(N), xdist_old(N), xdist(N), noise(M), z(N), z_old(N), dual(N), dual_old(N), auxv, distance_ist(T), distance_admm(T), cumulative_distance_ist(T), cumulative_distance_admm(T);
    MatrixXd A(M,N), AA(M,N), A_T(N,M), AS, A_T_A(N,N), QI, ID,Q(N,N), ORT(N,M), ORT_T(M,N), grid, sensors_pos, d(M,N), A_A_T(M,M), A_PINV(M,N), B(N,N); //X= dynamic d=double
    VectorXi index_nonzero(k), est_ist_index_nonzero(k), est_admm_index_nonzero(k), est_dist_index_nonzero(k);
    VectorXd::Index maxind;
    ID.setIdentity(N,N);
    // FOR DISTA
    VectorXd xmdista[numNodes], xdista[numNodes], xdista_old[numNodes], TAUs(numNodes), aux(N), final_mean(N), distance_dist(T), cumulative_distance_dist(T);
    MatrixXd As[numNodes], As_T[numNodes];

    set<int> inz, einz;
    set<int>::iterator it, oit, uno, due;
    admitted_run_time = 50;  // ms
    /****/

        //set grid
 	grid.resize(2,N);
 	
 	sensors_pos.resize(2,M);
 	for (i = 0; i < N; i++) {
		grid(0,i) = W/2+ (i%L)*W; //cell's coordinates (center), is it?
		grid(1,i) = W/2+ (i/L)*W;
	}
	shu=0; // different from zero if you want to add some randomness to the sensor deployment
	for (i = 0; i < M-measuresNode+1; i = i+measuresNode) {
			sensors_pos(0,i) = L*W/16 + L*W/(6*measuresNode)*(i%(6*measuresNode))+shu*(((double)(rand() % 1000))/1000-0.5);
			sensors_pos(1,i) = H*W/16 + H*W/6*(i/(6*measuresNode))+shu*(((double)(rand() % 1000))/1000-0.5);
		for (h = 1; h < measuresNode; h++) {
			sensors_pos(0,i+h) = sensors_pos(0,i)+shu*(((double)(rand() % 1000))/1000-0.5);
			sensors_pos(1,i+h) = sensors_pos(1,i)+shu*(((double)(rand() % 1000))/1000-0.5);
		}
	}
	
	
   	for (j = 0; j < M; j++) {
		cerr << sensors_pos(0,j) << "  "  << sensors_pos(1,j) << endl;
	}
 	ofstream sensorfile;
  	sensorfile.open ("sensors.dat");
  	for (j = 0; j < M; j++) {
  		sensorfile << sensors_pos(0,j) << "  "  << sensors_pos(1,j) << endl;
  	}
  	sensorfile.close();
  	ofstream connectionfile;
  	connectionfile.open ("connections.dat");
  	for (j=0; j< 6; j++) {
  		connectionfile << sensors_pos(0,j*6*measuresNode) << "  "  << sensors_pos(1,j*6*measuresNode) << endl;
  		connectionfile << sensors_pos(0,(j*6+5)*measuresNode) << "  "  << sensors_pos(1,(j*6+5)*measuresNode) << endl << endl << endl;

  		connectionfile << sensors_pos(0,j*measuresNode) << "  "  << sensors_pos(1,j*measuresNode) << endl;
  		connectionfile << sensors_pos(0,(j+30)*measuresNode) << "  "  << sensors_pos(1,(j+30)*measuresNode) << endl << endl << endl;
  	}
  	
  	connectionfile.close();
	

	// matrix of the distances between cells and sensors (standard: see paper Feng 2009 GLOBECOM)
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			d(i,j)=sqrt( (sensors_pos(0,i) - grid(0,j))*(sensors_pos(0,i) - grid(0,j)) + (sensors_pos(1,i) - grid(1,j))*(sensors_pos(1,i) - grid(1,j)));
			if (d(i,j)< 800)
				A(i,j)=25-40.2-20*log(d(i,j))+Gaussian_Noise(np);
			else
				A(i,j)=25-58.5-33*log(d(i,j))+Gaussian_Noise(np);
		}
	}
	A_T=A.transpose();

    //** ORTHOGONALIZATION see Feng 2009 **//
    JacobiSVD<MatrixXd> svd(A_T, ComputeThinU | ComputeThinV);
    ORT=svd.matrixU();
    ORT_T = ORT.transpose();
    A_A_T=A*A_T;
    A_PINV= A_T*A_A_T.inverse();
	A=ORT_T;
	A_T=A.transpose();

    B=A_T*A+(1+MU)* MatrixXd::Identity(N, N);
	B=B.inverse();
	TAU=2./(A.squaredNorm()); 

 	for (j = 0; j < numNodes; j++) {
 		As[j].resize(measuresNode, N);
 		xdista[j].resize(N);
 		xmdista[j].resize(N);
 		xdista_old[j].resize(N);
 	}


 	for (j = 0; j < numNodes; j++) {
 		for (i = 0; i < measuresNode; i++) {
			for (h=0; h< N; h++) {
				As[j](i,h)= A(j*measuresNode+i, h) ;
			}
		}
		As_T[j]=As[j].transpose();
		TAUs(j)=2./(As[j].squaredNorm());
 	}   
	
	auto t1 = Clock::now();
	auto t2 = Clock::now();

	xOriginal.setZero();
    xOriginal(0) = 1;	
   	row = 0;
   	col = 0; 
   	est_ist_row = 0;
    est_ist_col = 0;

    est_admm_row = 0;
    est_admm_col = 0;
	
	est_dist_row = 0;
    est_dist_col = 0;
    
    z_old.setZero();
    dual_old.setZero();
    xadmm_old.setZero();
    xist_old.setZero();
    for (v = 0; v < numNodes; v++) {
    	xdista_old[v].setZero();
	}
	cumulative_distance_ist(0)=0;
	cumulative_distance_dist(0)=0;
	cumulative_distance_admm(0)=0;

	xOriginal.setZero();
    xOriginal(0) = 1;	
    index_nonzero.setZero();
    index_nonzero(0) = 0;
	
    for (j = 0; j < M; j++) {
		noise(j)=Gaussian_Noise(np);
	}
        
	y = A*xOriginal + noise;
	
	for (t = 1; t < T; t++)    {
		if (t < T/4) {
			new_row = row + ( ((double) rand() / (RAND_MAX)) > 0.4); 
			if (new_row > L-1)
				new_row = L-1;
			else if (new_row < 0)
				new_row = 0;
			new_col = col + (((double) rand() / (RAND_MAX)) > 0.6 );
			if (new_col > L-1)
				new_col = L-1;
			else if (new_col < 0)
				new_col = 0;
			}
			else if (t<2*T/4) {
				new_row = row + (((double) rand() / (RAND_MAX)) > 0.6); 
			if (new_row > L-1)
				new_row = L-1;
			else if (new_row < 0)
				new_row = 0;
			new_col = col + (((double) rand() / (RAND_MAX)) > 0.4);
			if (new_col > L-1)
				new_col = L-1;
			else if (new_col < 0)
				new_col = 0;
			} 
			else if (t<3*T/4) {
				new_row = row - (((double) rand() / (RAND_MAX)) > 0.4); // solo verso destra per un po'
			if (new_row > L-1)
				new_row = L-1;
			else if (new_row < 0)
				new_row = 0;
			new_col = col - (((double) rand() / (RAND_MAX)) > 0.6);
			if (new_col > L-1)
				new_col = L-1;
			else if (new_col < 0)
				new_col = 0;
			} 
			else  {
				new_row = row - (((double) rand() / (RAND_MAX)) > 0.6); // solo verso destra per un po'
			if (new_row > L-1)
				new_row = L-1;
			else if (new_row < 0)
				new_row = 0;
			new_col = col - (((double) rand() / (RAND_MAX)) > 0.4);
			if (new_col > L-1)
				new_col = L-1;
			else if (new_col < 0)
				new_col = 0;
		} 

		cout << row << "  "  << col << "  " << new_row-row << "  " << new_col-col <<  "  "; 
  		
		xist=xist_old;
		t1 = Clock::now();

		//** IST (Iterative Soft Thresholding) **//
		for (t_in=0; t_in < MAX_ITERATIONS; t_in++) {
			xist= xist+TAU*A_T*(y-A*xist)-TAU*MU*xist;


			for (j = 0; j < N; j++) {
				xist(j)=soft_thresholding(xist(j), LAMBDA*TAU);
			}
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /1e6 > admitted_run_time) {
				break;
			}
		}
		xist_old=xist;



		/** DR (Douglas Rachoford) - equivalent to ADMM**/
 		t1 = Clock::now();
 		
        xadmm = xadmm_old;
        dual = dual_old;
        z = z_old;
        for (t_in = 0; t_in < MAX_ITERATIONS; t_in++)  {   
        	xadmm=B*(A_T*y+dual);			
            
            for (j = 0; j < N; j++) {
				z(j)=soft_thresholding(2*xadmm(j)-dual(j), LAMBDA);
			}
			dual=dual+2*(-xadmm+z);
			
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /1e6 > admitted_run_time) {
				break;
			} 		
		}
		xadmm_old = xadmm;
		dual_old = dual;
		z_old = z;

		/** DISTA (Distributed Iterative Soft Thresholding) for "GRID" TOPOLOGY **/
		for (v = 0; v < numNodes; v++) {
    		xdista[v] = xdista_old[v];
    	}
        t1 = Clock::now();

		for (t_in = 0; t_in < MAX_ITERATIONS; t_in++)  { 
			// FIRST ROUND: COMPUTE THE LOCAL MEANS
			
			// central
			for (i = 0; i < 4; i++) {
				for (v = 7+6*i; v < 11+6*i; v++) {	
  					xmdista[v]=(xdista[v]+ xdista[v+1]+xdista[v-1]+xdista[v+6]+xdista[v-6])/5;
				}
			}
			//last row
			for (v = 1; v < 5; v++) {	
  				xmdista[v]=(xdista[v]+ xdista[v+1]+xdista[v-1]+xdista[v+6])/4;
			}
			//first row
			for (v = 31; v < 35; v++) {	
  				xmdista[v]=(xdista[v]+ xdista[v+1]+xdista[v-1]+xdista[v-6])/4;
			}
			//last col
			for (v = 11; v < 30; v=v+6) {	
  				xmdista[v]=(xdista[v]+ xdista[v-1]+xdista[v+6]+xdista[v-6])/4;
			}
			//first col
			for (v = 6; v < 25; v=v+6) {	
  				xmdista[v]=(xdista[v]+ xdista[v+1]+xdista[v-6]+xdista[v+6])/4;
			}

			xmdista[0]=(xdista[0]+xdista[1]+xdista[6])/3;
			xmdista[numNodes-1]=(xdista[numNodes-1]+xdista[numNodes-2]+xdista[numNodes-6])/3;
			xmdista[30]=(xdista[30]+xdista[24]+xdista[31])/3;
			xmdista[5]=(xdista[5]+xdista[4]+xdista[11])/3;
			
			
			// central
			for (i = 0; i < 4; i++) {
				for (v = 7+6*i; v < 11+6*i; v++) {	
  					aux=(xmdista[v]+ xmdista[v+1]+xmdista[v-1]+xmdista[v+6]+xmdista[v-6])/5;
  					xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
				}
			}
			//last row
			for (v = 1; v < 5; v++) {	
  				aux=(xmdista[v]+ xmdista[v+1]+xmdista[v-1]+xmdista[v+6])/4;
  				xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
			}
			//first row
			for (v = 31; v < 35; v++) {	
  				aux=(xmdista[v]+ xmdista[v+1]+xmdista[v-1]+xmdista[v-6])/4;
  				xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
			}
			//last col
			for (v = 11; v < 30; v=v+6) {	
  				aux=(xmdista[v]+ xmdista[v-1]+xmdista[v+6]+xmdista[v-6])/4;
  				xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
			}
			//first col
			for (v = 6; v < 25; v=v+6) {	
  				aux=(xmdista[v]+ xmdista[v+1]+xmdista[v-6]+xmdista[v+6])/4;
  				xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
			}

			v=0;
			aux=(xmdista[0]+xmdista[1]+xmdista[6])/3;
			xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
		
			v=numNodes-1;
			aux=(xmdista[numNodes-1]+xmdista[numNodes-2]+xmdista[numNodes-6])/3;
			xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
		
			v=30;
			aux=(xmdista[30]+xmdista[24]+xmdista[31])/3;
			xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
		
			v=5;
			aux=(xmdista[5]+xmdista[4]+xmdista[11])/3;
			xdista[v]=(1-q)/(1)*aux+q/(1)*(xdista[v]+TAUs(v)*As_T[v]*(y.segment(v*measuresNode,measuresNode)-As[v]*xdista[v])-TAUs[v]*MU*xdista[v]);
		


			for (v=0; v< numNodes; v++) {
				for (j = 0; j < N; j++) {
					xdista[v](j)=soft_thresholding(xdista[v](j), LAMBDA*TAUs(v)/(2));
				}
			}

			aux=(xmdista[0]+xmdista[1]+xmdista[11])/3;
			xdista[0]=(1-q)/(1)*aux+q/1*(xdista[0]+TAUs(0)*As_T[0]*(y.segment(0*measuresNode,measuresNode)-As[0]*xdista[0])-TAUs[0]*MU*xdista[0]);
			
			for (j = 0; j < N; j++) {
				xdista[0](j)=soft_thresholding(xdista[0](j), LAMBDA*TAUs(0)/(2));
			}

			aux=(xmdista[numNodes-1]+xmdista[numNodes-2]+xmdista[24])/3;
			xdista[numNodes-1]=(1-q)/(1)*aux+q/(1)*(xdista[numNodes-1]+TAUs(numNodes-1)*As_T[numNodes-1]*(y.segment((numNodes-1)*measuresNode,measuresNode)-As[numNodes-1]*xdista[numNodes-1])-TAUs[numNodes-1]*MU*xdista[numNodes-1]);
			
			for (j = 0; j < N; j++) {
				xdista[numNodes-1](j)=soft_thresholding(xdista[numNodes-1](j), LAMBDA*TAUs(numNodes-1)/(2));
			}
			
			t2 = Clock::now();
			if (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /(1e6*numNodes) > admitted_run_time) {
				break;
			}
    	}
    	for (v=0; v< numNodes; v++) {
    		xdista_old[v] = xdista[v];
    	}

    	final_mean.setZero();
		for (v = 0; v < numNodes; v++) {
			final_mean=final_mean+xdista[v];
		}
		final_mean=final_mean/numNodes;
		
		/***/
		
		for (i = 0; i < k; i++) {
			nonserve = xist.maxCoeff(&maxind);
			est_ist_index_nonzero(i) = maxind;
		
			nonserve = xadmm.maxCoeff(&maxind);
			est_admm_index_nonzero(i) = maxind;
		
			nonserve = final_mean.maxCoeff(&maxind);
			est_dist_index_nonzero(i) = maxind;
		}
	
		new_est_ist_row= est_ist_index_nonzero(0)/L;
		new_est_ist_col = est_ist_index_nonzero(0)%L;
	
		new_est_admm_row= est_admm_index_nonzero(0)/L;
		new_est_admm_col = est_admm_index_nonzero(0)%L;
	
		new_est_dist_row= est_dist_index_nonzero(0)/L;
		new_est_dist_col = est_dist_index_nonzero(0)%L;
	
		// IST
   		for (i=0; i<k; i++) {
   			inz.insert(new_row*L+ new_col);
   			einz.insert(est_ist_index_nonzero(i));
   		}

		for (it=inz.begin(); it!=inz.end(); ++it) {
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (*it == *oit) {
					inz.erase(it);
					einz.erase(oit);
				}
			}	
		}

   	
   		for (it=inz.begin(); it!=inz.end(); ++it) {
   			distance_ist(t)=1e10;
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*oit) - grid(1,*it), 2.0)<distance_ist(t) ) {
					distance_ist(t)=sqrt(pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*it) - grid(1,*oit), 2.0))/100;
				
					uno=it;
					due=oit;
				}
			}
		
			inz.erase(uno);
			einz.erase(due);
		}
		cumulative_distance_ist(t) = cumulative_distance_ist(t-1) + distance_ist(t);
    
		//DR - ADMM
		for (i=0; i<k; i++) {
   			inz.insert(new_row*L+ new_col);
   			einz.insert(est_admm_index_nonzero(i));
   		}

		for (it=inz.begin(); it!=inz.end(); ++it) {
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (*it == *oit) {
					inz.erase(it);
					einz.erase(oit);
				}
			}	
		}

   	
   		for (it=inz.begin(); it!=inz.end(); ++it) {
   			distance_admm(t)=1e10;
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*oit) - grid(1,*it), 2.0)<distance_admm(t) ) {
					distance_admm(t)=sqrt(pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*it) - grid(1,*oit), 2.0))/100;
					uno=it;
					due=oit;
				}
			}
			inz.erase(uno);
			einz.erase(due);
		}
    	cumulative_distance_admm(t) = cumulative_distance_admm(t-1) + distance_admm(t);
	

		//DISTA
		for (i=0; i<k; i++) {
   			inz.insert(new_row*L+ new_col);
   			einz.insert(est_dist_index_nonzero(i));
   		}

		for (it=inz.begin(); it!=inz.end(); ++it) {
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (*it == *oit) {
					inz.erase(it);
					einz.erase(oit);
				}
			}	
		}

   	
   		for (it=inz.begin(); it!=inz.end(); ++it) {
   			distance_dist(t)=1e10;
			for (oit=einz.begin(); oit!=einz.end(); ++oit) {
				if (pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*oit) - grid(1,*it), 2.0)<distance_dist(t) ) {
					distance_dist(t)=sqrt(pow( grid(0,*oit) - grid(0,*it), 2.0)+ pow( grid(1,*it) - grid(1,*oit), 2.0))/100;
					uno=it;
					due=oit;
				}
			}
			inz.erase(uno);
			einz.erase(due);
		}
		cumulative_distance_dist(t) = cumulative_distance_dist(t-1) + distance_dist(t);
	
		cout << est_ist_row << "  "  << est_ist_col << "  " << new_est_ist_row-est_ist_row << "  " << new_est_ist_col-est_ist_col <<  "  "<< est_admm_row << "  "  << est_admm_col << "  " << new_est_admm_row-est_admm_row << "  " << new_est_admm_col-est_admm_col << "  " << est_dist_row << "  "  << est_dist_col << "  " << new_est_dist_row-est_dist_row << "  " << new_est_dist_col-est_dist_col << " "<<cumulative_distance_ist(t) << " " << cumulative_distance_admm(t) << "  " << cumulative_distance_dist(t) << "   "<<distance_ist(t) << " " << distance_admm(t) << "  " << distance_dist(t) <<endl;
	
		est_ist_row = new_est_ist_row;
		est_ist_col = new_est_ist_col;

		est_admm_row = new_est_admm_row;
		est_admm_col = new_est_admm_col;

		est_dist_row = new_est_dist_row;
		est_dist_col = new_est_dist_col;

		xOriginal.setZero();
		xOriginal(new_row*L+ new_col)= 1;
	
    	for (j = 0; j < M; j++) {
			noise(j)=Gaussian_Noise(np);
		}
        
		y = A*xOriginal + noise;
		row = new_row;
		col = new_col;
	}  /** end t **/

	return EXIT_SUCCESS;
}














