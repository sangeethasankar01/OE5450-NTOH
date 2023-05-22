/********************************************************************
******Lid Driven Cavity with Obstacle using lattice Boltzmann method**
**********************************************************************/
//// D2Q9 LBM model.....

#include <iostream>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <cmath>
#define pi 3.14159

using namespace std;

int main()
{
	clock_t ct;
	ct=clock(); //Variable to measure CPU time
// Numerical Parameters
	int t_max=10000,t; //#Total number of time-steps : change if needed(to ensure convergence)
	int N=100,M=100;	//No. of grid points
	double Lx=1.0;	//Length in x-axis
	double Ly=1.0;	//Lenght in y-axis
	double U= 1.0;  //Velocity of top lid
	double t1,t2,temp;
    double error=1.0; //// Initial error (greater than tolerance)
    double toler=1e-6; //Tolerance

// Physical Parameters
	double nu=0.0025; // Kinematic Viscosity
	const double Re = U*Lx/nu; // Solution run for Re=300, 500, 1000
	const double dx = 0.01, dy= 0.01;
    double dt = 0.01;
	double rho_0 = 1.0; // refrence density
	int c=(1/sqrt(3))*dx/dt, Q=9; // Speed of lattice sound (equal to 1) and Q=direction index of D2Q9
    cout<<"Physical Parameters assigned and Re.no:"<< Re <<endl;

// Parameters for LBM distributions
	double omega=2.0/(1.0 + 6.0*(nu*dt/(dx*dx))),tauf= 1.0/omega ;  //Relaxation time;
	cout << "The value of tau is:"<< tauf << endl;
	double rho[N+1][M+1], ux[N+1][M+1], vy[N+1][M+1]; //2d array for storing macroscopic density and velocity
    double f[Q][N+1][M+1], feq[Q][N+1][M+1];  // distribution function and equilibrium distribution function

// Define the vorticity and stream function
    double zeta[N+1][M+1]; // vorticity
    double psi[N+1][M+1],psi_new[N+1][M+1];   // stream function

//Geometry of the obstacle
    double L = 0.2; // length of square obstacle
    double x_obst = 0.5 - L/2; // x-coordinate of centre corner of obstacle
    double y_obst = 0.5 - L/2; // y-coordinate of centre corner of obstacle
    cout<<"Geometry of the obstacle specified"<<endl;

//Lattice Weights in equilibrium distribution function
	double w[9]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

//assigning lattice particle velocities for D2Q9
	int ex[9]={0,1,0,-1,0,1,-1,-1,1}, ey[9]={0,0,1,0,-1,1,1,-1,-1};
	cout<<"weights and particle velocity for model assigned"<<endl;

// Boundary conditions (Dirichlet)
    const double u1 = 0.0; // bottom boundary condition
    const double u2 = 0.1; // top boundary condition
    const double u3 = 0.0; // right boundary condition
    const double u4 = 0.0; // left boundary condition

    const double v1 = 0.0;
    const double v2 = 0.0;
    const double v3 = 0.0;
    const double v4 = 0.0;

//Check for model stability
    double CFL1, CFL2;
    CFL1= (U*dt)/dx; // CFL number
    CFL2= (U*dt)/dt;
    if (CFL1 <=1 && CFL2<=1)
    {
        cout<<"Model is Stable = 1"<<endl;
    }
    else
    {
        cout<<"Model is unstable"<<endl;
    }

// initialize density and velocity distribution functions
	for (int j = 0; j < M+1; j++)
        {
            for (int i = 0; i< N+1; i++)
            {
					rho[i][j] = rho_0;
					ux[i][j] = 0.0;
					vy[i][j] = 0.0;
					psi[i][j] =0.0;
					psi_new[i][j] =0.0;
            }
        }
  
	//velocity and stream function at wall boundary
		for (int i=0; i<N+1; i++)
	    {
             //Top lid boundary condition
            ux[i][M] = u2;
            vy[i][M] = v2;
            psi[i][M] = 0.0;
            // Bottom wall boundary condition
            ux[i][0] = u1;
            vy[i][0] = v1;
            psi[i][0] = 0.0;
            // Left walls boundary condition
            ux[0][i] = u3;
            vy[0][i] = v3;
            psi[0][i] = 0.0;
            // Right walls boundary condition
            ux[M][i] = u4;
            vy[M][i] = v4;
            psi[M][i] = 0.0;
	    }
		for (int j = 0; j < M+1; j++) 
		{
			for (int i = 0; i < N+1; i++) 
			{
				if (i*dx >= x_obst && i*dx <= x_obst + L && j*dy >= y_obst && j*dy <= y_obst + L) 
				{
					ux[i][j] = 0.0; // Set velocity to zero at obstacle boundary
					vy[i][j] = 0.0; // Set velocity to zero at obstacle boundary
					zeta[i][j] = 0.0; // Set stream function to zero at obstacle boundary
				}
			}
		}
//Initializing the f with feq at t=0
	for(int j=0; j<M+1; j++)
	{
			for(int i=0; i<N+1; i++)
			{
					t1=ux[i][j]*ux[i][j]+vy[i][j]*vy[i][j];
					for(int k=0; k<Q; k++)
					{
						t2=ux[i][j]*ex[k]+vy[i][j]*ey[k];
						f[k][i][j]=rho[i][j]*w[k]*(1.0+3.0*t2+4.5*t2*t2-1.50*t1);
					}
			}
	}
//Main Loop Starts to compute field variables
	for(t=0; t<t_max && error > toler; t++)
	{
       	//COLLISION STEP (evolution Lattice Boltzmann Equation)
		for(int j=0; j<M+1; j++)
        {
			for(int i=0; i<N+1; i++)
				{
					t1=ux[i][j]*ux[i][j]+vy[i][j]*vy[i][j];
					for(int k=0; k<Q; k++)
					{
						t2=ux[i][j]*ex[k]+vy[i][j]*ey[k];
						feq[k][i][j]=rho[i][j]*w[k]*(1.0+3.0*t2+4.5*t2*t2-1.50*t1);
                        f[k][i][j]=omega*feq[k][i][j]+(1.0-omega)*f[k][i][j];
					}
				}
		}
        
		//STREAMING STEP
		for(int j=0; j<M+1; j++)
		{
			for(int i=N; i>=1; i--)
            {
				f[1][i][j]=f[1][i-1][j]; //Right to Left
            }
			for(int i=0; i<=N-1; i++) //Left to Right
            {
				f[3][i][j]=f[3][i+1][j];
            }
		}
		for(int j=M; j>=1; j--) //Top to Bottom
		{
			for(int i=0; i<=N; i++)
            {
				f[2][i][j]=f[2][i][j-1];
            }

			for(int i=N; i>=1; i--)
            {
				f[5][i][j]=f[5][i-1][j-1];
            }

			for(int i=0; i<=N-1; i++)
            {
				f[6][i][j]=f[6][i+1][j-1];
            }
		}
		for(int j=0; j<=M-1; j++) //Bottom to Top
		{
			for(int i=0; i<=N; i++)
            {
				f[4][i][j]=f[4][i][j+1];
            }

			for(int i=0; i<=N-1; i++)
            {
				f[7][i][j]=f[7][i+1][j+1];
            }

			for(int i=N; i>=1; i--)
            {
				f[8][i][j]=f[8][i-1][j+1];
            }
		}
		//BOUNDARY CONDITIONS
		double rh;
		for(int j=0; j<M+1; j++)
		{
			//Bounce Back on West Boundary
			f[1][0][j]=f[3][0][j];
			f[5][0][j]=f[7][0][j];
			f[8][0][j]=f[6][0][j];

			//Bounce Back on East Boundary
			f[3][N][j]=f[1][N][j];
			f[7][N][j]=f[5][N][j];
			f[6][N][j]=f[8][N][j];
		}
		//Bounce Back on South Boundary
		for(int i=0; i<N+1; i++)
		{
			f[2][i][0]=f[4][i][0];
			f[5][i][0]=f[7][i][0];
			f[6][i][0]=f[8][i][0];
		}
		//Moving Lid, North boundary
		for(int i=1; i<N; i++)
		{
			rh=f[0][i][M]+f[1][i][M]+f[3][i][M]+2.0*(f[2][i][M]+f[6][i][M]+f[5][i][M]);
			f[4][i][M]=f[2][i][M];
			f[7][i][M]=f[5][i][M] + 0.5*(f[1][i][M]-f[3][i][M]) - 0.5*rh*0.1;
            f[8][i][M]=f[6][i][M] + 0.5*(f[3][i][M]-f[1][i][M]) + 0.5*rh*0.1;
		}
		//Bounce Back on obstacle boundary
		for (int j = 0; j < M+1; j++) {
			for (int i = 0; i < N+1; i++) {
				if (i*dx >= x_obst && i*dx <= x_obst + L && j*dy >= y_obst && j*dy <= y_obst + L) 
				{
					// Inside obstacle, set bounce back condition
					f[1][i][j] = f[3][i][j];
					f[3][i][j] = f[1][i][j];
					f[2][i][j] = f[4][i][j];
					f[4][i][j] = f[2][i][j];
					f[5][i][j] = f[7][i][j];
					f[7][i][j] = f[5][i][j];
					f[6][i][j] = f[8][i][j];
					f[8][i][j] = f[6][i][j];
				}
			}
		}
	//COMPUTATION OF RHO, U, V
		double ssum, usum, vsum;
		//Computation of the density
		for(int j=0; j<=M; j++)
		{
			for(int i=0; i<=N; i++)
			{
				ssum=0.0;
				for(int k=0; k<=8; k++)
				{
                        ssum+=f[k][i][j];
                }
				rho[i][j]=ssum;
            }
		}
		//Computation of u, v velocities
		for(int i=0; i<N+1; i++)
		{
			for(int j=0; j<M+1; j++)
			{
				usum=0.0;
				vsum=0.0;
				for(int k=0; k<Q; k++)
				{            
                       usum+=f[k][i][j]*ex[k];
                       vsum+=f[k][i][j]*ey[k];
				}
				ux[i][j]=usum/rho[i][j];
				vy[i][j]=vsum/rho[i][j];
			}
		}
        // Calculate Stream function
            for(int j=0; j<M+1; j++)
                {
                    for(int i=0; i<N+1; i++)
                    {
                       zeta[i][j] = (vy[i+1][j]-vy[i-1][j])/(2*dx)-(ux[i][j+1]-ux[i][j-1])/(2*dy);
                    }
                }
		for (int i=1; i<M; i++)
	    {
             //Top lid boundary condition
            ux[i][M] = u2;
            vy[i][M] = v2;
            psi[i][M] = 0.0;
            // Bottom wall boundary condition
            ux[i][0] = u1;
            vy[i][0] = v1;
            psi[i][0] = 0.0;
            // Left walls boundary condition
            ux[0][i] = u3;
            vy[0][i] = v3;
            psi[0][i] = 0.0;
            // Right walls boundary condition
            ux[M][i] = u4;
            vy[M][i] = v4;
            psi[M][i] = 0.0;
        }
		//Velocity and stream function at obstacle boundary
		for (int j = 1; j < M; j++) 
		{
			for (int i = 1; i < N; i++) 
			{
				if (i*dx >= x_obst && i*dx <= x_obst + L && j*dy >= y_obst && j*dy <= y_obst + L) 
				{
					ux[i][j] = 0.0; // Set velocity to zero at obstacle boundary
					vy[i][j] = 0.0; // Set velocity to zero at obstacle boundary
					zeta[i][j] = 0.0; // Set stream function to zero at obstacle boundary
				}
			}
		}

		if(t%1000==0)
        cout << "Timestep = " << t << endl;
	}//End of Main Loop
	cout<<"Lattice Boltmaan method solved"<<endl;

// Printing Results
		ofstream oStream;
		//u velocity output
    	oStream.open("output_u.txt");
    	for (int j = M;  j >= 0; j--)
    	{
        	for (int i = 0; i < N+1; i++)
        	{
            	oStream << ux[i][j] << " ";
        	}
        	oStream << endl;
   		 }
    	oStream.close();

		//v velocity output
		oStream.open("output_v.txt");
    	for (int j = M;  j >= 0; j--)
    	{
        	for (int i = 0; i < N+1; i++)
        	{
            	oStream << vy[i][j] << " ";
        	}
        	oStream << endl;
   		 }
    	oStream.close();

        //d.Total velocity output
        oStream.open("output_Total_vel.txt");
        for (int j = M;  j >= 0; j--)
        {
            for (int i = 0; i < N+1; i++)
            {
                oStream << sqrt(pow(ux[i][j],2) + pow(vy[i][j],2)) << " ";
            }
            oStream << endl;
        }
        oStream.close();

        // Vorticity output
        oStream.open("output_zeta.txt");
    	for (int j = M;  j >= 0; j--)
    	{
        	for (int i = 0; i < N+1; i++)
        	{
            	oStream << zeta[i][j] << " ";
        	}
        	oStream << endl;
   		 }
   		   oStream.close();

   		 //g.Writing center line data to file, Yc vs Uc
            oStream.open("YcUc.txt");
            for(int j=0; j<M+1; j++)
            {
                oStream << j*dy << "\t" << (ux[M/2][j]+ux[-1+M/2][j])/2 << endl;
            }
            oStream.close();

            //h.Writing Xc vs Vc data to file
            oStream.open("XcVc.txt");
            for(int i=0; i<N+1; i++)
            {
                oStream<< i*dx << "\t" << (vy[i][N/2] + vy[i][-1+ N/2])/2 << endl;
            }
            oStream.close();

	//Displaying CPU time
	ct=clock()-ct;
	cout << "\n\nTime Taken = " << (double)(ct/CLOCKS_PER_SEC)/60.00 << " mins." << endl;
	return 0;
}
