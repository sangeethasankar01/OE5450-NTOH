/********************************************************
******Lid Driven Cavity using Finite Difference Method***
*********************************************************/

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{

// Numerical Parameters
    int N;
    std::cout << "Enter a value for N: ";
    std::cin >> N;
    std::cout << "The value of N is: " << N << std::endl;
    double zeta[N][N], zetat[N][N], psi[N][N], psi0[N][N], psit[N][N], u[N][N], v[N][N];
    double z_resi, p_resi, nu, error1, toler=0.05; //Toler=(d*d/2*nu)
    std::cout << "Enter a value for nu: ";
    std::cin >> nu;
    std::cout << "The value of kinematic viscosity is: " << nu << std::endl;
    const double d = 0.0078;
    const double dt = 0.001;
    double u1,u2,u3,u4,v1,v2,v3,v4; // Velocity components at boundary
    cout<<"Numerical Parameters declared"<<endl;

// Physical Parameters
    double U;
    std::cout << "Enter a value for U: ";
    std::cin >> U;
    std::cout << "Velocity of Lid is: " << U << std::endl;
	const double X = 1.0;
	const double Re = U*X/nu; // Solution run for Re=300, 500, 1000
    cout<<"Physical Parameters assigned and Re.no:"<< Re <<endl;

//Geometry of the obstacle
    double L = 0.2; // length of square obstacle
    double x_obst = 0.5 - L/2; // x-coordinate of centre corner of obstacle
    double y_obst = 0.5 - L/2; // y-coordinate of centre corner of obstacle
    cout<<"Geometry of the obstacle specified"<<endl;

//Boundary conditions (Dirichlet)
    u1=0.0; //bottom boundary condition
    u2=1.0; //top boundary condition
    u3=0.0; //right boundary condition
    u4=0.0; //left boundary condition

    v1=0.0;
    v2=0.0;
    v3=0.0;
    v4=0.0;
//Check for model stability
    double CFL1, CFL2;
    CFL1= (U*dt)/d; // CFL number
    CFL2= (U*dt)/dt;

    if (CFL1 <=1 && CFL2<=1)
    {
        cout<<"Model is Stable = 1"<<endl;
    }
    else
    {
        cout<<"Model is unstable"<<endl;
    }

//1. Initializing Variables
    for (int i = 0;  i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            zeta[i][j] = 0.0;
            zetat[i][j] = 0.0;
            psi[i][j] = 0.0;
            psit[i][j] = 0.0;
            psi0[i][j] = 0.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
    z_resi = 0.0;
    p_resi = 0.0;
    error1 = 0.0;
    cout<<"Initial conditions assigned"<<endl;

// Define Output Streams
    ofstream pStream, zStream;
    pStream.open("output_psi_residual.txt");
    zStream.open("output_vorticity_residual.txt");

// 2. Assigning Boundary Conditions
    // a. Stream Function:
    for (int i = 0; i < N; i++)
    {
        psi[i][0] = 0.0;  // Bottom
        psi[i][N-1] = 0.0;  // Top
        psi[0][i] = 0.0;  // Left
        psi[N-1][i] = 0.0;  // Right
    }
     // Add boundary condition for square obstacle
    for (int i = 1; i < N-1; i++)
    {
        for (int j = 1; j < N-1; j++)
        {
            if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
            {
                psi[i][j] = 0.0; // Inside the obstacle

            }
        }
    }
    // b. Velocities: (No Slip Condition)
    for (int i = 0; i < N; i++)
    {

        u[0][i] = u3;
        u[N-1][i] = u4;
        u[i][0] = u1;
        u[i][N-1] = u2;

        v[0][i] = v3;
        v[N-1][i] = v4;
        v[i][0] = v1;
        v[i][N-1] = v2;
    }
    cout<<"Boundary Conditions assigned"<<endl;
    // Add boundary condition for square obstacle
    for (int i = 1; i < N-1; i++)
    {
        for (int j = 1; j < N-1; j++)
        {
        if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
            {
                u[i][j] = 0.0; // Inside the obstacle
                v[i][j] = 0.0; // Inside the obstacle
            }
        }
    }

// Solution
    do
    {
        //c.Defining Vorticity Boundary Conditions
        for (int i = 0; i < N; i++)
        {
            zeta[i][N-1] = (2*(psi[i][N-1] - U*d -psi[i][N-2]))/(d*d);  // Lid
            zeta[i][0] = -(2*(psi[i][1] - psi[i][0]))/(d*d); // Bottom
            zeta[0][i] = -(2*(psi[1][i] - psi[0][i]))/(d*d); // Left
            zeta[N-1][i] = -(2*(psi[N-2][i] - psi[N-1][i]))/(d*d); // Right
        }

        // Add vorticity boundary condition for square obstacle
        for (int i = 1; i < N-1; i++)
        {
            for (int j = 1; j < N-1; j++)
                {
                if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                    {
                        zeta[i][j] = 0.0; // Inside the obstacle
                    }
                }
        }
        // 3.Solving Vorticity Transport Equation
        for (int i = 1;  i < N-1; i++)
        {
            for (int j = 1; j < N-1; j++)
            {
                zetat[i][j] = zeta[i][j] + dt*((1/Re)*(((zeta[i+1][j] -2*zeta[i][j] + zeta[i-1][j])/(d*d)) +
                ((zeta[i][j+1] -2*zeta[i][j] + zeta[i][j-1])/(d*d)))- ((v[i][j+1]*zeta[i][j+1] - v[i][j-1]*zeta[i][j-1])/(2*d))
                 -((u[i+1][j]*zeta[i+1][j] - u[i-1][j]*zeta[i-1][j])/(2*d)));
                 if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                {
                    zetat[i][j] = 0.0; // Inside the obstacle
                }
            }
        }

        // 4. Solving Stream Poisson
        do{
           error1 = 0;
            // Jacobi Iteration:
            for (int i = 1;  i < N-1; i++)
            {
                for (int j = 1; j < N-1; j++)
                {
                 psit[i][j] = (0.25)*(psi[i+1][j] + psi[i-1][j] + psi[i][j-1] + psi[i][j+1] + zetat[i][j]*(d*d));
                 if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                    {
						psit[i][j] = 0.0; // Inside the obstacle
                    }
                }
            }
            // Error:
            for (int i = 1;  i < N-1; i++)
            {
                for (int j = 1; j < N-1; j++)
                {
                    error1 += pow((psit[i][j] - psi[i][j]),2);
                }
            }
            error1 = sqrt(error1);
            // Assign to Next step:
            for (int i = 1;  i < N-1; i++)
            {
                for (int j = 1; j < N-1; j++)
                {
                    psi[i][j] = psit[i][j];
                    if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                    {
                        psi[i][j] = 0.0; // Inside the obstacle
                    }

                }
            }
        } while (error1> toler);

        //5. Calculation of Velocities
        for (int i = 1; i < N-1; i++)
        {
            for( int j = 1; j < N-1; j++)
            {
                u[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*d);
                v[i][j] = -((psi[i+1][j] - psi[i-1][j])/(2*d));
                if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                    {
                        u[i][j] = 0.0;
                        v[i][j] = 0.0; // Inside the obstacle
                    }
            }
        }
        //6. Convergence Criteria
        for (int i = 1; i < N-1; i++)
        {
            for( int j = 1; j < N-1; j++)
            {
                z_resi += abs(zetat[i][j] - zeta[i][j]);
                p_resi += abs(psit[i][j] - psi0[i][j]);
            }
        }
        z_resi = (z_resi)/(N*N);
        p_resi = (p_resi)/(N*N);

        pStream << p_resi << endl;
        zStream << z_resi << endl;
        // Assign t+1 to t
        for (int i = 0;  i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                zeta[i][j] = zetat[i][j];
                psi[i][j] = psit[i][j];
                psi0[i][j] = psit[i][j];
                if (i*d > x_obst && i*d < x_obst + L && j*d > y_obst && j*d < y_obst + L)
                    {
                        zeta[i][j] = 0.0;
                        psi0[i][j] = 0.0;
                        psi[i][j] = 0.0; // Inside the obstacle
                    }
            }
        }
    } while (z_resi > 1e-9 || p_resi > 1e-9);
    cout<<"Vorticity Stream function solved"<< endl;
    pStream.close();
    zStream.close();

// Printing results
    ofstream oStream;
    //b.u velocity output
    oStream.open("output_u.txt");

    for (int j = N-1;  j >= 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            oStream << u[i][j] << " ";
        }
        oStream << endl;
    }
    oStream.close();

    //c.v velocity output
    oStream.open("output_v.txt");

    for (int j = N-1;  j >= 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            oStream << v[i][j] << " ";
        }
        oStream << endl;
    }
    oStream.close();

    //d.Total velocity output
    oStream.open("output_Total_vel.txt");

    for (int j = N-1;  j >= 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            oStream << sqrt(pow(u[i][j],2) + pow(v[i][j],2)) << " ";
        }
        oStream << endl;
    }
    oStream.close();

    //e.Stream Function output:
    oStream.open("output_psi.txt");

    for (int j = N-1;  j >= 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            oStream << psi[i][j] << " ";
        }
        oStream << endl;
    }
    oStream.close();

    //f.Vorticity output:
    oStream.open("output_zeta.txt");

    for (int j = N-1;  j >= 0; j--)
    {
        for (int i = 0; i < N; i++)
        {
            oStream << zeta[i][j] << " ";
        }
        oStream << endl;
    }
    oStream.close();

	//g.Writing center line data to file, Yc vs Uc
	oStream.open("YcUc.txt");
	for(int j=0; j<N; j++)
	{
		oStream << j*d << "\t" << (u[N/2][j]+u[-1+N/2][j])/2 << endl;
	}
	oStream.close();

	//h.Writing Xc vs Vc data to file
	oStream.open("XcVc.txt");
	for(int i=0; i<N; i++)
	{
		oStream<< i*d << "\t" << (v[i][N/2] + v[i][-1+ N/2])/2 << endl;
	}
	oStream.close();
    return 0;
}
