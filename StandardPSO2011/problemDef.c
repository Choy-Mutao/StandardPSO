
//===================================================
struct problem problemDef(int functionCode)
{
	int d;
	struct problem pb;

	int nAtoms; // For Lennard-Jones problem
	static double lennard_jones[14]=
	{-1, -3, -6, -9.103852, -12.7121, -16.505384,-19.821489,-24.113360,-28.422532,
		-32.77,-37.97,-44.33,-47.84,-52.32};


	pb.function=functionCode;	
	// Default values
	// Can be modified below for each function
	pb.epsilon = 0.00000;	// Acceptable error (default). May be modified below
	pb.objective = 0;       // Objective value (default). May be modified below
//	pb.constraintNb=0;
	pb.SS.quantisation=0;		// No quantisation needed (all variables are continuous)
	pb.SS.normalise=0; // Set to a value x.
											//  x>0 => Normalisation is applied (search space => [0,x]^D) 

	// ------------------ Search space
	switch (pb.function)
	{   
		case 0:			// Parabola
			pb.SS.D =30;//  Dimension							

		for (d = 0; d < pb.SS.D; d++)
		{   
			pb.SS.min[d] = -100; // -100
			pb.SS.max[d] = 100;	// 100
			pb.SS.q.q[d] = 0;	// Relative quantisation, in [0,1].   
		}

		pb.evalMax = 75000; //100000;// Max number of evaluations for each run
		pb.epsilon= 0.01; //0.000001; // 1e-3;	
		pb.objective = 0;

		break;
		// Cases 100 etc.
#include "cec2005pb.c"

		case 1:		// Griewank
			pb.SS.D = 30; //30;	

		// Boundaries
		for (d = 0; d < pb.SS.D; d++) 
		{	
			pb.SS.min[d] = -600; 
			pb.SS.max[d] = 600;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 200000;	 
		pb.epsilon=0.01; //0.001; //0.15; //0.001;
		pb.objective=0;
		break;

		case 2:		// Rosenbrock
			pb.SS.D = 30;	// 30

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{	
			pb.SS.min[d] =-30; // -30; 
			pb.SS.max[d] =30; // 30;			
			pb.SS.q.q[d] = 0;	      
		}
		pb.epsilon = 100; //0.05; //.0001;		
		pb.evalMax =75000; //200000; //2.e6;  // 40000 
		pb.objective= 0;


		break;


		case 3:		// Rastrigin
			pb.SS.D =30; // 10;	

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] =-5.12; 
			pb.SS.max[d] =5.12; 	 
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax =75000; //3200; 
		pb.epsilon=50; //0.001;
		pb.objective=0;


		break;

		case 4:		// Tripod
			pb.SS.D = 2;	// Dimension

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 10000; 	
		pb.epsilon=0.0001;
		
	//pb.epsilon=log(1+pb.epsilon);
	
		break;

		case 5: // Ackley
			pb.SS.D = 30; //20;	
		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -32; //-32.768; // 32
			pb.SS.max[d] = 32; //32.768; 
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax =80000; //30000; 
		pb.epsilon=0.000; 
		pb.objective=0;


		break;

		case 6: // Schwefel. Min on (A=420.8687, ..., A)
			pb.SS.D=30;
		//pb.objective=-pb.SS.D*420.8687*sin(sqrt(420.8687));
		pb.objective=-12569.5;
		pb.epsilon=2569.5;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -500;
			pb.SS.max[d] = 500;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 300000;	


		break;

		case 7: // Schwefel 1.2
			pb.SS.D=40;
		pb.objective=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 40000;	


		break;

		case 8: // Schwefel 2.22
			pb.SS.D=30;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 100000;
		pb.objective=0;
		pb.epsilon=0.0001;
		break;

		case 9: // Neumaier 3
			pb.SS.D=40;
		pb.objective=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -pb.SS.D*pb.SS.D;
			pb.SS.max[d] = -pb.SS.min[d];
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 40000;


		break;

		case 10: // G3 (constrained)
			pb.SS.D=10;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 340000;
		pb.objective=0;
		pb.epsilon=1.e-6;
		pb.SS.quantisation=1;
		break;

		case 11: // Network
		//	btsNb=5; bcsNb=2;
		btsNb=19; bcsNb=2;
		pb.SS.D=bcsNb*btsNb+2*bcsNb;
		pb.objective=0;
		for (d = 0; d < bcsNb*btsNb; d++) // Binary representation. 1 means: there is a link
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 1;	
		}
		pb.SS.quantisation=1;
		for (d = bcsNb*btsNb; d < pb.SS.D; d++) // 2D space for the BSC positions
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 20; //15;
			pb.SS.q.q[d] = 0;
		}
		pb.SS.normalise=1;
		pb.evalMax = 5000;
		pb.objective=0;
		pb.epsilon=0;

		break;

		case 12: // Schwefel
			pb.SS.D=3;
		pb.objective=-418.98288727243369*pb.SS.D;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -500;
			pb.SS.max[d] = 500;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 60000; //60000;

		break;

		case 13:		  // 2D Goldstein-Price function (f_min=3, on (0,-1))
			pb.SS.D = 2;	// Dimension
		pb.objective=0;

		pb.SS.min[0] = -100;
		pb.SS.max[0] =100;
		pb.SS.q.q[0] = 0;	
		pb.SS.min[1] = -100;
		pb.SS.max[1] = 100;
		pb.SS.q.q[1] = 0;	
		pb.evalMax = 720;

		break;
		case 14: // Schaffer f6	 
			pb.SS.D = 2;	// Dimension
		pb.objective=0;
		pb.epsilon=0.0001;
		pb.SS.min[0] = -100;
		pb.SS.max[0] =100;
		pb.SS.q.q[0] = 0;	
		pb.SS.min[1] = -100;
		pb.SS.max[1] = 100;
		pb.SS.q.q[1] = 0;	

		pb.evalMax = 30000;


		break;

		case 15: // Step
			pb.SS.D=10;
		pb.objective=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 2500;

		break;

		case 16: // Schwefel 2.21
			pb.SS.D=30;
		pb.objective=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 100000;

		break;

		case 17: // Lennard-Jones
			nAtoms=6; // in {2, ..., 15}
		pb.SS.D=3*nAtoms; 
		pb.objective=lennard_jones[nAtoms-2];
		pb.evalMax =5000+3000*nAtoms*(nAtoms-1) ; // Empirical rule
		pb.epsilon=1.e-6;

pb.SS.normalise=1;
		//pb.SS.D=3*21; pb.objective=-81.684;	
		//pb.SS.D=3*27; pb.objective=-112.87358;
		//pb.SS.D=3*38; pb.objective=-173.928427;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -2;
			pb.SS.max[d] = 2;
			pb.SS.q.q[d] = 0;	
		}

		break;

		case 18: // Gear train
			// solution (16,19,43,49) and equivalent ones (like (19,16,49,43)
			// Success rate 9% is reasonable
			pb.SS.D=4;
		pb.objective=2.7e-12;
		pb.epsilon=1.e-13;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 12;
			pb.SS.max[d] = 60;
			pb.SS.q.q[d] = 1;	
		}
		pb.evalMax = 20000;
		pb.SS.quantisation=1;
		break;

		case 19: // Sine sine function
			pb.SS.D=10;
		pb.objective=-10;// Arbitrary large negative number
		// Remember that the error is abs(f - objective), though
		// Best known (2010-09: -9.5983769). 
		pb.epsilon=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = pi;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 60000;

		break;

		case 20: // Perm function
			pb.SS.D=5;
		pb.objective=0;
		pb.epsilon=0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -pb.SS.D;
			pb.SS.max[d] = pb.SS.D;
			pb.SS.q.q[d] = 1;	
		}
		pb.evalMax = 10000;
		pb.SS.quantisation=1;
		break;

		case 21 : // Compression spring
			pb.SS.D=3;

		pb.SS.min[0] = 1; pb.SS.max[0] = 70; pb.SS.q.q[0] = 1; // N
		pb.SS.min[1] = 0.6; pb.SS.max[1] = 3; pb.SS.q.q[1] = 0; // D
		pb.SS.min[2] = 0.207; pb.SS.max[2] = 0.5; pb.SS.q.q[2] = 0.001; // d

	/*
		pb.SS.min[0] = 2; pb.SS.max[0] = 10; pb.SS.q.q[0] = 1; // N
		pb.SS.min[1] = 0.25; pb.SS.max[1] = 1.3; pb.SS.q.q[1] = 0; // D
		pb.SS.min[2] = 0.05; pb.SS.max[2] = 2; pb.SS.q.q[2] = 0.001; // d
	*/
		
		pb.SS.quantisation=1;
		pb.SS.normalise=1;
		pb.evalMax = 20000; 
		pb.epsilon = 1.e-10;			pb.objective = 2.6254214578; 
		break;
		
		case 22:// Cellular phone
			pb.SS.D=2*5; //2*10  2*nb_of_stations

		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=0;
			pb.SS.max[d]=100; // Warning: hard coded in cellular_phone.c
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 2000*pb.SS.D; 
		pb.epsilon = 1e-9;	
		pb.objective =  0.005530517; // Best known result (2010-01-03)
		// pb.epsilon=0; pb.objective=0;
		break;
		
		case 23:		// Penalized
			pb.SS.D = 30; //30;	

		// Boundaries
		for (d = 0; d < pb.SS.D; d++) 
		{	
			pb.SS.min[d] = -50; 
			pb.SS.max[d] = 50;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 50000; 
		pb.epsilon=0;
		pb.objective=0;

		break;

case 24:// Repulsion (2D)
			pb.SS.D=2*40; // 2*nb_of_charged_points

		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=0;
			pb.SS.max[d]=100; // Warning: hard coded in repulsion.c
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 3000*pb.SS.D; 
		pb.epsilon = 0;	
		pb.objective =  0;
		break;
// case 25:
#include "pressure_vessel.c"

case 26: // Ellipsoidal
			pb.SS.D=30; 

		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-30;
			pb.SS.max[d]=30; // Warning: hard coded in repulsion.c
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 3000; 
		pb.epsilon = 50;	
		pb.objective =  0;

break;

case 999:// for tests
						pb.SS.D=20;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-5;
			pb.SS.max[d]=5;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 10000; 
		pb.epsilon = 0.0001;	
		pb.objective =  -78.3323;

			break; 
		
			pb.SS.D=2;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-512;
			pb.SS.max[d]=512;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 10000; 
		pb.epsilon = 0;	
		pb.objective =  -511.7;

			break; 
			pb.SS.D=2;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-100;
			pb.SS.max[d]=100;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 2000; 
		pb.epsilon = 0;	
		pb.objective =  -1.031628;

break;

		
// Goldstein & Price function
pb.SS.D=2;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-2;
			pb.SS.max[d]=2;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 20000; 
		pb.epsilon = 1.e-6;	
		pb.objective =  3;

break;


break;
// Six-Hump Camel Back
pb.SS.D=2;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-5;
			pb.SS.max[d]=5;
			pb.SS.q.q[d] = 0;	
		}
		pb.evalMax = 20000; 
		pb.epsilon = 1.e-6;	
		pb.objective =  -1.0316;

break;
// Periodic
			pb.SS.D=2;
		for (d = 0; d < pb.SS.D; d++)                  		{
			pb.SS.min[d]=-10;
			pb.SS.max[d]=10;
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 4000; 
		pb.epsilon = 0.001;	
		pb.objective =  0.9;
		break;
}

	if(pb.SS.normalise>0) // Normalise the quanta
	{ 	
			for (d = 0; d < pb.SS.D; d++)
				pb.SS.q.q[d]=pb.SS.q.q[d]*pb.SS.normalise/(pb.SS.max[d]-pb.SS.min[d]);
	}

	pb.SS.q.size = pb.SS.D;

	return pb;
}
#include "lennard_jones.c"