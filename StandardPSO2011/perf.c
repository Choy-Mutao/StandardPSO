double perf (struct position x, int function, struct SS SS, double objective) 
{				// Evaluate the fitness value for the particle of rank s 
	double beta;  
	double c;
	int d;
	double DD;
	double dx1,dx2;
	int grid;
	int i,j;
	int  k;
	double min,max;
	int n;
		struct fitness ff={0};
	double f, p, xd, x1, x2,x3,x4,x5,x6;
	double s11, s12, s21, s22;
	double sum1,sum2;
	double t0, tt, t1;
	double theta;
	double u;
	struct position xs; 
	double y1,y2;
	double z;
	#include "cec2005data.c"
	
	//--------- Test
	struct problem pb;
	struct param param;
	
	// See function 11 (Network)
/*
	 static float bts[5][2]=
	 {
		 {6.8, 9.0},
		 {8.3, 7.9},
		 {6.6, 5.6},
		 {10, 5.4},
		 {8, 3} 
	 };
*/
	static float bts [19][2]=
	{
		{6, 9},
		{8, 7},
		{6, 5},
		{10, 5},
		{8, 3} ,
		{12, 2},
		{4, 7},
		{7, 3},
		{1, 6},
		{8, 2},
		{13, 12},
		{15, 7},
		{15, 11},
		{16, 6},
		{16, 8},
		{18, 9},
		{3, 7},
		{18, 2},
		{20, 17}
	};
	
	float btsPenalty= 100;

	double z1,z2;

	if(SS.normalise>0) 
	{
		// Back to the real search space
			xs.size=x.size;
		for(d=0;d<xs.size;d++)
			xs.x[d]=SS.min[d]+(SS.max[d]-SS.min[d])*x.x[d]/SS.normalise;
	}
	else xs=x;

	switch (function)
	{
	#include "cec2005.c"
	
				case -1:		// Constant. For test of biases
				f = 0; 
			break;

			case 0:		// Parabola (Sphere)
				f = 0;

			for (d = 0; d < xs.size; d++) 
		{    
			xd = xs.x[d]-d;   
			f = f + xd * xd;    
		}	  
			break;

			case 1:		// Griewank
				f = 0; 
			p = 1;

			for (d = 0; d < xs.size; d++)
		{      
			xd = xs.x[d];
			f = f + xd * xd;	      
			p = p * cos (xd / sqrt ((double) (d + 1)));	    
		} 
			f = f / 4000 - p + 1;	  
			break;

			case 2:		// Rosenbrock
				f = 0;  
			t0 = xs.x[0] + 1;	// Solution on (0,...0) when
			// offset=0
			for (d = 1; d < xs.size; d++)
		{     

			t1 = xs.x[d]  + 1;	      
			tt = 1 - t0;	      
			f += tt * tt;      
			tt = t1 - t0 * t0;      
			f += 100 * tt * tt;	      
			t0 = t1;    
		}  
			break;

			case 3:		// Rastrigin
				k = 10;  
			f = 0;

			for (d = 0; d < xs.size; d++)    
		{     
			xd = xs.x[d];
			f =f+ xd * xd - k * cos (2 * pi * xd);	    
		}	  
			f =f+ xs.size * k;  
			break;

			case 4:		// 2D Tripod function
				// Note that there is a big discontinuity right on the solution
				// point. 
				x1 = xs.x[0] ; 
			x2 = xs.x[1];  
			s11 = (1.0 - sign (x1)) / 2;
			s12 = (1.0 + sign (x1)) / 2; 
			s21 = (1.0 - sign (x2)) / 2;
			s22 = (1.0 + sign (x2)) / 2;

			//f = s21 * (fabs (x1) - x2); // Solution on (0,0)
			f = s21 * (fabs (x1) +fabs(x2+50)); // Solution on (0,-50)  
			f = f + s22 * (s11 * (1 + fabs (x1 + 50) +
			                      fabs (x2 - 50)) + s12 * (2 +
			                                               fabs (x1 - 50) +
			                                               fabs (x2 - 50)));	  
			
	//f=log(1+f);	
			break;

			case 5:  // Ackley
				sum1=0;
			sum2=0;
			DD=x.size;
			pi=acos(-1);
			for (d=0;d<x.size;d++)
		{
			xd=xs.x[d];
			sum1=sum1+xd*xd;
			sum2=sum2+cos(2*pi*xd);
		}
			f=-20*exp(-0.2*sqrt(  sum1/DD  ))-exp(sum2/DD)+20+exp(1);

			break;

			case 6: // Schwefel
				f=0;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d];
			f=f-xd*sin(sqrt(fabs(xd)));
		}
		break;

			case 7: // Schwefel 1.2
				f=0;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d];
			sum1=0;
			for(k=0;k<=d;k++) sum1=sum1+xd;
			f=f+sum1*sum1;
		}
			break;

			case 8: // Schwefel 2.22
				sum1=0; sum2=1;
			for (d=0;d<x.size;d++)
		{
			xd = fabs(xs.x[d]);
			sum1=sum1+xd;
			sum2=sum2*xd;
		}
			f=sum1+sum2;
			break;

			case 9: // Neumaier 3
				sum1=0; sum2=1;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d]-1;
			sum1=sum1+xd*xd;
		}
			for (d=1;d<x.size;d++)
		{
			sum2=sum2+ xs.x[d]* xs.x[d-1];
		}	

			f=sum1+sum2;
			break;

			case 10: // G3 (constrained) 
							// min =0 on (1/sqrt(D), ...)
				f=1;
			sum1=0;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d];
			f=f*xd;
			sum1=sum1+xd*xd;
		}
			f=fabs(1-pow(x.size,x.size/2)*f) + x.size*fabs(sum1-1);
			break;

			case 11: // Network  btsNb BTS, bcdNb BSC

				f=0;
			// Constraint: each BTS has one link to one BSC 
			for(d=0;d<btsNb;d++)
		{
			sum1=0;
			for(k=0;k<bcsNb;k++) sum1=sum1+xs.x[d+k*btsNb];
			if(sum1<1-zero || sum1>1+zero) f=f+btsPenalty;	

		}
			// Distances
			for(d=0;d<bcsNb;d++) //For each BCS d
		{	
			for(k=0;k<btsNb;k++) // For each BTS k
			{
				if(xs.x[k+d*btsNb]<1) continue;
				// There is a link between BTS k and BCS d
				n=bcsNb*btsNb+2*d;
				z1=bts[k][0]-xs.x[n];
				z2=bts[k][1]-xs.x[n+1];		
				f=f+sqrt(z1*z1+z2*z2);
			}
		}
			break;

		case 12: // Schwefel
			f=0;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d];
			f=f-xd*sin(sqrt(fabs(xd)));
		}	
			break;

			case 13: // 2D Goldstein-Price function
				x1=xs.x[0]; x2=xs.x[1];

			f= (1 + pow(x1 + x2 + 1, 2) *(19-14 *x1 + 3*x1*x1-14* x2 + 6* x1* x2 + 3*x2*x2 ))
				* (30 + pow(2* x1 - 3*x2 ,2)*
				   (18 -32 *x1 + 12 *x1*x1 + 48* x2 - 36* x1 *x2 + 27* x2*x2 ));
			break;

			case 14:  //Schaffer F6
				x1=xs.x[0]; x2=xs.x[1];
			f= 0.5 + (pow(sin(sqrt(x1*x1 + x2*x2)),2) - 0.5)/pow(1.0 + 0.001*(x1*x1 + x2*x2),2); 

			break;

			case 15: // Step
				f=0;
			for (d=0;d<x.size;d++)
		{
			xd = (int)(xs.x[d]+0.5);
			f=f+xd*xd;
		}	
			break;

 	case 16: // Schwefel 2.21
				f=0;
			for (d=0;d<x.size;d++)
		{
			xd = fabs(xs.x[d]);
			if(xd>f) f=xd;
		}
			break;
			
			case 17: // Lennard-Jones
			f=lennard_jones(xs);
			break;
			
			case 18: // Gear train
			f=pow(1./6.931 -x.x[0]*x.x[1]/(x.x[2]*x.x[3]),2);
		//	f=pow(fabs(1./6.0 -x.x[0]*x.x[1]/(x.x[2]*x.x[3])),2);

			break;
			
			case 19: // Sine-sine function
			f=0;
			for (d=0;d<x.size;d++)
		{
			xd = xs.x[d];
			f=f-sin(xd)*pow(sin((d+1)*xd*xd/pi),20);
	
		}
		break;
		
		case 20: // Perm function
		beta=10;
		f=0;
		for (k=0;k<x.size;k++)
		{
			sum1=0; 
			for (d=0;d<x.size;d++)
			{
				xd = xs.x[d];
				sum1=sum1+  ( pow(d+1,k)+beta)*(pow(xd/(d+1),k)-1);
			}
			sum1=sum1*sum1;
		f=f+sum1;	
		}
			
		break;
		
case 21: // Coil compression spring  (penalty method)
			// Ref New Optim. Tech. in Eng. p 644

		x1=xs.x[0]; // {1,2, ... 70}
		x2= xs.x[1];//[0.6, 3]
		x3= xs.x[2];// relaxed form [0.207,0.5]  dx=0.001
		// In the original problem, it is a list of
		// acceptable values
		// {0.207,0.225,0.244,0.263,0.283,0.307,0.331,0.362,0.394,0.4375,0.5}

		f=pi*pi*x2*x3*x3*(x1+2)*0.25;
		//	f=x2*x3*x3*(x1+2);
		// Constraints
		ff=constraint(xs,function);

			if (ff.f[1]>0) {c=1+ff.f[1]; f=f*c*c*c;}
			if (ff.f[2]>0) {c=1+ff.f[2]; f=f*c*c*c;}
			if (ff.f[3]>0) {c=1+ff.f[3]; f=f*c*c*c;}
			if (ff.f[4]>0) {c=1+pow(10,10)*ff.f[4]; f=f*c*c*c;}
			if (ff.f[5]>0) {c=1+pow(10,10)*ff.f[5]; f=f*c*c*c;}
		break;
		
 case 22: //Cellular phone
#include "cellular_phone.c"
	break;
		
		case 23: // Penalized
		f=pow(sin(pi*xs.x[0]),2);
		for (d=1;d<x.size-1;d++)
		{
			f=f+pow(xs.x[d],2)*(1+pow(sin(3*pi*xs.x[d+1]) ,2));
		}
		f=0.1*(f+pow(xs.x[x.size-2],2)*(1+pow(sin(2*pi*xs.x[x.size-1]),2)));
		
		for (d=0;d<x.size;d++)
		{
			xd=xs.x[d];
			if(xd>5) {u=100*pow(xd-5,4); f=f+u;}
			if(xd<-5) {u=100*pow(-xd-5,4); f=f+u;}
		}
		
		break;

case 24: // Repulsion
#include "repulsion.c"
break;

	case 25: // Pressure vessel (penalty method)
		case 1025: // confinement method
#include "perf_pressure_vessel.c"
	break;
	
case 26: // Ellipsoidal
				f = 0;

			for (d = 0; d < xs.size; d++) 
		{    
			xd = xs.x[d]-d-1;   
			f = f + xd * xd;    
		}	 
break;

			case 27:		// Quadric
				f = 0;
			for (d = 0; d < xs.size; d++) 
		{    
			xd=xs.x[0];
			if(d>0)
			  for(j=1;j<d;j++) xd=xd+xs.x[j];
	   
			f = f + xd * xd;    
		}	  
			break;
			
			case 28:// Frequency modulation sound parameter identification
			theta=2*pi/100;
			f=0;
			x1=xs.x[0];x2=xs.x[1];x3=xs.x[2];x4=xs.x[3];x5=xs.x[4];x6=xs.x[5];
			
			for(d=1;d<=100;d++)
			{
			  z=x1*sin(x2*d*theta+x3*sin(x4*d*theta+x5*sin(x6*d*theta)))
			  -sin(5*d*theta+1.5*sin(4.8*d*theta+2*sin(4.9*d*theta)));
			  f=f+z*z;
			}
			break;	
				
		case 999: // for tests
f=0;
			for (d = 0; d < xs.size; d++) 
		{
				xd = xs.x[d];
			f=f+pow(xd,4)-16*xd*xd+5*xd;
		}
f=f/xs.size;
			break;

// Rana
	x1=xs.x[0]; x2=xs.x[1];
		f=x1*sin(sqrt(fabs(x2+1-x1)))*cos(sqrt(fabs(x2+1+x1)))
			+(x2+1)*cos(sqrt(fabs(x2+1-x1)))*sin(sqrt(fabs(x2+1+x1)));
			break;

	// Goldstein & Price function	
		x1=xs.x[0]; x2=xs.x[1];
		f=(1+pow(x1+x2+1,2)*(19-14*x1+13*x1*x1-14*x2	+6*x1*x2+3*x2*x2));
		break;
		
	// Six-Hump Camel Back
	x1=xs.x[0]; x2=xs.x[1]; 
	f=4*x1*x1-2.1*pow(x1,4)+pow(x1,6)/3+x1*x2-4*x2*x2+4*pow(x2,4);

break;
// Periodic
x1=xs.x[0]; x2=xs.x[1];
f=1+pow(sin(x1),2) + pow(sin(x2),2) -0.1*exp(-x1*x1-x2*x2);

break;

case 1004: // For Tripod
pb.function=4;
			pb.SS.D = 2;	// Dimension

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;	
		}

		pb.evalMax = 100; 	 // 10000
		pb.epsilon=0.0001;
		pb.objective=0;
//--
param.BW[0]=0; 	
param.BW[1]=0; 	
param.BW[2]=4;
param.BW[3]=0;	
param.confin=0; 
param.distrib=0; 
param.mean=0.5; 
param.sigma=1./12; // Default: 1./12 (standard deviation of U(0,1))			

param.S=40; 
param.K=3; 	
param.p=1-pow(1-1./param.S,param.K); 

param.w = 1. / (2 * log ((double) 2)); // 0.721
param.c = 0.5 + log ((double) 2); // 1.193
param.topology = 0; 
param.trace=0;

nCycle=0;
nCycleMax=3;

for(d=0;d<3;d++) 
{randNumber[d]=xs.x[d];
printf("%f ", randNumber[d]);
}
f=1;
for(i=0;i<2;i++) // 100
{
  if(PSO(param,pb).error<pb.epsilon) f=f+1; 
}
f=1./f; // Minimise the inverse of the success rate

break;
		
	}
//------------------
f=fabs(f-objective);  
if(f<errMin) errMin=f; // For information
if(f>errMax) {if(f<infinity) errMax=f; else errMax=infinity;} // For information

return f;
/*
switch(function)
{
	default:
	return  f;
	case 102: // Rosenbrock
	return f/220793652865; // 
}
*/ 
}

//==========================================================
struct fitness constraint(struct position x, int functCode)
{
	// ff[0] is defined in perf()
	// Variables specific to Coil compressing spring
	static double	Fmax=1000.0;
	static double	Fp=300;
	double Cf;
	double K;
	double sp;
	double lf;

	static double	S=189000.0;
	static double	lmax=14.0;
	static double	spm=6.0;
	static double	sw=1.25;
	static double	G=11500000;
	struct fitness ff={0};
	ff.size=1; // Default value

	switch(functCode)
	{
		case 21: // Compression Spring
			Cf=1+0.75*x.x[2]/(x.x[1]-x.x[2]) + 0.615*x.x[2]/x.x[1];
			K=0.125*G*pow(x.x[2],4)/(x.x[0]*x.x[1]*x.x[1]*x.x[1]);
			sp=Fp/K;
			lf=Fmax/K + 1.05*(x.x[0]+2)*x.x[2];

			ff.f[1]=8*Cf*Fmax*x.x[1]/(pi*x.x[2]*x.x[2]*x.x[2]) -S;
			ff.f[2]=lf-lmax;
			ff.f[3]=sp-spm;			
			ff.f[4]=sw- (Fmax-Fp)/K;
			break;
			
		case 25: // Pressure vessel
			ff.f[1]=0.0193*x.x[2]-x.x[0];
			ff.f[2]=0.00954*x.x[2]-x.x[1];
			ff.f[3]=750*1728-pi*x.x[2]*x.x[2]*(x.x[3]+(4.0/3)*x.x[2]); 
			break;	

	}

	return ff;
	}
