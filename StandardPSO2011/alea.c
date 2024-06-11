
double alea (double a, double b, int option) 
{				// random number (uniform distribution) in  [a b]
	int close;
	double r;
	int i;
	int ir;
	char list[nBit]; 
	int listBit[nBit];

	double pw;
	int scanNb;
	float zz;

	switch(option)
	{
		default:
			case 0: // KISS	
			r=(double)rand_kiss()/RAND_MAX_KISS;
			break;
			
			case 10:  // Mersenne
			r=genrand64_real2();
			break;

      case -2: // Truncated KISS (simulated)
      r=(double)alea_integer(0,rMax,0)/rMax;
      break;
      
      case -1: // Standard C	
			r=(double)rand()/RAND_MAX;
			break;

			case 3: // Read on a file (list of 0 and 1)
			// Attempt to read in nBit characters 
			scanNb = fread( list, sizeof( char ), nBit, f_rand_bin );
			if(scanNb!=nBit || nCycle>=nCycleMax)
		{
			nCycle=0;
			close=fclose(f_rand_bin);
			if(close!=0) ERROR("\n Can't close f_rand_bin");
			if( (f_rand_bin = fopen( "f_rand_bin.txt", "r+t" )) != NULL ) 
				fread( list, sizeof( char ), nBit, f_rand_bin );
			else	ERROR("\n Can't reopen f_rand_bin");
		}

			for(i=0;i<nBit;i++) listBit[i]=list[i]-48;

			// Convert into double
			r=listBit[0];
			if(nBit>1)
			{
				pw=1;
				for(i=1;i<nBit;i++)
				{
					pw=pw*2;
					if(listBit[i]==0) continue;
					r=r+pw;
				}
				r=r/rMax;
			}
				nCycle++; 
			break;
			
			case 4: // False random numbers have been read on a file
	 if(nCycle>=nCycleMax) nCycle=0;
	  r=randNumber[nCycle];
	  nCycle++;
			break;
	}

	nbRand++;
	 r=a+r*(b-a);
	return r; 
}

//==================================================
void aleaIndex(int index[], int S,int option)
{
	int indexTemp[S_max];
	int length;
	int rank;
	int s;
	int t;

	length=S;
	for (s=0;s<S;s++) indexTemp[s]=s; //=index[s];

	for (s=0;s<S;s++)
	{	
			rank=alea_integer(0,length-1,option);
			index[s]=indexTemp[rank];
			if (rank<length-1)	// Compact
			{	
				for (t=rank;t<length-1;t++)
					indexTemp[t]=indexTemp[t+1];
			}					
			length=length-1;
	}
}
// ===========================================================
int alea_integer (int a, int b,int option) 
{				// Integer random number in [a b]
	int ir;
	double r;

	//r = alea (0, 1);
	//ir = (int) (a + r * (b + 1 - a));
	//if (ir > b)	ir = b;

	if(a==b) return a;
	r=alea((double)a,(double)b+1,option); ir=floor(r);

	if(ir>b) ir=b;

	return ir;  
}
// ===========================================================
double alea_stable (double alpha, double beta, double nu, double delta,int option)
{	// CMS algorithm (J.M. Chambers, C.L. Mallows and B.W. Stuck)
	// alpha 	:= 	stability parameter. Tail thickness (the higher the thinner)
	//						 Must be in ]0,2]. 
	//						For normal distribution, alpha=2
	// beta		:=	skweness.  0 => symmetric distribution
	// nu			:=	dispersion. Must be >0. For normal distribution, nu=standard dev.
	// delta	:=	mean (mesure of centrality)
	// WARNING: doesn't work if a random number drawn in [0,1] is precisely 0 or 1
	double betaPrime; 
	double d;
	double eps;
	double kAlpha;
	double min=zero; //0; // To avoid to have to compute ln(0) ...
	double max=1;
	double phi0, phi;
	double r;
	double s;
	double t1,t2,t3;
	double tau;
	double u;
	double temp;
	double w;
	double z;

	if(alpha<0 || alpha>2)
	{ printf("\n alpha %f ",alpha);
		ERROR("alea_levy. alpha must be in ]0,2]");
	}
	if(nu<0)
	{
		printf("\n nu %f ",nu);
		ERROR("alea_levy. nu must be positive");
	}
	//--------------------------------------------
	// Define k(α) = 1 − |1−α|. Thus if α≤1 then k(α)=α and if α≥1 then k(α)=2−α.   
	if(alpha<1) kAlpha=alpha; else kAlpha=2-alpha;

	// Compute φ0 = −½β(k(α)/α).  
	phi0=0.5*beta*(kAlpha/alpha);

	temp=tan(alpha*phi0);

	// Transform β to β' by β' = β if α=1 and β'= −tan(½π(1−α))tan(αφ0) otherwise.

	if(fabs(alpha-1)<zero) betaPrime=beta;
	else
		betaPrime=-tan(0.5*pi*(1-alpha))*temp;

	// Generate a random variable u uniformly distributed on the interval [0, 1] and compute φ = π(u−½).
	u= alea(min,max,option);			
	phi=pi*(u-0.5);

	//Compute ε=1−α and then τ = −εtan(αφ0) 		
	eps=1-alpha; 
	tau=-eps*temp;

	//Compute tan(½φ), tan(½εφ) and tan(½εφ)/(½εφ). 
	t1=tan(0.5*phi);
	t2=tan(0.5*eps*phi);
	t3=2*t2/(eps*phi);

	// Generate a random variable v which has a uniform distribution 
	//	on the interval [0,1] and then compute w=−ln(v)
	w=-log(alea(min,max,option));

	// Compute z = (cos(εφ)−tan(αφ0)sin(εφ)/(wcos(φ))
	z=cos(eps*phi)-tan(alpha*phi0)*sin(eps*phi)/(w*cos(phi));

	// Compute d = z^ε/α /ε 
	temp=pow(z,eps/alpha);
	d=temp/eps;

	// Compute s = tan(αφ0) + z^ε/α (sin(αφ)−tan(αφ0)cos(αφ))/cos(φ) 
	s=tan(alpha*phi0) + temp*(sin(alpha*phi) - tan(alpha*phi0)*cos(alpha*phi))/cos(phi);	

	// Multiply by the dispersion, and add the mean
	r=s*nu + delta;

	return r;
}
// ===========================================================
double alea_normal (double mean, double std_dev,int option) 
{ 
	// Use the polar form of the Box-Muller transformation to obtain a pseudo
	// random number from a Gaussian distribution 

	double x1, x2, w, y1; 
	// WARNING. This method is valid only if alea (...) defines a (more or less)
	// uniform distribution. In particular, if "option" is so that the numbers are
	// read on a list, whose distribution is, say, Gaussian, then w will never be
	// smaller than 1 

int count=0;
	do  
	{
		x1 = 2.0 * alea (0, 1,option) - 1.0;
		x2 = 2.0 * alea (0, 1,option) - 1.0;
		w = x1 * x1 + x2 * x2; 
		count++;
		if(count>10000) //RAND_MAX)
		{
		if(option==4) printf("\n Random numbers read on a list");
		ERROR("alea_normal. Probable non uniform distribution for alea(...)");
		}
	}
	while (w >= 1.0); // WARNING. Infinite loop possible.  w may be _never_ < 1 !!
	  

	w = sqrt (-2.0 * log (w) / w);
	y1 = x1 * w;

	if(alea(0,1,option)<0.5) y1=-y1; 
	y1 = y1 * std_dev + mean;
	return y1;  
}

// ===========================================================
struct vector alea_sphere(int D, double radius, int distrib, double mean,
    double sigma,int option)
{
	/*  ******* Random point in a hypersphere ********
	 Maurice Clerc 2003-07-11
	 Last update: 2011-01-01

	 Put  a random point inside the hypersphere S(center 0, radius 1), 
	 or on its surface
		 */

	int 	j;
	double   length;
	double      pw;
	double      r;
	struct	vector	x;

	x.size=D;
	pw=1./(double)D;

	// ----------------------------------- Step 1.  Direction
	length=0;
	for (j=0;j<D;j++)
	{
		// Here a Gaussian distribution is needed
		//if(distrib<2 && option!=4) x.v[j]=alea_normal(0,1,option); 
		if(option!=4) 
		x.v[j]=alea_normal(0,1,option);// Gaussian (Box-Muller method)
		else // When the "random" numbers are read on a file, the BM method 
		      // may loop infinitly
		x.v[j]=alea_stable(2,0,1,0,option); // Gaussian (CMS method)
		length=length+  x.v[j]*x.v[j];
	}

	length=sqrt(length);
	//----------------------------------- Step 2. Random radius

	switch(distrib) 
	{
		default: // Uniform distribution
							// Note that the distribution in the hypersphere is NOT uniform
		case 0:
			r=alea(0,1,option); 
		break;

		case -1: // So that the final distribution is uniform in the sphere
		  r=alea(0,1,option); 
			r=pow(r,pw);
		break;
		
		case 1: //Gaussian (Box-Muller)
			r=fabs(alea_normal(mean,sigma,option)); 
		break;

		case 2: //Gaussian (CMS)
			r=fabs(alea_stable(2,0,sigma,mean,option)); 
		break; 

		case 3: //TEST
			//r=fabs(alea_stable(0.2,0,sigma,alea(0,1))); 
			r=fabs(alea_stable(0.2,0,sigma,0.5*(alea(0, 1,option)+alea(0,1,option)),option)); 
		break; 

		case 4: // Slash distribution
			r=fabs(alea_normal(mean,sigma,option))/fabs(alea_normal(0,1,option));
		break;

		case 99: // Constant radius, for some specific uses (random ON the sphere)
			r=1; 
		break;
	} 


	for (j=0;j<D;j++)
	{
		x.v[j]=radius*r*x.v[j]/length;
	}
	return x;
}
//==========================================================================    
struct vectorList quasiRand(int D, int nRand, int option)
{   
	/*
	 Generate nRand vectors of size D that are the coordinates of
	 nRand quasi-random points in [0,1]^D

		 For C under Linux, you need to use the GSL library:
#include <gsl/gsl_qrng.h> // Do not forget to link to gsl and gslcblast

		 */

	int i;	
	void *qrng_q;    
	struct vectorList qRand;

	switch(option)
	{
		case 1: // Sobol
			if(D>40) 
		{
			printf("\nThe embedded Sobol sequences generator can not be used for dimensions greater than 40");
			printf("\nYou should use Halton sequences (for dimensions up to 1229)");
			ERROR("\n I stop here");

		}
		qrng_q=gsl_qrng_alloc (sobol, D);
		break;

		default: // Halton
			if(D>1229) 
		{printf("\nThe embedded Halton sequences generator can not be used for dimensions greater than 1229");
			printf("\nSorry");
			ERROR("\n I stop here");
		}
		qrng_q=gsl_qrng_alloc (halton, D);
		break;
	}

	gsl_qrng * q = qrng_q;

	for (i = 0; i < nRand; i++)
	{
		gsl_qrng_get (q, qRand.V[i].v);       
	}

	gsl_qrng_free (q); 

	return qRand;
}
//================================================
#include "KISS.c"
#include "mersenne.c"
