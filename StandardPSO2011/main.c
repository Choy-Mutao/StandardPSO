/*
Standard PSO 2011 (from the Particle Swarm Central http://particleswarm.info)
+ some options. In particular
List Based Optimiser 
Mersenne RNG
Quasi-random numbers

 Contact for remarks, suggestions etc.:
Maurice.Clerc@WriteMe.com

	For more details, see ReadMe.txt		
*/

#include "main.h"

// =================================================
int main () 
{ 
	struct position bestBest; // Best position over all runs
	int d;			// Current dimension
	double D;
	double error;			// Current error
	double errorMean;		// Average error
	double errorMin;		// Best result over all runs
	double errorMeanBest[R_max]; 
	double evalMean;		// Mean number of evaluations
	int functionCode;
	int func[funcMax]; // List of functions (codes) to optimise
	int indFunc;

	int nbFunc;
	int nFailure;		// Number of unsuccessful runs
	double logProgressMean;
	struct param param;
	struct problem pb; 
int randCase;
	int run, runMax; 
	struct result result; 

	time_t seconds; 
	
	int scanNb;
	double Smean;
	double success[funcMax];
	double successRate;
	int t;
	double variance;
	float z;
	double zz;

	E = exp ((long double) 1);  
	pi = acos ((long double) -1); 
	errMax=0;
	nbRand=0;
	
	// Files
	f_run = fopen ("f_run.txt", "w");  
	f_synth = fopen ("f_synth.txt", "w");
	f_trace = fopen ("f_trace.txt", "w"); // For information


	//------------------------------------------------ PARAMETERS
			// Bells and Whistles
				// Not really part of the standard
				// May improve the performance (not always)
				// Sometimes not well mathematically founded (rules of thumbs)
				// * => suggested value
				
	param.BW[0]=0; 	//	0 => same swarm size for each run 	
									//	1 => random swarm size around the given mean
									
	param.BW[1]=0; 	//*	0 => when P=G, use the "standard" method	
									// 	1 => when P=G, a specific probabilistic method 
									//	2 => when P=G, a more conservative method
									//  3 => when P=G, just look around
									// 4 =>  different weights for X, P, G  (TEST)
									
	param.BW[2]=0;	// Randomness options
	                // -2nn => Truncated KISS (simulated). 
									//			nn is the number of bits you use to define
									//			each random number/ Example: -207 for 7 bits
	                // -1 => "native" rand() of the C language
									//* 0 => pseudo-random number generator KISS
									//* 10 => pseudo-random number generator Mersenne 64 bits
									// 1 => quasi-random numbers for initialisation, Sobol sequences
									//      KISS after that
									// 2 => quasi-random numbers for initialisation, Halton sequences
									//      KISS after that
									// 3nn => Read on a list of bits (f_rand_bin.txt). 
									//      Normally coming from a "true" (physical) random number generator
									//      (quantic system, atmospheric noise ...)
									//			nn is the number of bits you use to define
									//			each random number/ Example: 307 for 7 bits
									//      Warning: nn must be >=2
									// 4 => Read on a list (f_rand_quasi.txt)
	f_rand_bin= fopen ("f_rand_bin.txt","r"); // A sequence of bits, if BW[2]=3
	f_rand= fopen ("Ltest.txt","r"); //A list of real numbers in ]0,1], if BW[2]=4
									
	param.BW[3]=0;	// 1 => random numbering of the particles before each iteration
									// 0 => always the same loop "particle 0 to S-1"

	//--------
	
	param.confin=0; 	// 0 => keep inside the search space (supposed to be a D-rectangle)
										// 1 => no confinement 
										//   WARNING: may be very slow (and bad) for discrete problems
										
	param.distrib=0; // -1 => uniform in the hypersphere
										//* 0 => in the hypersphere, uniform along the radius 
										// 			(and therefore NOT uniform in the sphere)							
										// 1 => Gaussian (Box-Muller method). Warning: infinite loop possible
										// 2 =>	Gaussian (CMS method)
										// 3 => Other stable (CMS, experimental parameters)
										// 4 => Slash distribution (Gaussian BM/Gaussian BM)
				// Useful only if param.distrib>0;
				param.mean=0.5; //Default: 0.5. For some functions 0 is better, though 
												//	Example: shifted Rosenbrock (code 102)
				param.sigma=1./12; // Default: 1./12 (standard deviation of U(0,1))	
		// WARNING: the CMS method may not work with randomness option >=2
	if(param.BW[2]>=2) param.distrib=	0;

	Smean=40; //Swarm size or Mean swarm size (if BW[0]=1). Suggested: 40 

	param.K=3; 	// Parameter to compute the probability p for a particle to be an
							// external informant. You may also directly define p (see below),
							// but K is about the mean number of the informants of a particle. 
							// Default: 3

	// Confidence coefficients. Default:
	param.w = 1. / (2 * log ((double) 2)); // 0.721
	param.c = 0.5 + log ((double) 2); // 1.193
	param.topology = 0; // 0 => information links as in SPSO 2007 (quasi-random)
											// 1 => variable random ring (EXPERIMENTAL)
											
  //-------------------------------------------------- False randomnesses
  switch (param.BW[2]) //
  {
    default: break;
    case 4: // Prepare a list of false random number, read on a file
    t=0;
    readRand:
    scanNb=fscanf(f_rand,"%f",&z);
    if(scanNb!=EOF) 
    {
      randNumber[t]=z;
      t=t+1;
    goto readRand;
    } 
    nCycleMax=t; 
    printf("\n%i false random numbers read on a file",nCycleMax);

    break; 
  }
	// ----------------------------------------------- PROBLEM
	param.trace=0; // If >0 more information is displayed/saved (f_trace.txt)
	// Functions to optimise
	nbFunc=1; // Number of functions
	func[0]=17; // 4
	func[1]=11;  // 11
	func[2]=15;  // 15
	func[3]=17;  // 17
	func[4]=18; // 18
	func[5]=20;
	func[6]=21;
	func[7]=100;
	func[8]=102;
	func[9]=103;
	func[10]=104;
	func[11]=105;
	func[12]=106;	
 
	/* (see problemDef( ) for precise definitions)
		-1  Constant. For test of biases
	0 Parabola (Sphere)
	1 Griewank
	2 Rosenbrock (Banana)
	3 Rastrigin
	4 Tripod (dimension 2)
	5 Ackley
	6 Schwefel
	7 Schwefel 1.2
	8 Schwefel 2.2
	9 Neumaier 3
	10 G3
	11 Network optimisation (Warning: see problemDef() and also perf() for
			problem elements (number of BTS and BSC)
	12 Schwefel
	13 2D Goldstein-Price
	14 Schaffer f6
	15 Step	
	16 Schwefel 2.21
	17 Lennard-Jones
	18 Gear train
	19 Sine_sine function
	20 Perm function
	21 Compression Spring
	22 Cellular phone (2D) 
	23 Penalized
	24 Repulsion
	25 Pressure Vessel (penalty method)
	26 Ellipsoidal
	27 Quadric
	28 Frequency modulation sound parameter identification

	CEC 2005 benchmark  (no more than 30D. See cec2005data.c)
	100 F1 (shifted Parabola/Sphere) 
	102 F6 (shifted Rosenbrock) 
	103 F9 (shifted Rastrigin) 
	104 F2 Schwefel 
	105 F7 Griewank  (NOT rotated)
	106 F8 Ackley  (NOT rotated)
	107 F4 Schwefel + noise
	
	999 for tests
 

*/ 

	runMax = 100; // Numbers of runs
	if (runMax > R_max) 
	{
		runMax = R_max;
		printf("\nWARNING. I can perform only %i runs. See R_max in main.h",R_max);
	}

	for(indFunc=0;indFunc<nbFunc;indFunc++) // Loop on problems
	{
		functionCode =func[indFunc];

		// Some information
		printf ("\n Function %i ", functionCode);

		// Define the problem
		pb=problemDef(functionCode);
		if(pb.SS.D>DMax) ERROR ("Can't solve it. You should increase DMax");			

		// ----------------------------------------------- RUNS	
		errorMean = 0;	    
		evalMean = 0;	    
		nFailure = 0;	
		D=pb.SS.D;
		randCase=param.BW[2];
		if(randCase>300) {nBit=randCase-300; randCase=3;} // Decode the number of bits
		if(randCase<-200) {nBit=-randCase-200; randCase=-2;}

		switch(randCase)
		{
		default:
		 break;
		 
		case 0:
		seed_rand_kiss(1294404794); // Initialise the RNG KISS for reproducible results
		break;
		
		case 10: // Mersenne 64 bits
		init_genrand64(1294404794);
		//init_genrand64(1234567890);		
		break;
		
		case -2: // Truncated KISS (simulated)
		rMax=pow(2,nBit)-1;
		break;

		case 3: // The file is a string of bits
    //nBit=3; 
		rMax=pow(2,nBit)-1; // For conversion of a nBit string into a number in [0,1]
		break;
		
		case 4: // The file directly contains the numbers
    nCycle=0;	
    break;	
		}
		
		
randRank=0; randChaos=0.02;

	/*
		seconds=time(NULL); // Initialise the RNG KISS more randomly
		 printf("\n time %ld",seconds);
		seed_rand_kiss(time(NULL)); 
*/
	switch(randCase) // "Warm up" the RNG for pseudo-random numbers
	{
		default:
		for (t=0;t<10000;t++) zz=alea(0,1,randCase); 
		break;
		
		case 3:
		case 4:
		break;
	}
//nCycle=4; 
		for (run = 0; run < runMax; run++)  
		{
			if(param.BW[0]==0) param.S=Smean; // Constant swarm size
			else // Random swarm size "around" the mean value
			param.S=(int)(0.5*(0.5+alea(Smean-D/2,Smean+D/2,0)+alea(Smean-D/2,Smean+D/2,0))); 
			
			param.p=1-pow(1-1./param.S,param.K); 
//printf("\n p %f",param.p);      
			// (for a "global best" PSO, directly set param.p=1)
			
			printf("\n Swarm size %i", param.S);
			result = PSO (param, pb);
			error = result.error;

			if (error > pb.epsilon) // Failure
				nFailure = nFailure + 1;

			if(pb.SS.normalise>0)
			{
				for(d=0;d<pb.SS.D;d++)
					result.SW.P[result.SW.best].x[d]=
					pb.SS.min[d]+(pb.SS.max[d]-pb.SS.min[d])*result.SW.P[result.SW.best].x[d]
					/pb.SS.normalise;
			}		

			// Memorize the best (useful if more than one run)
			if(run==0) bestBest=result.SW.P[result.SW.best];
			else
				if(error<bestBest.f) bestBest=result.SW.P[result.SW.best];

			// Result display
			errorMean=errorMean+error;
			printf ("\nRun %i. S %i,  Eval %f. Error %e ", run+1, param.S, result.nEval, error);
			printf(" Mean %e",errorMean/(run+1));
			zz=100*(1-(double)nFailure/(run+1));
			printf("  Success  %.2f%%",zz);

			// Best position display
			//for (d=0;d<pb.SS.D;d++) printf(" %f",result.SW.P[result.SW.best].x[d]);

			// Save result
					fprintf( f_run, "\n%i %.1f %.0f %e %e ", run+1, zz, result.nEval,  error, errorMean/(run+1) );
			// Save best position
					for ( d = 0; d < pb.SS.D; d++ ) fprintf( f_run, " %f",  result.SW.P[result.SW.best].x[d] );

			// Compute/save some statistical information
			if (run == 0)
				errorMin = error;
			else if (error < errorMin)
				errorMin = error;
				
			evalMean = evalMean + result.nEval;	
			errorMeanBest[run] = error;
			logProgressMean  = logProgressMean - log(error);		
		}		// End loop on "run"

		// ---------------------END 
		
		// Display some statistical information
		evalMean = evalMean / (double) runMax;   
		errorMean = errorMean / (double) runMax;
		logProgressMean = logProgressMean/(double) runMax;

		printf ("\n Eval. (mean)= %f", evalMean);	
		printf ("\n Error (mean) = %e", errorMean);

		// Variance
		variance = 0;

		for (run = 0; run < runMax; run++)
			variance = variance + pow (errorMeanBest[run] - errorMean, 2);

		variance = sqrt (variance / runMax);	    
		printf ("\n Std. dev. %e", variance); 
		printf("\n Log_progress (mean) = %f", logProgressMean);	

		// Success rate and minimum value
		printf("\n Failure(s) %i  ",nFailure);

		successRate = 100 * (1 - nFailure / (double) runMax);			
		printf ("Success rate = %.2f%%", successRate);
		success[indFunc]=successRate;

		printf ("\n Best min value = %1.20e", errorMin);
		printf ("\nPosition of the optimum: ");
		for (d=0;d<pb.SS.D;d++) printf(" %.20f",bestBest.x[d]);
		

		// Save	
		fprintf (f_synth, "%i %i %.0f %e %e %f %f", functionCode, pb.SS.D,successRate,
		    errorMean, variance, evalMean,bestBest.f);
		for (d=0;d<pb.SS.D;d++) fprintf(f_synth, " %1.20e",bestBest.x[d]);  
		fprintf(f_synth,"\n"); 

// Specific save for Repulsive problem
if(pb.function==24)
{		
	for (d=0;d<pb.SS.D-1;d=d+2) 
	{
	fprintf(f_synth, " %1.20e %1.20e",bestBest.x[d],bestBest.x[d+1]);  
	fprintf(f_synth,"\n"); 
	}	
}	     	   

} // End "for ind[..."

printf("\n errMax : %f",errMax);

		// Repeat informations
	printf("\n---------");
	printf("\n Function(s):");
		for(indFunc=0;indFunc<nbFunc;indFunc++) // Loop on problems
		{
			functionCode =func[indFunc];
			printf(" %i",functionCode);
		}
		printf("\n Confinement: "); 
		if(param.confin==0) printf("YES"); else printf("NO");
		printf("\n Distribution: ");
		switch(param.distrib)
		{
			case 0: printf(" uniform"); break;
			case 1: printf(" Gaussian (%f,%f), Box-Muller",param.mean,param.sigma); break;
			case 2: printf(" Gaussian (%f,%f), CMS",param.mean,param.sigma); break;
			case 3: printf(" Stable (%f,%f)",param.mean,param.sigma); break;
			case 4: printf(" Slash (%f,%f)",param.mean,param.sigma); break;
		}

		printf("\n BW = (%i, %i, %i, %i)",param.BW[0],param.BW[1],
		       param.BW[2],param.BW[3]);
		printf("\n Swarm size: ");
		if(param.BW[0]==0) printf("%i",(int)Smean); else printf(" mean %i",(int)Smean);
		printf("\n K = %i",param.K);
		printf("\n w = %f",param.w);
		printf("\n c = %f",param.c);
		printf("\n %e random numbers have been used",nbRand);
		fprintf(f_run,"\nnbRand %e",nbRand);
	return 0; // End of main program
}
// ===============================================================
#include "alea.c"
#include "perf.c"
#include "problemDef.c"
#include "PSO.c"
#include "tools.c"


