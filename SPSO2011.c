#include "spso2011.h"

static int SPSO2011()
{
	printf_s("Structure SPSO %i", 2011);

	struct position bestBest; // Best position over all runs
	bestBest = *(struct position*)malloc(sizeof(struct position));
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
	double logProgressMean = 0;
	struct param param;
	param = *(struct param*)malloc(sizeof(struct param));
	struct problem pb;
	pb = *(struct problem*)malloc(sizeof(struct problem));
	int randCase;
	int run, runMax;
	struct result* result;
	//result = *(struct result*)malloc(sizeof(struct result));

	time_t seconds;

	int scanNb;
	double Smean;
	double success[funcMax];
	double successRate;
	int t;
	double variance;
	float z;
	double zz;

	E = exp((long double)1);
	pi = acos((long double)-1);
	errMax = 0;
	nbRand = 0;

	// Files
	f_run = fopen("f_run.txt", "w");
	f_synth = fopen("f_synth.txt", "w");
	f_trace = fopen("f_trace.txt", "w"); // For information


	//------------------------------------------------ PARAMETERS
			// Bells and Whistles
				// Not really part of the standard
				// May improve the performance (not always)
				// Sometimes not well mathematically founded (rules of thumbs)
				// * => suggested value

	param.BW[0] = 0; 	//	0 => same swarm size for each run 	
	//	1 => random swarm size around the given mean

	param.BW[1] = 0; 	//*	0 => when P=G, use the "standard" method	
	// 	1 => when P=G, a specific probabilistic method 
	//	2 => when P=G, a more conservative method
	//  3 => when P=G, just look around
	// 4 =>  different weights for X, P, G  (TEST)

	param.BW[2] = 0;	// Randomness options
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

	f_rand_bin = fopen("f_rand_bin.txt", "r"); // A sequence of bits, if BW[2]=3
	f_rand = fopen("Ltest.txt", "r"); //A list of real numbers in ]0,1], if BW[2]=4

	param.BW[3] = 0;	// 1 => random numbering of the particles before each iteration
	// 0 => always the same loop "particle 0 to S-1"

	//--------

	param.confin = 0; 	// 0 => keep inside the search space (supposed to be a D-rectangle)
	// 1 => no confinement 
	//   WARNING: may be very slow (and bad) for discrete problems

	param.distrib = 0; // -1 => uniform in the hypersphere
	//* 0 => in the hypersphere, uniform along the radius 
	// 			(and therefore NOT uniform in the sphere)							
	// 1 => Gaussian (Box-Muller method). Warning: infinite loop possible
	// 2 =>	Gaussian (CMS method)
	// 3 => Other stable (CMS, experimental parameters)
	// 4 => Slash distribution (Gaussian BM/Gaussian BM)
// Useful only if param.distrib>0;
	param.mean = 0.5; //Default: 0.5. For some functions 0 is better, though 
	//	Example: shifted Rosenbrock (code 102)
	param.sigma = 1. / 12; // Default: 1./12 (standard deviation of U(0,1))	
	// WARNING: the CMS method may not work with randomness option >=2
	if (param.BW[2] >= 2) param.distrib = 0;

	Smean = 40; //Swarm size or Mean swarm size (if BW[0]=1). Suggested: 40 

	param.K = 3; 	// Parameter to compute the probability p for a particle to be an
	// external informant. You may also directly define p (see below),
	// but K is about the mean number of the informants of a particle. 
	// Default: 3

	// Confidence coefficients. Default:
	param.w = 1. / (2 * log((double)2)); // 0.721
	param.c = 0.5 + log((double)2); // 1.193
	param.topology = 0; // 0 => information links as in SPSO 2007 (quasi-random)
	// 1 => variable random ring (EXPERIMENTAL)

	//-------------------------------------------------- False randomnesses
	switch (param.BW[2])
	{
	case 4: // Prepare a list of false random number, read on a file
		t = 0;
	readRand:
		scanNb = fscanf(f_rand, "%f", &z);
		if (scanNb != EOF)
		{
			randNumber[t] = z;
			t = t + 1;
			goto readRand;
		}
		nCycleMax = t;
		printf("\n%i false random numbers read on a file", nCycleMax);
		break;
	default:
		break;
	}

	// ----------------------------------------------- PROBLEM
	param.trace = 0; // If >0 more information is displayed/saved(f_trace.txt)
	// Functions to optimise
	nbFunc = 1; // Number of functions
	func[0] = 17; // 4
	func[1] = 11;  // 11
	func[2] = 15;  // 15
	func[3] = 17;  // 17
	func[4] = 18; // 18
	func[5] = 20;
	func[6] = 21;
	func[7] = 100;
	func[8] = 102;
	func[9] = 103;
	func[10] = 104;
	func[11] = 105;
	func[12] = 106;

	runMax = 100; // Numbers of runs
	if (runMax > R_max)
	{
		runMax = R_max;
		printf("\nWARNING. I can perform only %i runs. See R_max in main.h", R_max);
	}

	for (indFunc = 0; indFunc < nbFunc; indFunc++) // Loop on problems
	{
		functionCode = func[indFunc];

		// Some information
		printf("\n Function %i ", functionCode);

		// Define the problem
		pb = problemDef(functionCode);
		if (pb.SS.D > DMax) ERROR("Can't solve it. You should increase DMax");

		// ----------------------------------------------- RUNS	
		errorMean = 0;
		evalMean = 0;
		nFailure = 0;
		D = pb.SS.D;
		randCase = param.BW[2];
		if (randCase > 300) { nBit = randCase - 300; randCase = 3;} // Decode the number of bits
		if (randCase < -200) { nBit = -randCase - 200;randCase = -2;}

		switch (randCase)
		{
		case 0:
			seed_rand_kiss(1294404794); // Initialise the RNG KISS for reproducible results
			break;
		case 10:						// Mersenne 64 bits
			init_genrand64(1294404794);
			//init_genrand64(1234567890);
			break;
		case -2:						// Truncated KISS (simulated)
			rMax = pow(2, nBit) - 1;
			break;
		case 3:							// The file is a string of bits
			// nBit = 3;
			rMax = pow(2, nBit) - 1;	// For conversion of a nBit string int a number in [0,1]
			break;
		case 4:							// The file directly contains the numbers
			nCycle = 0;
		default:
			break;
		}

		randRank = 0; randChaos = 0.02;

		switch (randCase)				// "Warm up" the RNG for pseudo-random numbers
		{
		case 3:
		case 4:
		default:
			for (int t = 0; t < 10000; t++)
				zz = alea(0, 1, randCase);
			break;
		}

		// nCycle = 4;
		for (int run = 0; run < runMax; run++)
		{
			if (param.BW[0] == 0) param.S = Smean; // Constant swarm size
			else
				param.S = (int)(0.5 * (0.5 + alea(Smean - D / 2, Smean + D / 2, 0) + alea(Smean - D / 2, Smean + D / 2, 0)));

			param.p = 1 - pow(1 - 1. / param.S, param.K);
			printf("\n p = %f", param.p); //(for a "global best" PSO, directly set param.p = 1)

			printf("\n Swarm size %i", param.S);
			result = PSO(&param, &pb);
			error = result->error;

			if (error > pb.epsilon) // Failure
				nFailure += 1;

			if (pb.SS.normalise > 0)
			{
				for ( d = 0; d < pb.SS.D; d++)
				{
					result->SW.P[result->SW.best].x[d] =
						pb.SS.min[d] + (pb.SS.max[d] - pb.SS.min[d]) * result->SW.P[result->SW.best].x[d]/ pb.SS.normalise;
				}
			}

			// Memorize the best (useful if more than one run)
			if (run == 0) bestBest = result->SW.P[result->SW.best];
			else
				if (error < bestBest.f) bestBest = result->SW.P[result->SW.best];

			// Result display
			errorMean = errorMean + error;
			printf("\nRun %i. S %i,  Eval %f. Error %e ", run + 1, param.S, result->nEval, error);
			printf(" Mean %e", errorMean / (run + 1));
			zz = 100 * (1 - (double)nFailure / (run + 1));
			printf("  Success  %.2f%%", zz);

			// Best position display
			for (d=0;d<pb.SS.D;d++) printf(" %f",result->SW.P[result->SW.best].x[d]);

			// Save result
			fprintf(f_run, "\n%i %.1f %.0f %e %e ", run + 1, zz, result->nEval, error, errorMean / (run + 1));
			// Save best position
			for (d = 0; d < pb.SS.D; d++) fprintf(f_run, " %f", result->SW.P[result->SW.best].x[d]);

			// Compute/save some statistical information
			if (run == 0)
				errorMin = error;
			else if (error < errorMin)
				errorMin = error;

			evalMean = evalMean + result->nEval;
			errorMeanBest[run] = error;
			logProgressMean = logProgressMean - log(error);
		}		// End loop on "run"

		// ---------------------END
		
		// Display some statistical information
		evalMean = evalMean / (double)runMax;
		errorMean = errorMean / (double)runMax;
		logProgressMean = logProgressMean / (double)runMax;

		printf("\n Eval. (mean)= %f", evalMean);
		printf("\n Error (mean) = %e", errorMean);
		// Variance
		variance = 0;

		for (run = 0; run < runMax; run++)
			variance = variance + pow(errorMeanBest[run] - errorMean, 2);

		variance = sqrt(variance / runMax);
		printf("\n Std. dev. %e", variance);
		printf("\n Log_progress (mean) = %f", logProgressMean);

		// Success rate and minimum value
		printf("\n Failure(s) %i  ", nFailure);

		successRate = 100 * (1 - nFailure / (double)runMax);
		printf("Success rate = %.2f%%", successRate);
		success[indFunc] = successRate;

		printf("\n Best min value = %1.20e", errorMin);
		printf("\nPosition of the optimum: ");
		for (d = 0; d < pb.SS.D; d++) printf(" %.20f", bestBest.x[d]);


		// Save	
		fprintf(f_synth, "%i %i %.0f %e %e %f %f", functionCode, pb.SS.D, successRate,
			errorMean, variance, evalMean, bestBest.f);
		for (d = 0; d < pb.SS.D; d++) fprintf(f_synth, " %1.20e", bestBest.x[d]);
		fprintf(f_synth, "\n");

		// Specific save for Repulsive problem
		if (pb.function == 24)
		{
			for (d = 0; d < pb.SS.D - 1; d = d + 2)
			{
				fprintf(f_synth, " %1.20e %1.20e", bestBest.x[d], bestBest.x[d + 1]);
				fprintf(f_synth, "\n");
			}
		}
	}	// End "for ind[..."

	printf("\n errMax : %f", errMax);

	// Repeat informations
	printf("\n---------");
	printf("\n Function(s):");
	for (indFunc = 0; indFunc < nbFunc; indFunc++) // Loop on problems
	{
		functionCode = func[indFunc];
		printf(" %i", functionCode);
	}
	printf("\n Confinement: ");
	if (param.confin == 0) printf("YES"); else printf("NO");
	printf("\n Distribution: ");
	switch (param.distrib)
	{
	case 0: printf(" uniform"); break;
	case 1: printf(" Gaussian (%f,%f), Box-Muller", param.mean, param.sigma); break;
	case 2: printf(" Gaussian (%f,%f), CMS", param.mean, param.sigma); break;
	case 3: printf(" Stable (%f,%f)", param.mean, param.sigma); break;
	case 4: printf(" Slash (%f,%f)", param.mean, param.sigma); break;
	}

	printf("\n BW = (%i, %i, %i, %i)", param.BW[0], param.BW[1],
		param.BW[2], param.BW[3]);
	printf("\n Swarm size: ");
	if (param.BW[0] == 0) printf("%i", (int)Smean); else printf(" mean %i", (int)Smean);
	printf("\n K = %i", param.K);
	printf("\n w = %f", param.w);
	printf("\n c = %f", param.c);
	printf("\n %e random numbers have been used", nbRand);
	fprintf(f_run, "\nnbRand %e", nbRand);
	return 0; // End of main program
}

// ===============================================================

// ========================Probdef.c==============================
struct problem problemDef(int functionCode)
{
	int d;
	struct problem pb;

	int nAtoms; // For Lennard-Jones problem
	static double lennard_jones[14] = { -1, -3, -6, -9.103852, -12.7121, -16.505384,-19.821489,-24.113360,-28.422532,
		-32.77,-37.97,-44.33,-47.84,-52.32 };

	pb.function = functionCode;
	// Default values
	// Can be modified below for each function
	pb.epsilon = 0.00000;	// Acceptable error(default). May be modified below
	pb.objective = 0;		// Objective value (default). May be modified below
	// pb.constraintNb = 0;
	pb.SS.quantisation = 0;	// No quantisation needed (all variables are continuous)
	pb.SS.normalise = 0;	// Set to a value x.
							// x>0 => Normalisation is applied(search space => [0,x]^D)

	// --------------- Search space
	switch (pb.function)
	{
	case 0:					// Parabola
		pb.SS.D = 30;		// Dimension
		
		for (int d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100; // -100
			pb.SS.max[d] = 100;	 // 100
			pb.SS.q.q[d] = 0;	 // Relative quantisation, in [0,1]
		}
		pb.evalMax = 75000;		 // 100000; Max number of evaluations for each run
		pb.epsilon = 0.01;		 // 0.0000001; // 1e-3
		pb.objective = 0;
		break;

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
		pb.epsilon = 0.01; //0.001; //0.15; //0.001;
		pb.objective = 0;
		break;
	case 2:		// Rosenbrock
		pb.SS.D = 30;	// 30

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -30; // -30; 
			pb.SS.max[d] = 30; // 30;			
			pb.SS.q.q[d] = 0;
		}
		pb.epsilon = 100; //0.05; //.0001;		
		pb.evalMax = 75000; //200000; //2.e6;  // 40000 
		pb.objective = 0;
		break;
	case 3:		// Rastrigin
		pb.SS.D = 30; // 10;	

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -5.12;
			pb.SS.max[d] = 5.12;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 75000; //3200; 
		pb.epsilon = 50; //0.001;
		pb.objective = 0;
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
		pb.epsilon = 0.0001;
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
		pb.evalMax = 80000; //30000; 
		pb.epsilon = 0.000;
		pb.objective = 0;


		break;
	case 6: // Schwefel. Min on (A=420.8687, ..., A)
		pb.SS.D = 30;
		//pb.objective=-pb.SS.D*420.8687*sin(sqrt(420.8687));
		pb.objective = -12569.5;
		pb.epsilon = 2569.5;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -500;
			pb.SS.max[d] = 500;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 300000;


		break;
	case 7: // Schwefel 1.2
		pb.SS.D = 40;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 40000;
		break;
	case 8: // Schwefel 2.22
		pb.SS.D = 30;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 100000;
		pb.objective = 0;
		pb.epsilon = 0.0001;
		break;
	case 9: // Neumaier 3
		pb.SS.D = 40;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -pb.SS.D * pb.SS.D;
			pb.SS.max[d] = -pb.SS.min[d];
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 40000;


		break;
	case 10: // G3 (constrained)
		pb.SS.D = 10;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 340000;
		pb.objective = 0;
		pb.epsilon = 1.e-6;
		pb.SS.quantisation = 1;
		break;
	case 11: // Network
		//	btsNb=5; bcsNb=2;
		btsNb = 19; bcsNb = 2;
		pb.SS.D = bcsNb * btsNb + 2 * bcsNb;
		pb.objective = 0;
		for (d = 0; d < bcsNb * btsNb; d++) // Binary representation. 1 means: there is a link
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 1;
		}
		pb.SS.quantisation = 1;
		for (d = bcsNb * btsNb; d < pb.SS.D; d++) // 2D space for the BSC positions
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 20; //15;
			pb.SS.q.q[d] = 0;
		}
		pb.SS.normalise = 1;
		pb.evalMax = 5000;
		pb.objective = 0;
		pb.epsilon = 0;

		break;
	case 12: // Schwefel
		pb.SS.D = 3;
		pb.objective = -418.98288727243369 * pb.SS.D;

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
		pb.objective = 0;

		pb.SS.min[0] = -100;
		pb.SS.max[0] = 100;
		pb.SS.q.q[0] = 0;
		pb.SS.min[1] = -100;
		pb.SS.max[1] = 100;
		pb.SS.q.q[1] = 0;
		pb.evalMax = 720;

		break;
	case 14: // Schaffer f6	 
		pb.SS.D = 2;	// Dimension
		pb.objective = 0;
		pb.epsilon = 0.0001;
		pb.SS.min[0] = -100;
		pb.SS.max[0] = 100;
		pb.SS.q.q[0] = 0;
		pb.SS.min[1] = -100;
		pb.SS.max[1] = 100;
		pb.SS.q.q[1] = 0;

		pb.evalMax = 30000;
		break;
	case 15: // Step
		pb.SS.D = 10;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 2500;

		break;
	case 16: // Schwefel 2.21
		pb.SS.D = 30;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 100000;

		break;
	case 17: // Lennard-Jones
		nAtoms = 6; // in {2, ..., 15}
		pb.SS.D = 3 * nAtoms;
		pb.objective = lennard_jones[nAtoms - 2];
		pb.evalMax = 5000 + 3000 * nAtoms * (nAtoms - 1); // Empirical rule
		pb.epsilon = 1.e-6;

		pb.SS.normalise = 1;
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
		pb.SS.D = 4;
		pb.objective = 2.7e-12;
		pb.epsilon = 1.e-13;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 12;
			pb.SS.max[d] = 60;
			pb.SS.q.q[d] = 1;
		}
		pb.evalMax = 20000;
		pb.SS.quantisation = 1;
		break;
	case 19: // Sine sine function
		pb.SS.D = 10;
		pb.objective = -10;// Arbitrary large negative number
		// Remember that the error is abs(f - objective), though
		// Best known (2010-09: -9.5983769). 
		pb.epsilon = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = pi;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 60000;

		break;
	case 20: // Perm function
		pb.SS.D = 5;
		pb.objective = 0;
		pb.epsilon = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -pb.SS.D;
			pb.SS.max[d] = pb.SS.D;
			pb.SS.q.q[d] = 1;
		}
		pb.evalMax = 10000;
		pb.SS.quantisation = 1;
		break;
	case 21: // Compression spring
		pb.SS.D = 3;

		pb.SS.min[0] = 1; pb.SS.max[0] = 70; pb.SS.q.q[0] = 1; // N
		pb.SS.min[1] = 0.6; pb.SS.max[1] = 3; pb.SS.q.q[1] = 0; // D
		pb.SS.min[2] = 0.207; pb.SS.max[2] = 0.5; pb.SS.q.q[2] = 0.001; // d

		/*
			pb.SS.min[0] = 2; pb.SS.max[0] = 10; pb.SS.q.q[0] = 1; // N
			pb.SS.min[1] = 0.25; pb.SS.max[1] = 1.3; pb.SS.q.q[1] = 0; // D
			pb.SS.min[2] = 0.05; pb.SS.max[2] = 2; pb.SS.q.q[2] = 0.001; // d
		*/

		pb.SS.quantisation = 1;
		pb.SS.normalise = 1;
		pb.evalMax = 20000;
		pb.epsilon = 1.e-10;
		pb.objective = 2.6254214578;
		break;
	case 22:// Cellular phone
		pb.SS.D = 2 * 5; //2*10  2*nb_of_stations

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 100; // Warning: hard coded in cellular_phone.c
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 2000 * pb.SS.D;
		pb.epsilon = 1e-9;
		pb.objective = 0.005530517; // Best known result (2010-01-03)
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
		pb.epsilon = 0;
		pb.objective = 0;

		break;
	case 24:// Repulsion (2D)
		pb.SS.D = 2 * 40; // 2*nb_of_charged_points

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 100; // Warning: hard coded in repulsion.c
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 3000 * pb.SS.D;
		pb.epsilon = 0;
		pb.objective = 0;
		break;

		// -------------------pressure_vessel.c---------------------
	case 25: //  Pressure vessel (penalty method)		
		pb.SS.D = 4;

		pb.SS.min[0] = 1.125;
		pb.SS.max[0] = 12.5;
		pb.SS.q.q[0] = 0.0625;

		pb.SS.min[1] = 0.625;
		pb.SS.max[1] = 12.5;
		pb.SS.q.q[1] = 0.0625;

		pb.SS.min[2] = 0; pb.SS.max[2] = 240;
		pb.SS.q.q[2] = 0;
		pb.SS.min[3] = 0; pb.SS.max[3] = 240;
		pb.SS.q.q[3] = 0;

		pb.objective = 7197.72893;
		pb.SS.quantisation = 1;
		pb.SS.normalise = 1;	 // Very important here.
		pb.evalMax = 50000;
		pb.epsilon = 0; // 0.00001; 
		break;
		// ---------------------------------------------------------

	case 26: // Ellipsoidal
		pb.SS.D = 30;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -30;
			pb.SS.max[d] = 30; // Warning: hard coded in repulsion.c
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 3000;
		pb.epsilon = 50;
		pb.objective = 0;

		break;

		// ----------------------cec2005pb.c------------------------
	case 100:				// CEC 2005 F1
		pb.SS.D = 30;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = pb.SS.D * 10000;
		pb.epsilon = 0.000001;		// Acceptable error
		pb.objective = -450;		// Objective value
		break;
	case 102:				// Rosenbrock. CEC 2005 F6
		pb.SS.D = 10;
		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = pb.SS.D * 10000;
		pb.epsilon = 0.01;		// Acceptable error
		pb.objective = 390;
		break;
	case 103:					// CEC 2005 F9, Rastrigin
		pb.SS.D = 30;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -5.12;
			pb.SS.max[d] = 5.12;
			pb.SS.q.q[d] = 0;
		}
		pb.epsilon = 0.01;		// Acceptable error
		pb.objective = -330;	// Objective value
		pb.evalMax = pb.SS.D * 10000;
		break;
	case 104:// CEC 2005 F2  Schwefel
		break;

	case 107: // F4 Schwefel + noise
		pb.SS.D = 10;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;

		}
		pb.epsilon = 0.00001;	//0.00001 Acceptable error
		pb.objective = -450;       // Objective value
		pb.evalMax = pb.SS.D * 10000;
		break;

	case 105:// CEC 2005 F7  Griewank (NON rotated)
		pb.SS.D = 10;	 // 10 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -600;
			pb.SS.max[d] = 600;
			pb.SS.q.q[d] = 0;

		}
		pb.epsilon = 0.01;	//0.01 Acceptable error
		pb.objective = -180;       // Objective value
		pb.evalMax = pb.SS.D * 10000;
		break;

	case 106:// CEC 2005 F8 Ackley (NON rotated)
		pb.SS.D = 10; // 10;	 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -32;
			pb.SS.max[d] = 32;
			pb.SS.q.q[d] = 0;
		}
		pb.epsilon = 0.0001;	//.0001 Acceptable error
		pb.objective = -140;       // Objective value
		pb.evalMax = pb.SS.D * 10000;
		break;

		// -------------------------------------------------------
	case 999:// for tests
		pb.SS.D = 20;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -5;
			pb.SS.max[d] = 5;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 10000;
		pb.epsilon = 0.0001;
		pb.objective = -78.3323;

		break;

		pb.SS.D = 2;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -512;
			pb.SS.max[d] = 512;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 10000;
		pb.epsilon = 0;
		pb.objective = -511.7;

		break;
		pb.SS.D = 2;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 2000;
		pb.epsilon = 0;
		pb.objective = -1.031628;
		break;
	default:
		break;
	}

	return pb;
}

// ===============================================================

// =====================lennard_jones.c===========================
double lennard_jones(struct position x)
{
	/*
		This is for black-box optimisation. Therefore, we are not supposed to know
		that there are some symmetries. That is why the dimension of the problem is
		3*nb_of_atoms, as it could be 3*nb_of_atoms-6
	*/
	int dim = 3;
	double dist;
	double f;
	int nPoints = x.size / dim;
	int alpha = 6;
	struct position x1, x2;
	double zz;

	x1.size = dim; x2.size = dim;

	f = 0;
	for (int i = 0; i < nPoints - 1; i++)
	{
		for (int d = 0; d < dim; d++)
			x1.x[d] = x.x[3 * i + d];
		for (int j = i+1; j < nPoints; j++)
		{
			for (int d = 0; d < dim; d++)
				x2.x[d] = x.x[3 * j + d];
			dist = distanceL(x1, x2, 2);
			if (dist < zero)
				zz = infinity;
			else
				zz = pow(dist, -alpha);
			f = f + zz * (zz - 1);
		}
	}
	f = 4 * f;
	return f;
}
// ===============================================================

// ==========================alea.c===============================

double alea(double a, double b, int option)
{
	/* random number (uniform distribution) in [a b]*/
	int close;
	double r;
	int i;
	int ir;

	double pw;
	int scanNb;
	float zz;

	switch (option)
	{
	case 0:						// KISS
		r = (double)rand_kiss() / RAND_MAX_KISS;
		break;
	case 10:					// Mersenne
		r = genrand64_real2();
		break;
	case -2:					// Truncated KISS(simulated)
		r = (double)alea_integer(0, rMax, 0) / rMax;
	case -1:
		r = (double)rand() / RAND_MAX;
		break;
	default:
		break;
	}
	return r;
}

// ---------------------------------------------------------------
int alea_integer(int a, int b, int option)
{
	// Integer random number in [a b]
	int ir;
	double r;

	//r = alea (0, 1);
	//ir = (int) (a + r * (b + 1 - a));
	//if (ir > b)	ir = b;

	if(a == b) return a;
	r = alea((double)a, (double)b + 1, option);
	ir = floor(r);

	if (ir > b) ir = b;
	return ir;
}

// ---------------------------------------------------------------
void aleaIndex(int index[], int S, int option)
{
	int indexTemp[S_max];
	int length;
	int rank;
	int s;
	int t;

	length = S;
	for (s = 0; s < S; s++) indexTemp[s] = s; //=index[s];

	for (s = 0; s < S; s++)
	{
		rank = alea_integer(0, length - 1, option);
		index[s] = indexTemp[rank];
		if (rank < length - 1)	// Compact
		{
			for (t = rank; t < length - 1; t++)
				indexTemp[t] = indexTemp[t + 1];
		}
		length = length - 1;
	}
}

// -----------------------------------------------------------
double alea_stable(double alpha, double beta, double nu, double delta, int option)
{
	// CMS algorithm (J.M. Chambers, C.L. Mallows and B.W. Stuck)
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
	double min = zero; //0; // To avoid to have to compute ln(0) ...
	double max = 1;
	double phi0, phi;
	double r;
	double s;
	double t1, t2, t3;
	double tau;
	double u;
	double temp;
	double w;
	double z;

	if (alpha < 0 || alpha>2)
	{
		printf("\n alpha %f ", alpha);
		ERROR("alea_levy. alpha must be in ]0,2]");
	}
	if (nu < 0)
	{
		printf("\n nu %f ", nu);
		ERROR("alea_levy. nu must be positive");
	}
	//--------------------------------------------
	if (alpha < 1) kAlpha = alpha; else kAlpha = 2 - alpha;

	phi0 = 0.5 * beta * (kAlpha / alpha);

	temp = tan(alpha * phi0);

	if (fabs(alpha - 1) < zero) betaPrime = beta;
	else
		betaPrime = -tan(0.5 * pi * (1 - alpha)) * temp;

	u = alea(min, max, option);
	phi = pi * (u - 0.5);

	eps = 1 - alpha;
	tau = -eps * temp;

	t1 = tan(0.5 * phi);
	t2 = tan(0.5 * eps * phi);
	t3 = 2 * t2 / (eps * phi);

	w = -log(alea(min, max, option));

	z = cos(eps * phi) - tan(alpha * phi0) * sin(eps * phi) / (w * cos(phi));

	temp = pow(z, eps / alpha);
	d = temp / eps;

	s = tan(alpha * phi0) + temp * (sin(alpha * phi) - tan(alpha * phi0) * cos(alpha * phi)) / cos(phi);
	r = s * nu + delta;
	return r;
}

// -----------------------------------------------------------
double alea_normal(double mean, double std_dev, int option)
{
	// Use the polar form of the Box-Muller transformation to obtain a pseudo
	// random number from a Gaussian distribution 

	double x1, x2, w, y1;
	// WARNING. This method is valid only if alea (...) defines a (more or less)
	// uniform distribution. In particular, if "option" is so that the numbers are
	// read on a list, whose distribution is, say, Gaussian, then w will never be
	// smaller than 1 

	int count = 0;
	do
	{
		x1 = 2.0 * alea(0, 1, option) - 1.0;
		x2 = 2.0 * alea(0, 1, option) - 1.0;
		w = x1 * x1 + x2 * x2;
		count++;
		if (count > 10000) //RAND_MAX)
		{
			if (option == 4) printf("\n Random numbers read on a list");
			ERROR("alea_normal. Probable non uniform distribution for alea(...)");
		}
	} while (w >= 1.0); // WARNING. Infinite loop possible.  w may be _never_ < 1 !!


	w = sqrt(-2.0 * log(w) / w);
	y1 = x1 * w;

	if (alea(0, 1, option) < 0.5) y1 = -y1;
	y1 = y1 * std_dev + mean;
	return y1;
}


// -----------------------------------------------------------
struct vector alea_sphere(int D, double radius, int distrib, double mean,
	double sigma, int option)
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

	x.size = D;
	pw = 1. / (double)D;

	// ----------------------------------- Step 1.  Direction
	length = 0;
	for (j = 0; j < D; j++)
	{
		// Here a Gaussian distribution is needed
		//if(distrib<2 && option!=4) x.v[j]=alea_normal(0,1,option); 
		if (option != 4)
			x.v[j] = alea_normal(0, 1, option);// Gaussian (Box-Muller method)
		else // When the "random" numbers are read on a file, the BM method 
			  // may loop infinitly
			x.v[j] = alea_stable(2, 0, 1, 0, option); // Gaussian (CMS method)
		length = length + x.v[j] * x.v[j];
	}

	length = sqrt(length);
	//----------------------------------- Step 2. Random radius

	switch (distrib)
	{
	default: // Uniform distribution
		// Note that the distribution in the hypersphere is NOT uniform
	case 0:
		r = alea(0, 1, option);
		break;

	case -1: // So that the final distribution is uniform in the sphere
		r = alea(0, 1, option);
		r = pow(r, pw);
		break;

	case 1: //Gaussian (Box-Muller)
		r = fabs(alea_normal(mean, sigma, option));
		break;

	case 2: //Gaussian (CMS)
		r = fabs(alea_stable(2, 0, sigma, mean, option));
		break;

	case 3: //TEST
		//r=fabs(alea_stable(0.2,0,sigma,alea(0,1))); 
		r = fabs(alea_stable(0.2, 0, sigma, 0.5 * (alea(0, 1, option) + alea(0, 1, option)), option));
		break;

	case 4: // Slash distribution
		r = fabs(alea_normal(mean, sigma, option)) / fabs(alea_normal(0, 1, option));
		break;

	case 99: // Constant radius, for some specific uses (random ON the sphere)
		r = 1;
		break;
	}


	for (j = 0; j < D; j++)
	{
		x.v[j] = radius * r * x.v[j] / length;
	}
	return x;
}


//==========================================================================    
//struct vectorList quasiRand(int D, int nRand, int option)
//{
//	/*
//	 Generate nRand vectors of size D that are the coordinates of
//	 nRand quasi-random points in [0,1]^D
//
//		 For C under Linux, you need to use the GSL library:
//#include <gsl/gsl_qrng.h> // Do not forget to link to gsl and gslcblast
//
//		 */
//
//	int i;
//	gsl_qrng* qrng_q;
//	struct vectorList qRand;
//
//	switch (option)
//	{
//	case 1: // Sobol
//		if (D > 40)
//		{
//			printf("\nThe embedded Sobol sequences generator can not be used for dimensions greater than 40");
//			printf("\nYou should use Halton sequences (for dimensions up to 1229)");
//			ERROR("\n I stop here");
//
//		}
//		qrng_q = gsl_qrng_alloc(sobol, D);
//		break;
//
//	default: // Halton
//		if (D > 1229)
//		{
//			printf("\nThe embedded Halton sequences generator can not be used for dimensions greater than 1229");
//			printf("\nSorry");
//			ERROR("\n I stop here");
//		}
//		qrng_q = gsl_qrng_alloc(halton, D);
//		break;
//	}
//
//	gsl_qrng* q = qrng_q;
//
//	for (i = 0; i < nRand; i++)
//	{
//		gsl_qrng_get(q, qRand.V[i].v);
//	}
//
//	gsl_qrng_free(q);
//
//	return qRand;
//}


// ===============================================================

// ==========================PSO.c================================
struct result* PSO(struct param* param, struct problem* pb)
{
	int d;
	double error = 0;
	struct result* R;
	R = (struct result*)malloc(sizeof(struct result));
	double errorPrev;
	int g;
	struct position* Gr;
	Gr = (struct position*)malloc(sizeof(struct position));
	struct vector* GX;
	GX = (struct vector*)malloc(sizeof(struct vector));
	int index[S_max];
	int initLinks;	// Flag to (re)init or not the information links
	int iter; 		// Iteration number (time step)
	int iterBegin;
	int LINKS[S_max][S_max];	// Information links
	int m;
	int noStop;
	double out;
	double p;
	struct vector* PX;
	PX = (struct vector*)malloc(sizeof(struct vector));
	struct vectorList* qRand;
	qRand = (struct vectorList*)malloc(sizeof(struct vectorList));
	double rad;
	int randCase;
	int s0, s, s1;
	struct vector *V1, *V2;
	V1 = (struct vector*)malloc(sizeof(struct vector));
	V2 = (struct vector*)malloc(sizeof(struct vector));
	double w1, w2, w3;
	double xMax, xMin;
	struct position* XPrev;
	XPrev = (struct position*)malloc(sizeof(struct position));
	double zz;

	if (PX == NULL) return NULL;
	PX->size = pb->SS.D;
	if (GX == NULL) return NULL;
	GX->size = pb->SS.D;
	if (V1 == NULL) return NULL;
	V1->size = pb->SS.D;
	if (V2 == NULL) return NULL;
	V2->size = pb->SS.D;
	if (Gr == NULL) return NULL;
	Gr->size = pb->SS.D;

	// -----------------------------------------------------
	// INITIALISATION
	p = param->p; // Probability threshold for random topology
	R->SW.S = param->S; // Size of the current swarm
	randCase = param->BW[2];
	if (randCase > 300)  randCase = 3;
	if (randCase < -200) randCase = -2;

	// Position and velocity
	for (s = 0; s < R->SW.S; s++)
	{
		R->SW.X[s].size = pb->SS.D;
		R->SW.V[s].size = pb->SS.D;
	}

	if (pb->SS.normalise > 0) { xMin = 0; xMax = pb->SS.normalise; } // [0,normalise]^D

	switch (randCase)
	{

	default:
		for (s = 0; s < R->SW.S; s++)
		{

			for (d = 0; d < pb->SS.D; d++)
			{
				if (pb->SS.normalise == 0) { xMin = pb->SS.min[d]; xMax = pb->SS.max[d]; }
				R->SW.X[s].x[d] = alea(xMin, xMax, randCase);
				R->SW.V[s].v[d] = alea(xMin - R->SW.X[s].x[d], xMax - R->SW.X[s].x[d], randCase);
				// So that  xMin<= x+v <= xMax
			}
		}
		break;
	case 1:
	case 2:
		/**qRand = quasiRand(pb->SS.D, R->SW.S, randCase);

		for (s = 0; s < R->SW.S; s++)
		{
			for (d = 0; d < pb->SS.D; d++)
			{
				if (pb->SS.normalise == 0) { xMin = pb->SS.min[d]; xMax = pb->SS.max[d]; }
				R->SW.X[s].x[d] = xMin + (xMax - xMin) * qRand->V[s].v[d];
				R->SW.V[s].v[d] = alea(xMin - R->SW.X[s].x[d], xMax - R->SW.X[s].x[d], randCase);
			}
		}*/
		break;
	}

	// Take quantisation into account
	if (pb->SS.quantisation == 1)
	{
		for (s = 0; s < R->SW.S; s++)
		{
			R->SW.X[s] = quantis(R->SW.X[s], pb->SS);
		}
	}

	// First evaluations
	errMax = 0;
	errMin = infinity;
	for (s = 0; s < R->SW.S; s++)
	{
		R->SW.X[s].f = perf(R->SW.X[s], pb->function, pb->SS, pb->objective);
		R->SW.P[s] = R->SW.X[s];	// Best position = current one
	}

	// If the number max of evaluations is smaller than 
	// the swarm size, just keep evalMax particles, and finish
	if (R->SW.S > pb->evalMax) R->SW.S = pb->evalMax;
	R->nEval = R->SW.S;

	// Find the best
	R->SW.best = 0;

	errorPrev = R->SW.P[R->SW.best].f; // "distance" to the wanted f value (objective)

	for (s = 1; s < R->SW.S; s++)
	{
		zz = R->SW.P[s].f;
		if (zz < errorPrev)
		{
			R->SW.best = s;
			errorPrev = zz;
		}
	}

	// Display the best
	//	printf( "\nBest value after init. %1.20e ", errorPrev ); printf(" ");
	//		printf( "\n Position :\n" );
	//		for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );

	// ---------------------------------------------- ITERATIONS	
	noStop = 0;
	error = errorPrev;
	iter = 0; iterBegin = 0;
	initLinks = 1;		// So that information links will be initialized

	// Each particle informs itself
	for (m = 0; m < R->SW.S; m++) LINKS[m][m] = 1;

	while (noStop == 0)
	{
		iter = iter + 1;
		/*
		// Display the swarm
			printf("\n Positions (%i) \ Velocities (%i) after iter %i.",pb.SS.D,pb.SS.D, iter-1 );
			for (s = 0; s < R.SW.S; s++)
			{
			printf("\n");
			for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.X[s].x[d] );
			printf("\ ");
			for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.V[s].v[d] );
			}
		*/
		switch (param->BW[3])
		{
		case 0:
			for (s = 0; s < R->SW.S; s++)  index[s] = s;  // No random permutation
			break;
		case 1:
			aleaIndex(index, R->SW.S, randCase); // Random numbering of the particles
			break;
		default:
			break;
		}

		if (initLinks == 1)	// Modify topology
		{
			switch (param->topology)
			{
			case 1: // Random ring
				// Init to zero (no link)
				for (s = 0; s < R->SW.S; s++)
				{
					for (m = 0; m < R->SW.S; m++)    LINKS[m][s] = 0;
				}

				// Information links (bidirectional ring)
				for (s = 0; s < R->SW.S - 1; s++)
				{
					LINKS[index[s]][index[s + 1]] = 1;
					LINKS[index[s + 1]][index[s]] = 1;
				}

				LINKS[index[0]][index[R->SW.S - 1]] = 1;
				LINKS[index[R->SW.S - 1]][index[0]] = 1;

				// Each particle informs itself
				for (m = 0; m < R->SW.S; m++) LINKS[m][m] = 1;
				break;
			default: // case 0. As in SPSO 2007  Who informs who, at random
				for (s = 0; s < R->SW.S; s++)
				{
					for (m = 0; m < R->SW.S; m++)
					{
						if (m == s) continue;

						if (alea(0, 1, randCase) < p)
						{
							LINKS[m][s] = 1;	// Probabilistic method
						}
						else LINKS[m][s] = 0;
					}
				}
				break;
			}
		}
		/*
			// Display the links
			for(s=0;s<R.SW.S;s++)
			{
				printf("\n");
				for(m=0;m<R.SW.S;m++)
				{
				printf("%i ",LINKS[s][m]);
				}
			}
		*/

		// Loop on particles
		for (s0 = 0; s0 < R->SW.S; s0++) // For each particle ...
		{
			s = index[s0];
			// ... find the first informant
			s1 = 0;
			while (LINKS[s1][s] == 0)	s1++;
			if (s1 >= R->SW.S)	s1 = s;

			// Find the best informant			
			g = s1;
			for (m = s1; m < R->SW.S; m++)
			{
				if (LINKS[m][s] == 1 && R->SW.P[m].f < R->SW.P[g].f)
					g = m;
			}

			//.. compute the new velocity, and move

			// Exploration tendency
			for (d = 0; d < pb->SS.D; d++)
			{
				R->SW.V[s].v[d] = param->w * R->SW.V[s].v[d];
			}

			// Prepare Exploitation tendency  p-x
			for (d = 0; d < pb->SS.D; d++)
			{
				PX->v[d] = R->SW.P[s].x[d] - R->SW.X[s].x[d];
			}

			if (g != s) // If the particle is not its own local best, prepare g-x
			{
				for (d = 0; d < pb->SS.D; d++)
				{
					GX->v[d] = R->SW.P[g].x[d] - R->SW.X[s].x[d];
				}
			}
			else // this is the best particle. Define another random "previous best" May be used or not, though (see below)
			{

				if(param->BW[1]==1) // WARNING: un-comment this line largely modify the performances
									  // of list-based option BW[2]=4. Not the same lists are valid
				{
				s1:
					s1 = alea_integer(0, param->S - 1, randCase);

					if (s1 == s) goto s1;
					// *** WARNING, may be infinite
				}
			}
			// Gravity centre Gr
			w1 = 1; w2 = 1; w3 = 1;

			switch (param->BW[1])
			{
			default:if (g == s) w3 = 0; break; // Pure standard
			case 1: break; // Will use a specific method (see below)
			case 2: break;
			case 3: if (g == s) { w1 = 0; w3 = 0; } break;
			case 4: // TEST
				w1 = R->SW.P[g].f; w2 = R->SW.P[s].f;
				if (g != s)
				{
					w3 = R->SW.X[s].f;
				}
				else
				{
					w3 = 0;
				}
				/*
				w[0]=alea(0,1);w[1]=alea(0,1);
				if(g!=s)
				{
					w[2]=alea(0,1);
					qsort(w,3,sizeof(w[0]),compareDoubles);
				w1=w[0]; w2=w[1]; w3=w[2];

				}
				else
				{
					w3=0;
					w1=min(w[0],w[1]);
					w2=max(w[0],w[1]);
				}
				*/
				break;
			}
			zz = w1 + w2 + w3;
			w1 = w1 / zz; w2 = w2 / zz; w3 = w3 / zz;

			for (d = 0; d < pb->SS.D; d++)
			{
				Gr->x[d] = w1 * R->SW.X[s].x[d] + w2 * (R->SW.X[s].x[d] + param->c * PX->v[d]);

				if (g != s)
				{
					Gr->x[d] = Gr->x[d] + w3 * (R->SW.X[s].x[d] + param->c * GX->v[d]);
				}
				else // Here P=G
					switch (param->BW[1])
					{
					default: break; //0. If "pure standard", do nothing.
					case 1: //G is a random informant (see above). 
						// However, it is used only for high dimension
						// and not even always (probabilistic choice)
						if (alea(0, 1, randCase) < (double)param->S / (pb->SS.D * param->K))	// Approximation
						{ // Use P twice
							Gr->x[d] = Gr->x[d] + w3 * (R->SW.X[s].x[d] + param->c * PX->v[d]);
						}
						else
						{	// Random informant
							Gr->x[d] = Gr->x[d] + w3 * (R->SW.X[s].x[d]
								+ param->c * (R->SW.P[s1].x[d] - R->SW.X[s].x[d]));
						}
						break;
					case 2: // More conservative
						Gr->x[d] = R->SW.P[g].x[d];
						break;
					}
				V1->v[d] = Gr->x[d] - R->SW.X[s].x[d]; // Vector X-G
			}

			// Random point around
			switch (param->BW[1])
			{
			default:
				rad = distanceL(R->SW.X[s], *Gr, 2);
				break;
			case 3:
				rad = (param->c - 1) * distanceL(R->SW.X[s], R->SW.P[s], 2);
				break;
			}

			*V2 = alea_sphere(pb->SS.D, rad, param->distrib, param->mean, param->sigma, randCase);

			// New "velocity"
			for (d = 0; d < pb->SS.D; d++)
			{
				R->SW.V[s].v[d] = R->SW.V[s].v[d] + V1->v[d] + V2->v[d]; // New "velocity"
			}

			// New position

			for (d = 0; d < pb->SS.D; d++)
			{
				R->SW.X[s].x[d] = R->SW.X[s].x[d] + R->SW.V[s].v[d];
			}

			if (R->nEval >= pb->evalMax)  goto end;

			// --------------------------

			if (pb->SS.quantisation == 1)
				R->SW.X[s] = quantis(R->SW.X[s], pb->SS);

			// Confinement			
			out = 0;
			switch (param->confin)
			{
			default:
				for (d = 0; d < pb->SS.D; d++)
				{
					if (pb->SS.normalise == 0) { xMin = pb->SS.min[d]; xMax = pb->SS.max[d]; }

					if (R->SW.X[s].x[d] < xMin)
					{
						R->SW.X[s].x[d] = xMin;
						R->SW.V[s].v[d] = -0.5 * R->SW.V[s].v[d];
						out = 1;
					}
					else
					{
						if (R->SW.X[s].x[d] > xMax)
						{
							R->SW.X[s].x[d] = xMax;
							R->SW.V[s].v[d] = -0.5 * R->SW.V[s].v[d];
							out = 1;
						}
					}
				}
				break;
			case 1: // No confinement and no evaluation if outside ("let if fly")
				for (d = 0; d < pb->SS.D; d++)
				{
					if (pb->SS.normalise == 0) { xMin = pb->SS.min[d]; xMax = pb->SS.max[d]; }

					if (R->SW.X[s].x[d] < xMin || R->SW.X[s].x[d] > xMax) out = 1;
				}
				break;
			}

			if (pb->SS.quantisation == 1 && out > 0)
			{
				R->SW.X[s] = quantis(R->SW.X[s], pb->SS);
			}
			// If the position is inside
			if (param->confin == 0 || (param->confin == 1 && out < zero))
			{
				// Evaluation
				R->SW.X[s].f = perf(R->SW.X[s], pb->function, pb->SS, pb->objective);
				R->nEval = R->nEval + 1;
				// ... update the best previous position		
				if (R->SW.X[s].f < R->SW.P[s].f)	// Improvement
				{
					R->SW.P[s] = R->SW.X[s];

					// ... update the best of the bests
					if (R->SW.P[s].f < R->SW.P[R->SW.best].f)
					{
						R->SW.best = s;
					}
				}
			}
			if (param->trace > 0)
			{
				// Keep trace of every position, for further tests
				fprintf(f_trace, "%i %f ", s, R->SW.X[s].f);
				for (d = 0; d < pb->SS.D; d++)
				{
					fprintf(f_trace, "%f ", R->SW.X[s].x[d]);
				}
				fprintf(f_trace, "\n");
			}
		}			// End of "for (s0=0 ...  "	

		// Check if finished
		error = R->SW.P[R->SW.best].f;

		if (error < errorPrev)	// Improvement of the global best
		{
			initLinks = 0;
		}
		else			// No global improvement
		{
			initLinks = 1;	// Information links will be	reinitialized	
		}

		errorPrev = error;
	end:

		if (error > pb->epsilon && R->nEval < pb->evalMax)
		{
			noStop = 0;	// Won't stop
		}
		else
		{
			noStop = 1;	// Will stop
		}
	} // End of "while nostop ...

	// printf( "\n and the winner is ... %i", R.SW.best );			
	R->error = error;
	return R;
}

struct position quantis(struct position x, struct SS SS)
{
	/*
	 Quantisation of a position
	 Only values like x+k*q (k integer) are admissible
		 */
	int d;
	double qd;
	struct position quantx;

	quantx = x;
	for (d = 0; d < x.size; d++)
	{
		qd = SS.q.q[d];

		if (qd > zero)	// Note that qd can't be < 0
		{
			quantx.x[d] = qd * floor(0.5 + x.x[d] / qd);
		}
	}
	return quantx;
};
// ===============================================================

// ==========================Perf.c===============================
double perf(struct position x, int function, struct SS SS, double objective)
{
	// Evaluate the fitness value for the particle of rank s 
	double beta;
	double c;
	int d;
	double DD;
	double dx1, dx2;
	int grid;
	int i, j;
	int  k;
	double min, max;
	int n;
	struct fitness ff = { 0 };
	double f = 0, p, xd, x1, x2, x3, x4, x5, x6;
	double s11, s12, s21, s22;
	double sum1, sum2;
	double t0, tt, t1;
	double theta;
	double u;
	struct position xs;
	double y1, y2;
	double z;

	// ---------------------cec2005data.c------------------
	// Shifted Parabola/Sphere (CEC 2005 benchmark)		
	static double offset_0[30] =
	{
	-3.9311900e+001, 5.8899900e+001, -4.6322400e+001, -7.4651500e+001, -1.6799700e+001,
	-8.0544100e+001, -1.0593500e+001, 2.4969400e+001, 8.9838400e+001, 9.1119000e+000,
	-1.0744300e+001, -2.7855800e+001, -1.2580600e+001, 7.5930000e+000, 7.4812700e+001,
	 6.8495900e+001, -5.3429300e+001, 7.8854400e+001, -6.8595700e+001, 6.3743200e+001,
	 3.1347000e+001, -3.7501600e+001, 3.3892900e+001, -8.8804500e+001, -7.8771900e+001,
	-6.6494400e+001, 4.4197200e+001, 1.8383600e+001, 2.6521200e+001, 8.4472300e+001
	};
	// Shifted Rosenbrock (CEC 2005 benchmark)
	static double offset_2[30] =
	{
	8.1023200e+001, -4.8395000e+001,  1.9231600e+001, -2.5231000e+000,  7.0433800e+001,
	 4.7177400e+001, -7.8358000e+000, -8.6669300e+001,  5.7853200e+001, -9.9533000e+000,
	  2.0777800e+001,  5.2548600e+001,  7.5926300e+001,  4.2877300e+001, -5.8272000e+001,
	 -1.6972800e+001,  7.8384500e+001,  7.5042700e+001, -1.6151300e+001,  7.0856900e+001,
	 -7.9579500e+001, -2.6483700e+001,  5.6369900e+001, -8.8224900e+001, -6.4999600e+001,
	 -5.3502200e+001, -5.4230000e+001,  1.8682600e+001, -4.1006100e+001, -5.4213400e+001
	};
	// Shifted Rastrigin (CEC 2005)
	static double offset_3[30] =
	{
	1.9005000e+000, -1.5644000e+000, -9.7880000e-001, -2.2536000e+000,  2.4990000e+000,
	-3.2853000e+000,  9.7590000e-001, -3.6661000e+000,  9.8500000e-002, -3.2465000e+000,
	3.8060000e+000, -2.6834000e+000, -1.3701000e+000,  4.1821000e+000,  2.4856000e+000,
	-4.2237000e+000,  3.3653000e+000,  2.1532000e+000, -3.0929000e+000,  4.3105000e+000,
	-2.9861000e+000,  3.4936000e+000, -2.7289000e+000, -4.1266000e+000, -2.5900000e+000,
	 1.3124000e+000, -1.7990000e+000, -1.1890000e+000, -1.0530000e-001, -3.1074000e+000
	};

	// Shifted Schwefel (F2 CEC 2005. Also for F4)
	static double offset_4[30] =
	{
	  3.5626700e+001, -8.2912300e+001, -1.0642300e+001, -8.3581500e+001,  8.3155200e+001,
	  4.7048000e+001, -8.9435900e+001, -2.7421900e+001,  7.6144800e+001, -3.9059500e+001,
	  4.8885700e+001, -3.9828000e+000, -7.1924300e+001,  6.4194700e+001, -4.7733800e+001,
	 -5.9896000e+000 ,-2.6282800e+001, -5.9181100e+001,  1.4602800e+001, -8.5478000e+001,
	 -5.0490100e+001,  9.2400000e-001,  3.2397800e+001,  3.0238800e+001, -8.5094900e+001,
	  6.0119700e+001, -3.6218300e+001, -8.5883000e+000, -5.1971000e+000,  8.1553100e+001
	};

	// Shifted Griewank (CEC 2005)
	static double offset_5[30] =
	{
	 -2.7626840e+002, -1.1911000e+001, -5.7878840e+002, -2.8764860e+002, -8.4385800e+001,
	 -2.2867530e+002, -4.5815160e+002, -2.0221450e+002, -1.0586420e+002, -9.6489800e+001,
	 -3.9574680e+002, -5.7294980e+002, -2.7036410e+002, -5.6685430e+002, -1.5242040e+002,
	 -5.8838190e+002, -2.8288920e+002, -4.8888650e+002, -3.4698170e+002, -4.5304470e+002,
	 -5.0658570e+002, -4.7599870e+002, -3.6204920e+002, -2.3323670e+002, -4.9198640e+002,
	 -5.4408980e+002, -7.3445600e+001, -5.2690110e+002, -5.0225610e+002, -5.3723530e+002
	};

	// Shifted Ackley (CEC 2005)
	static double offset_6[30] =
	{
	 -1.6823000e+001,  1.4976900e+001,  6.1690000e+000,  9.5566000e+000,  1.9541700e+001,
	 -1.7190000e+001, -1.8824800e+001,  8.5110000e-001, -1.5116200e+001,  1.0793400e+001,
	  7.4091000e+000,  8.6171000e+000, -1.6564100e+001, -6.6800000e+000,  1.4543300e+001,
	  7.0454000e+000, -1.8621500e+001,  1.4556100e+001, -1.1594200e+001, -1.9153100e+001,
	 -4.7372000e+000,  9.2590000e-001,  1.3241200e+001, -5.2947000e+000,  1.8416000e+000,
	  4.5618000e+000, -1.8890500e+001,  9.8008000e+000, -1.5426500e+001,  1.2722000e+000
	};
	// ----------------------------------------------------


	static float bts[19][2] =
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

	float btsPenalty = 100;

	double z1, z2;

	if (SS.normalise > 0)
	{
		// Back to the real search space
		xs.size = x.size;
		for (d = 0; d < xs.size; d++)
			xs.x[d] = SS.min[d] + (SS.max[d] - SS.min[d]) * x.x[d] / SS.normalise;
	}
	else xs = x;

	switch (function)
	{

	case -1:		// Constant. For test of biases
		f = 0;
		break;
	case 0:		// Parabola (Sphere)
		f = 0;

		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d] - d;
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
			p = p * cos(xd / sqrt((double)(d + 1)));
		}
		f = f / 4000 - p + 1;
		break;
	case 2:		// Rosenbrock
		f = 0;
		t0 = xs.x[0] + 1;	// Solution on (0,...0) when
		// offset=0
		for (d = 1; d < xs.size; d++)
		{

			t1 = xs.x[d] + 1;
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
			f = f + xd * xd - k * cos(2 * pi * xd);
		}
		f = f + xs.size * k;
		break;
	case 4:		// 2D Tripod function
		// Note that there is a big discontinuity right on the solution
		// point. 
		x1 = xs.x[0];
		x2 = xs.x[1];
		s11 = (1.0 - sign(x1)) / 2;
		s12 = (1.0 + sign(x1)) / 2;
		s21 = (1.0 - sign(x2)) / 2;
		s22 = (1.0 + sign(x2)) / 2;

		//f = s21 * (fabs (x1) - x2); // Solution on (0,0)
		f = s21 * (fabs(x1) + fabs(x2 + 50)); // Solution on (0,-50)  
		f = f + s22 * (s11 * (1 + fabs(x1 + 50) +
			fabs(x2 - 50)) + s12 * (2 +
				fabs(x1 - 50) +
				fabs(x2 - 50)));

		//f=log(1+f);	
		break;
	case 5:  // Ackley
		sum1 = 0;
		sum2 = 0;
		DD = x.size;
		pi = acos(-1);
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			sum1 = sum1 + xd * xd;
			sum2 = sum2 + cos(2 * pi * xd);
		}
		f = -20 * exp(-0.2 * sqrt(sum1 / DD)) - exp(sum2 / DD) + 20 + exp(1);

		break;

	case 6: // Schwefel
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f - xd * sin(sqrt(fabs(xd)));
		}
		break;

	case 7: // Schwefel 1.2
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			sum1 = 0;
			for (k = 0; k <= d; k++) sum1 = sum1 + xd;
			f = f + sum1 * sum1;
		}
		break;

	case 8: // Schwefel 2.22
		sum1 = 0; sum2 = 1;
		for (d = 0; d < x.size; d++)
		{
			xd = fabs(xs.x[d]);
			sum1 = sum1 + xd;
			sum2 = sum2 * xd;
		}
		f = sum1 + sum2;
		break;

	case 9: // Neumaier 3
		sum1 = 0; sum2 = 1;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d] - 1;
			sum1 = sum1 + xd * xd;
		}
		for (d = 1; d < x.size; d++)
		{
			sum2 = sum2 + xs.x[d] * xs.x[d - 1];
		}

		f = sum1 + sum2;
		break;

	case 10: // G3 (constrained) 
		// min =0 on (1/sqrt(D), ...)
		f = 1;
		sum1 = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f * xd;
			sum1 = sum1 + xd * xd;
		}
		f = fabs(1 - pow(x.size, x.size / 2) * f) + x.size * fabs(sum1 - 1);
		break;

	case 11: // Network  btsNb BTS, bcdNb BSC

		f = 0;
		// Constraint: each BTS has one link to one BSC 
		for (d = 0; d < btsNb; d++)
		{
			sum1 = 0;
			for (k = 0; k < bcsNb; k++) sum1 = sum1 + xs.x[d + k * btsNb];
			if (sum1 < 1 - zero || sum1>1 + zero) f = f + btsPenalty;

		}
		// Distances
		for (d = 0; d < bcsNb; d++) //For each BCS d
		{
			for (k = 0; k < btsNb; k++) // For each BTS k
			{
				if (xs.x[k + d * btsNb] < 1) continue;
				// There is a link between BTS k and BCS d
				n = bcsNb * btsNb + 2 * d;
				z1 = bts[k][0] - xs.x[n];
				z2 = bts[k][1] - xs.x[n + 1];
				f = f + sqrt(z1 * z1 + z2 * z2);
			}
		}
		break;

	case 12: // Schwefel
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f - xd * sin(sqrt(fabs(xd)));
		}
		break;

	case 13: // 2D Goldstein-Price function
		x1 = xs.x[0]; x2 = xs.x[1];

		f = (1 + pow(x1 + x2 + 1, 2) * (19 - 14 * x1 + 3 * x1 * x1 - 14 * x2 + 6 * x1 * x2 + 3 * x2 * x2))
			* (30 + pow(2 * x1 - 3 * x2, 2) *
				(18 - 32 * x1 + 12 * x1 * x1 + 48 * x2 - 36 * x1 * x2 + 27 * x2 * x2));
		break;

	case 14:  //Schaffer F6
		x1 = xs.x[0]; x2 = xs.x[1];
		f = 0.5 + (pow(sin(sqrt(x1 * x1 + x2 * x2)), 2) - 0.5) / pow(1.0 + 0.001 * (x1 * x1 + x2 * x2), 2);

		break;

	case 15: // Step
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = (int)(xs.x[d] + 0.5);
			f = f + xd * xd;
		}
		break;

	case 16: // Schwefel 2.21
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = fabs(xs.x[d]);
			if (xd > f) f = xd;
		}
		break;

	case 17: // Lennard-Jones
		f = lennard_jones(xs);
		break;

	case 18: // Gear train
		f = pow(1. / 6.931 - x.x[0] * x.x[1] / (x.x[2] * x.x[3]), 2);
		//	f=pow(fabs(1./6.0 -x.x[0]*x.x[1]/(x.x[2]*x.x[3])),2);

		break;

	case 19: // Sine-sine function
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f - sin(xd) * pow(sin((d + 1) * xd * xd / pi), 20);

		}
		break;

	case 20: // Perm function
		beta = 10;
		f = 0;
		for (k = 0; k < x.size; k++)
		{
			sum1 = 0;
			for (d = 0; d < x.size; d++)
			{
				xd = xs.x[d];
				sum1 = sum1 + (pow(d + 1, k) + beta) * (pow(xd / (d + 1), k) - 1);
			}
			sum1 = sum1 * sum1;
			f = f + sum1;
		}

		break;

	case 21: // Coil compression spring  (penalty method)
		// Ref New Optim. Tech. in Eng. p 644

		x1 = xs.x[0]; // {1,2, ... 70}
		x2 = xs.x[1];//[0.6, 3]
		x3 = xs.x[2];// relaxed form [0.207,0.5]  dx=0.001
		// In the original problem, it is a list of
		// acceptable values
		// {0.207,0.225,0.244,0.263,0.283,0.307,0.331,0.362,0.394,0.4375,0.5}

		f = pi * pi * x2 * x3 * x3 * (x1 + 2) * 0.25;
		//	f=x2*x3*x3*(x1+2);
		// Constraints
		ff = constraint(xs, function);

		if (ff.f[1] > 0) { c = 1 + ff.f[1]; f = f * c * c * c; }
		if (ff.f[2] > 0) { c = 1 + ff.f[2]; f = f * c * c * c; }
		if (ff.f[3] > 0) { c = 1 + ff.f[3]; f = f * c * c * c; }
		if (ff.f[4] > 0) { c = 1 + pow(10, 10) * ff.f[4]; f = f * c * c * c; }
		if (ff.f[5] > 0) { c = 1 + pow(10, 10) * ff.f[5]; f = f * c * c * c; }
		break;

		// -------------------cellular_phone.c----------------
	case 22:
		// Cellular phone
				// Grid 100 x 100. You may modify it for more or less precision
				// (which implies, of course, more or less computation time)
		grid = 100;
		min = 0; max = 1; // Warning. Hard coded. Must be the same as in problemDef
		//f=0;

		dx1 = (max - min) / grid;
		dx2 = (max - min) / grid;

		// For each point of the grid, compute the maximum field generated 
		// by the stations. The aim is to maximize the smallest maximum
		f = infinity;

		for (i = 0; i <= grid; i++)
		{
			y1 = min + i * dx1;

			for (j = 0; j <= grid; j++)
			{
				y2 = min + j * dx2;
				z = 0; // Max known field
				for (d = 0; d < xs.size - 1; d = d + 2) // Loop on station positions
				{
					x1 = xs.x[d]; // First coordinate
					x2 = xs.x[d + 1]; // Second coordinate

					//	z2=1./((x1-i)*(x1-y1) +(x2-j)*(x2-y2)+1);
					// Field generated by the station (d, d+1)
					z2 = 1. / ((x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2) + 0.0001 * dx1);
					// If higher than already known, keep it
					if (z2 > z) z = z2;
				}

				//	f=f+z;
					// At this point, the maximum field generated is z
			// If it is smaller than in previous checked points, keep it
				if (z < f) f = z;

			}

		}
		// We want it as high as possible
		f = 1. / f; // In order to have something to minimise
		break;
		// ---------------------------------------------------

	case 23: // Penalized
		f = pow(sin(pi * xs.x[0]), 2);
		for (d = 1; d < x.size - 1; d++)
		{
			f = f + pow(xs.x[d], 2) * (1 + pow(sin(3 * pi * xs.x[d + 1]), 2));
		}
		f = 0.1 * (f + pow(xs.x[x.size - 2], 2) * (1 + pow(sin(2 * pi * xs.x[x.size - 1]), 2)));

		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			if (xd > 5) { u = 100 * pow(xd - 5, 4); f = f + u; }
			if (xd < -5) { u = 100 * pow(-xd - 5, 4); f = f + u; }
		}

		break;

		// --------------------Repulsion----------------------
	case 24:
		// Repulsion
		min = 0; max = 1; // Warning. Hard coded. Musted be the same asa in problemDef;
		n = (int)xs.size / Dim; // Number of points;
		u = 0.5 * (n - 1); // For the repulsive action of the boundaries The smaller u, the nearer to the boundaries canbe some points
		f = 0; // Total repulsive force to minimise
		for (int i = 0; i < n; i++) // For each charged point ...
		{
			// ... add the repulsive force between the point and the boundaries
			for (int d = 0; d < Dim; d++)
				f = f + u * (1. / pow(xs.x[d + Dim * i] - min, 2) + 1. / pow(max - xs.x[d + Dim * i], 2));

			// ... add the repulsive force between the charged points
			for (j = 0; j < n; j++)
			{
				if (j == i) continue;

				// Distance^2 between point i and point j1
				z = 0;
				for (d = 0; d < Dim; d++)
					z = z + pow(xs.x[d + Dim * i] - xs.x[d + Dim * j], 2);

				// Repulsion
				f = f + 1. / z;
			}
		}
		break;
		// ---------------------------------------------------


		// ---------------------pressure_vessel.c-------------
	case 25:
		// Ref New Optim. Tech. in Eng. p 638
		// D = 4  		

		x1 = xs.x[0]; // [1.1,12.5] granularity 0.0625
		x2 = xs.x[1];// [0.6,12.5] granularity 0.0625
		x3 = xs.x[2]; // [0,240]
		x4 = xs.x[3];// [0,240]

		f = 0.6224 * x1 * x3 * x4 + 1.7781 * x2 * x3 * x3 + 3.1611 * x1 * x1 * x4 + 19.84 * x1 * x1 * x3;
		//		f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1161*x1*x1*x4 + 19.84*x1*x1*x3;
		//  	f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1661*x1*x1*x4 + 19.84*x1*x1*x3;

		ff = constraint(xs, function);

		if (ff.f[1] > 0) { c = 1 + pow(10, 10) * ff.f[1];  f = f * c * c; }
		if (ff.f[2] > 0) { c = 1 + ff.f[2]; f = f * c * c; }
		if (ff.f[3] > 0) { c = 1 + ff.f[3]; f = f * c * c; }

		break;
		// ---------------------------------------------------
	case 26: // Ellipsoidal
		f = 0;

		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d] - d - 1;
			f = f + xd * xd;
		}
		break;
	case 27:		// Quadric
		f = 0;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[0];
			if (d > 0)
				for (j = 1; j < d; j++) xd = xd + xs.x[j];

			f = f + xd * xd;
		}
		break;

	case 28:// Frequency modulation sound parameter identification
		theta = 2 * pi / 100;
		f = 0;
		x1 = xs.x[0]; x2 = xs.x[1]; x3 = xs.x[2]; x4 = xs.x[3]; x5 = xs.x[4]; x6 = xs.x[5];

		for (d = 1; d <= 100; d++)
		{
			z = x1 * sin(x2 * d * theta + x3 * sin(x4 * d * theta + x5 * sin(x6 * d * theta)))
				- sin(5 * d * theta + 1.5 * sin(4.8 * d * theta + 2 * sin(4.9 * d * theta)));
			f = f + z * z;
		}
		break;

	case 999: // for tests
		f = 0;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d];
			f = f + pow(xd, 4) - 16 * xd * xd + 5 * xd;
		}
		f = f / xs.size;
		break;

	default:
		break;
	}

	//------------------
	f = fabs(f - objective);
	if (f < errMin) errMin = f; // For information
	if (f > errMax) { if (f < infinity) errMax = f; else errMax = infinity; } // For information

	return f;
}

//==========================================================
struct fitness constraint(struct position x, int functCode)
{
	// ff[0] is defined in perf()
	// Variables specific to Coil compressing spring
	static double	Fmax = 1000.0;
	static double	Fp = 300;
	double Cf;
	double K;
	double sp;
	double lf;

	static double	S = 189000.0;
	static double	lmax = 14.0;
	static double	spm = 6.0;
	static double	sw = 1.25;
	static double	G = 11500000;
	struct fitness ff = { 0 };
	ff.size = 1; // Default value

	switch (functCode)
	{
	case 21: // Compression Spring
		Cf = 1 + 0.75 * x.x[2] / (x.x[1] - x.x[2]) + 0.615 * x.x[2] / x.x[1];
		K = 0.125 * G * pow(x.x[2], 4) / (x.x[0] * x.x[1] * x.x[1] * x.x[1]);
		sp = Fp / K;
		lf = Fmax / K + 1.05 * (x.x[0] + 2) * x.x[2];

		ff.f[1] = 8 * Cf * Fmax * x.x[1] / (pi * x.x[2] * x.x[2] * x.x[2]) - S;
		ff.f[2] = lf - lmax;
		ff.f[3] = sp - spm;
		ff.f[4] = sw - (Fmax - Fp) / K;
		break;

	case 25: // Pressure vessel
		ff.f[1] = 0.0193 * x.x[2] - x.x[0];
		ff.f[2] = 0.00954 * x.x[2] - x.x[1];
		ff.f[3] = 750 * 1728 - pi * x.x[2] * x.x[2] * (x.x[3] + (4.0 / 3) * x.x[2]);
		break;

	default:
		break;
	}

	return ff;
}

// ===============================================================

// ==========================Tools.c==============================

// ---------------------------------------------------------------
double distanceL(struct position x1, struct position x2, double L)
{
	// Distance between two points
	// L = 2 => Euclidean
	double n;
	
	n = 0;

	for (int d = 0; d < x1.size; d++)
		n = n + pow(fabs(x1.x[d] - x2.x[d]), L);

	n = pow(n, 1 / L);
	return n;
}

// ---------------------------------------------------------------
int sign(double x)
{
	if (x == 0)	return 0;
	if (x < 0)	return -1;
	return 1;
}

// ===============================================================

// ==========================KISS.c===============================
/*
 A good pseudo-random numbers generator

 The idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
 x(n)=a*x(n-1)+1 mod 2^32
	 y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
	 z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
		 The y's are a shift register sequence on 32bit binary vectors
		 period 2^32-1;
	 The z's are a simple multiply-with-carry sequence with period
		 2^63+2^32-1.  The period of KISS is thus
		 2^32*(2^32-1)*(2^63+2^32-1) > 2^127
*/

static ulong kiss_x; //= 1;
static ulong kiss_y; //= 2;
static ulong kiss_z; //= 4;
static ulong kiss_w; //= 8;
static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;



void seed_rand_kiss(ulong seed)
{
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
}

ulong rand_kiss()
{
	kiss_x = kiss_x * 69069 + 1;
	kiss_y ^= kiss_y << 13;
	kiss_y ^= kiss_y >> 17;
	kiss_y ^= kiss_y << 5;
	kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
	kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
	kiss_z = kiss_w;
	kiss_w = kiss_m;
	kiss_carry = kiss_k >> 30;
	//printf("\n%f ",(double) (kiss_x + kiss_y + kiss_w));
	return kiss_x + kiss_y + kiss_w;
}
// ===============================================================

// ========================mersenne.c==============================
/*
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

	 1. Redistributions of source code must retain the above copyright
		notice, this list of conditions and the following disclaimer.

	 2. Redistributions in binary form must reproduce the above copyright
		notice, this list of conditions and the following disclaimer in the
		documentation and/or other materials provided with the distribution.

	 3. The names of its contributors may not be used to endorse or promote
		products derived from this software without specific prior written
		permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
	 ACM Transactions on Modeling and
	 Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
	 ``Mersenne Twister: a 623-dimensionally equidistributed
	   uniform pseudorandom number generator''
	 ACM Transactions on Modeling and
	 Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/
// initializes mt[NN] with a seed
void init_genrand64(unsigned long long seed)
{
	mt[0] = seed;
	for (mti = 1; mti < NN; mti++)
		mt[mti] = (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array64(unsigned long long init_key[],
	unsigned long long key_length)
{
	unsigned long long i, j, k;
	init_genrand64(19650218ULL);
	i = 1; j = 0;
	k = (NN > key_length ? NN : key_length);
	for (; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL))
			+ init_key[j] + j; /* non linear */
		i++; j++;
		if (i >= NN) { mt[0] = mt[NN - 1]; i = 1; }
		if (j >= key_length) j = 0;
	}
	for (k = NN - 1; k; k--) {
		mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL))
			- i; /* non linear */
		i++;
		if (i >= NN) { mt[0] = mt[NN - 1]; i = 1; }
	}

	mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0, 2^64-1]-interval */
unsigned long long genrand64_int64(void)
{
	int i;
	unsigned long long x;
	static unsigned long long mag01[2] = { 0ULL, MATRIX_A };

	if (mti >= NN) { /* generate NN words at one time */

		/* if init_genrand64() has not been called, */
		/* a default initial seed is used     */
		if (mti == NN + 1)
			init_genrand64(5489ULL);

		for (i = 0; i < NN - MM; i++) {
			x = (mt[i] & UM) | (mt[i + 1] & LM);
			mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
		}
		for (; i < NN - 1; i++) {
			x = (mt[i] & UM) | (mt[i + 1] & LM);
			mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];
		}
		x = (mt[NN - 1] & UM) | (mt[0] & LM);
		mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int)(x & 1ULL)];

		mti = 0;
	}

	x = mt[mti++];

	x ^= (x >> 29) & 0x5555555555555555ULL;
	x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
	x ^= (x << 37) & 0xFFF7EEE000000000ULL;
	x ^= (x >> 43);

	return x;
}

/* generates a random number on [0, 2^63-1]-interval */
long long genrand64_int63(void)
{
	return (long long)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
	return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
	return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void)
{
	return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
}
// ===============================================================
