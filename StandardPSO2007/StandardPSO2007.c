///*
//Standard PSO 2007
// Contact for remarks, suggestions etc.:
// Maurice.Clerc@WriteMe.com
//
// Last update
// 2011-01-08 Fixed seed for the RNG (for reproducible results)
// 2010-12-12 Compression spring
// 2010-12-11 fitness structure, to easier cope with constraints
// 2010-10-01 Perm function
// 2010-09-25 Sine sine function
// 2010-08-15 Gear train Problem
// 2010-06-15 Fixed a small bug in Position initialisation for discrete problems (thanks to Yue,Shuai)
// 2010-03-24 Lennard-Jones Problem
// 2010-02-06 A few functions of the CEC 2005 benchmark
// 2010-01-04 Fixed wrong fitness evaluation for G3 (function code 10)
// 2009-12-29 Random number generator KISS (option). For more reproducible results
// 2009-07-12 The initialisation space may be smaller than the search space (for tests)
// 2009-06-03 Fixed a small mistake about the first best Position
// 2009-05-05 Step function
// 2009-04-19 Schaffer f6, 2D Goldstein-Price
// 2009-03-31 A small network optimisation
// 2009-03-12 Schwefel 2.2, Neumaier 3, G3 (constrained)
// 2008-10-02 Two Schwefel functions
// 2008-08-12 For information: save the best Position over all runs
// 2007-12-10 Warning about rotational invariance (valid here only on 2D)
// 2007-11-22 stop criterion (option): distance to solution < epsilon
//			and log_progress evaluation
// 2007-11-21 Ackley function
//
//  -------------------------------- Contributors
// The works and comments of the following persons have been taken
// into account while designing this standard.  Sometimes this is for
// including a feature, and sometimes for leaving out one.
//
// Auger, Anne
// Blackwell, Tim
// Bratton, Dan
// Clerc, Maurice
// Croussette, Sylvain
// Dattasharma, Abhi
// Eberhart, Russel
// Hansen, Nikolaus
// Keko, Hrvoje
// Kennedy, James
// Krohling, Renato
// Langdon, William
// Li, Wentao
// Liu, Hongbo
// Miranda, Vladimiro
// Poli, Riccardo
// Serra, Pablo
// Stickel, Manfred
//
// -------------------------------- Motivation
//Quite often, researchers claim to compare their version of PSO
//with the "standard one", but the "standard one" itself seems to vary!
//Thus, it is important to define a real standard that would stay
//unchanged for at least one year.
//
//This PSO version does not intend to be the best one on the market
//(in particular, there is no adaptation of the swarm size nor of the
//coefficients).
//
//This is simply very near to the original version (1995),
//with just a few improvements based on some recent works.
//
// --------------------------------- Metaphors
//swarm: A team of communicating people (particles)
//At each time step
//	Each particle chooses a few informants at random, selects the best
//	one from this set, and takes into account the information given by
//	the chosen particle.
//	If it finds no particle better than itself, then the "reasoning" is:
//	"I am the best, so I just take my current velocity and my previous
//	best Position into account"
//
//----------------------------------- Parameters/Options
//clamping := true/false => whether to use clamping positions or not
//randOrder:= true/false => whether to avoid the bias due to the loop
//				on particles "for s = 1 to swarm_size ..." or not
//rotation := true/false => whether the algorithm is sensitive
//				to a rotation of the landscape or not
//You may also modify the following ones, although suggested values
//are either hard coded or automatically computed:
//S := swarm size
//K := maximum number of particles _informed_ by a given one
//w := first cognitive/confidence coefficient
//c := second cognitive/confidence coefficient
//
// ----------------------------------- Equations
//For each particle and each dimension
//Equation 1:	v(t+1) = w*v(t) + R(c)*(p(t)-x(t)) + R(c)*(g(t)-x(t))
//Equation 2:	x(t+1) = x(t) + v(t+1)
//where
//v(t) := velocity at time t
//x(t) := Position at time t
//p(t) := best previous Position of the particle
//g(t) := best Position amongst the best previous positions
//		of the informants of the particle
//R(c) := a number coming from a random distribution, which depends on c
//In this standard, the distribution is uniform on [0,c]
//
//Note 1:
//When the particle has no informant better than itself,
//it implies p(t) = g(t)
//Therefore, Equation 1 gets modified to:
//v(t+1) = w*v(t) + R(c)*(p(t)-x(t))
//
//Note 2:
//When the "non sensitivity to rotation" option is activated
//(p(t)-x(t)) (and (g(t)-x(t))) are replaced by rotated vectors,
//so that the final DNPP (Distribution of the Next Possible Positions)
//is not dependent on the system of co-ordinates.
//
// ----------------------------------- Information links topology
//A lot of work has been done about this topic. The main result is this:
//There is no "best" topology. Hence the random approach used here.
//
// ----------------------------------- Initialisation
//Initial positions are chosen at random inside the search space
//(which is supposed to be a hyperparallelepiped, and often even
//a hypercube), according to a uniform distribution.
//
//This is not the best way, but the one used in the original PSO.
//
//Each initial velocity is simply defined as the half-difference of two
//random positions. It is simple, and needs no additional parameter.
//However, again, it is not the best approach. The resulting distribution
//is not even uniform, as is the case for any method that uses a
//uniform distribution independently for each component.
//
//The mathematically correct approach needs to use a uniform
//distribution inside a hypersphere. It is not very difficult,
//and was indeed used in some PSO versions.
//
//However, it is quite different from the original one.
//
//Moreover, it may be meaningless for some heterogeneous problems,
//when each dimension has a different "interpretation".
//
//------------------------------------ From SPSO-06 to SPSO-07
//The main differences are:
//1. option "non sensitivity to rotation of the landscape"
//	Note: although theoretically interesting, this option is quite
//		computer time consuming, and the improvement in result may
//		only be marginal.
//2. option "random permutation of the particles before each iteration"
//	Note: same remark. Time consuming, no clear improvement
//3. option "clamping Position or not"
//	Note: in a few rare cases, not clamping positions may induce an
//	infinite run, if the stop criterion is the maximum number of
//	evaluations
//
//4. probability p of a particular particle being an informant of another
//	particle. In SPSO-06 it was implicit (by building the random infonetwork)
//	Here, the default value is directly computed as a function of (S,K),
//	so that the infonetwork is exactly the same as in SPSO-06.
//	However, now it can be "manipulated" ( i.e. any value can be assigned)
//
//5. The search space can be quantised (however this algorithm is _not_
//   for combinatorial problems)
//
//Also, the code is far more modular. It means it is slower, but easier
//to translate into another language, and easier to modify.
// ----------------------------------- Use
// Define the Problem (you may add your own one in problemDef() and perf())
// Choose your options
// Run and enjoy!
//
// */
//
//#include "stdio.h"
//#include "math.h"
//#include <stdlib.h>
//#include <time.h>
//
//#include "kiss_rng.h"
//#include "PSO.c"
//#include "Problem.c"
//
//#define R_max 500	// Max number of runs(Iterations)
//#define zero 0		// 1.0e-30 // To avoid numerical instablities
//
//
//#define ERROR(a) {printf("\n Error: %s\n", a);exit(-1);}
//
//// Global variables
//long double	E;
////long double	pi;
//int randOption;
//struct Matrix reflex1;
//long double sqrtD;
//
//double shift;
//
//// For network Problem
//int bcsNb;
//int btsNb;
//
//// File(s);
//FILE* f_run;
//FILE* f_synth;
//
//Param DefaultParam()
//{
//	Param param;
//
//	// -----------------------------------------------------
//	// PARAMETERS
//	// * means "suggested value"		
//
//	param.clamping = 1;
//	// 0 => no clamping AND no evaluation. WARNING: the program
//	// 				may NEVER stop (in particular with option move 20 (jumps)) 1
//	// *1 => classical. Set to bounds, and velocity to zero
//
//	param.initLink = 0; // 0 => re-init links after each unsuccessful iteration
//	// 1 => re-init links after each successful iteration
//
//	param.rand = 0; // 0 => Use KISS as random number generator. 
//	// Any other value => use the standard C one
//
//	param.randOrder = 1; // 0 => at each iteration, particles are modified
//	//     always according to the same order 0..S-1
//	//*1 => at each iteration, particles numbers are
//	//		randomly permutated
//	param.rotation = 0;
//	// WARNING. Experimental code, completely valid only for dimension 2
//	// 0 =>  sensitive to rotation of the system of coordinates
//	// 1 => non sensitive (except side effects), 
//	// 			by using a rotated hypercube for the probability distribution
//	//			WARNING. Quite time consuming!
//
//	param.stop = 0;	// Stop criterion
//	// 0 => error < pb.epsilon
//	// 1 => eval >= pb.evalMax		
//	// 2 => ||x-solution|| < pb.epsilon
//
//	return param;
//};
//
// // =================================================
//int main_2007()
//{
//	struct Position bestBest;	// Best Position over all runs
//	int d;						// Current dimension
//	double error;				// Current error
//	double errorMean;			// Average error
//	double errorMin;			// Best result over all runs
//	double errorMeanBest[R_max];
//	double evalMean;			// Mean number of evaluations
//	int functionCode;
//	int i, j;
//	int nFailure;				// Number of unsuccessful runs
//	double logProgressMean = 0;
//	struct Problem pb;
//	int run, runMax;
//	Result result;
//	double successRate;
//	double variance;
//
//	//f_run = fopen("f_run.txt", "w");
//	//f_synth = fopen("f_synth.txt", "w");
//	errno_t err1 = fopen_s(&f_run, "f_run.txt", "w");
//	if (err1 != 0)
//	{
//		perror("Error opening file");
//		return 1;
//	}
//	errno_t err2 = fopen_s(&f_synth, "f_synth.txt", "w");
//	if (err2 != 0)
//	{
//		perror("Error opening file");
//		return 1;
//	}
//
//	E = exp((long double)1);
//	//pi = acos((long double)-1);
//
//	// PROBLEM_CODE
//	functionCode = 102;
//
//	runMax = 30;				// Numbers of runs
//	if (runMax > R_max) runMax = R_max;
//
//	Param param = DefaultParam();
//	// -------------------------------------------------------
//	// Some information
//	printf("\n Function %i ", functionCode);
//	printf("\n (clamping, randOrder, rotation, stop_criterion) = (%i, %i, %i, %i)",
//		param.clamping, param.randOrder, param.rotation, param.stop);
//	if (param.rand == 0) printf("\n WARNING, I am using the RNG KISS");
//
//	// =========================================================== 
//	// RUNs
//
//	// Initialize some objects
//	pb = problemDef(functionCode);
//
//	// You may "manipulate" S, p, w and c
//	// but here are the suggested values
//	param.S = (int)(10 + 2 * sqrt(pb.SS.D));		// Swarm size
//	// param.S = 40;
//	if (param.S > S_max) param.S = S_max;
//
//	printf("\n Swarm size %i", param.S);
//
//	param.K = 3;
//	param.p = 1 - pow(1 - (double)1 / (param.S), param.K);
//	// (to simulate the global best PSO, set param.p = 1)
//	// param.P = 1;
//
//	// According to Clerc's Stagnation Analysis
//	param.w = 1 / (2 * log((double)2)); // 0.721
//	param.c = 0.5 + log((double)2); // 1.193
//
//	// According to Poli's sampling Distribution of PSOs analysis
//	// param.w = ?? // in [0,1]
//	// param.c = smaller than 12*(param.w*param.w-1)/(5*param.w - 7);
//
//	printf("\n c = %f, w = %f", param.c, param.w);
//	randOption = param.rand; // Global variable
//	// --------------------------------------------------------------
//	sqrtD = sqrt((long double)pb.SS.D);
//
//	// Define just once the first reflection Matrix;
//	if (param.rotation > 0)
//	{
//		reflex1.size = pb.SS.D;
//		for (i = 0; i < pb.SS.D; i++)
//		{
//			for (j = 0; j < pb.SS.D; j++) 
//			{
//				reflex1.v[i][j] = -2.0 / pb.SS.D;
//			}
//		}
//
//		for (d = 0; d < pb.SS.D; d++)
//		{
//			reflex1.v[d][d] = 1 + reflex1.v[d][d];
//		}
//	}
//
//	errorMean = 0;
//	evalMean = 0;
//	nFailure = 0;
//	//------------------------------------- RUNS
//
//	seed_rand_kiss(1294404794);		// For reproducible results, if using KISS
//	for (run = 0; run < runMax; run++)
//	{
//		// srand (clock()/100)		// May improve psedo-randomness
//		result = PSO(param, pb, 0);
//		error = result.error;
//
//		if (error > pb.epsilon)		// Failure
//		{
//			nFailure = nFailure + 1;
//		}
//
//		// Memorize the best (useful if more than one run)
//		if (run == 0) bestBest = result.SW.P[result.SW.best];
//		else
//			if (error < bestBest.f) bestBest = result.SW.P[result.SW.best];
//
//		// Result display
//		printf("\n Run %i. Eval %f. Error %e", run + 1, result.nEval, error);
//		printf("Success %.2f%% \n", 100 * (1 - (double)nFailure / (run + 1)));
//		// for(d=0; d<pb.SS.D; d++) printf("%f, bestBest.x[d]");
//
//		// Save result
//		fprintf(f_run, "\n %i %.0f %e", run + 1, result.nEval, error);
//		for (d = 0; d < pb.SS.D; d++) fprintf(f_run, "%f", bestBest.x[d]);
//
//		// Compute / store statictical information
//		if (run == 0)
//			errorMin = error;
//		else if (error < errorMin)
//			errorMin = error;
//		evalMean = evalMean + result.nEval;
//		errorMean = errorMean + error;
//		errorMeanBest[run] = error;
//		logProgressMean = logProgressMean - log(error);
//	}		// End loop on "run"
//
//	// ---------------------END 
//	// Display some statistical information
//	evalMean = evalMean / (double)runMax;
//	errorMean = errorMean / (double)runMax;
//	logProgressMean = logProgressMean / (double)runMax;
//
//	printf("\n Eval. (mean)= %f", evalMean);
//	printf("\n Error (mean) = %1.5e", errorMean);
//
//	// Variance
//	variance = 0;
//
//	for (run = 0; run < runMax; run++)
//		variance = variance + pow(errorMeanBest[run] - errorMean, 2);
//
//	variance = sqrt(variance / runMax);
//	printf("\n Log_progress (mean) = %f", logProgressMean);
//	// Success rate and minimum value
//	printf("\n Failure(s) %i", nFailure);
//	successRate = 100 * (1 - nFailure / (double)runMax);
//	printf("\n Success rate = %.2f%%", successRate);
//
//	//if (run > 1)
//	{
//		printf("\n Best min value = %1.5e", errorMin);
//		printf("\nPosition of the optimum: ");
//		for (d = 0; d < pb.SS.D; d++) printf(" %f", bestBest.x[d]);
//	}
//
//	// Save	
//	fprintf(f_synth, " %1.5e %1.5e %.0f%% %f   ",
//		errorMean, variance, successRate, evalMean);
//	for (d = 0; d < pb.SS.D; d++) fprintf(f_synth, " %f", bestBest.x[d]);
//	fprintf(f_synth, "\n");
//
//	return 0; // End of main program
//}