/*
Standard PSO 2007
 Contact for remarks, suggestions etc.:
 Maurice.Clerc@WriteMe.com

 Last update
 2011-01-08 Fixed seed for the RNG (for reproducible results)
 2010-12-12 Compression spring
 2010-12-11 fitness structure, to easier cope with constraints
 2010-10-01 Perm function
 2010-09-25 Sine sine function
 2010-08-15 Gear train problem
 2010-06-15 Fixed a small bug in position initialisation for discrete problems (thanks to Yue,Shuai)
 2010-03-24 Lennard-Jones problem
 2010-02-06 A few functions of the CEC 2005 benchmark
 2010-01-04 Fixed wrong fitness evaluation for G3 (function code 10)
 2009-12-29 Random number generator KISS (option). For more reproducible results
 2009-07-12 The initialisation space may be smaller than the search space (for tests)
 2009-06-03 Fixed a small mistake about the first best position
 2009-05-05 Step function
 2009-04-19 Schaffer f6, 2D Goldstein-Price
 2009-03-31 A small network optimisation
 2009-03-12 Schwefel 2.2, Neumaier 3, G3 (constrained)
 2008-10-02 Two Schwefel functions
 2008-08-12 For information: save the best position over all runs
 2007-12-10 Warning about rotational invariance (valid here only on 2D)
 2007-11-22 stop criterion (option): distance to solution < epsilon
			and log_progress evaluation
 2007-11-21 Ackley function

  -------------------------------- Contributors
 The works and comments of the following persons have been taken
 into account while designing this standard.  Sometimes this is for
 including a feature, and sometimes for leaving out one.

 Auger, Anne
 Blackwell, Tim
 Bratton, Dan
 Clerc, Maurice
 Croussette, Sylvain
 Dattasharma, Abhi
 Eberhart, Russel
 Hansen, Nikolaus
 Keko, Hrvoje
 Kennedy, James
 Krohling, Renato
 Langdon, William
 Li, Wentao
 Liu, Hongbo
 Miranda, Vladimiro
 Poli, Riccardo
 Serra, Pablo
 Stickel, Manfred

 -------------------------------- Motivation
Quite often, researchers claim to compare their version of PSO
with the "standard one", but the "standard one" itself seems to vary!
Thus, it is important to define a real standard that would stay
unchanged for at least one year.
This PSO version does not intend to be the best one on the market
(in particular, there is no adaptation of the swarm size nor of the
coefficients). This is simply very near to the original version (1995),
with just a few improvements based on some recent works.
 --------------------------------- Metaphors
swarm: A team of communicating people (particles)
At each time step
	Each particle chooses a few informants at random, selects the best
	one from this set, and takes into account the information given by
	the chosen particle.
	If it finds no particle better than itself, then the "reasoning" is:
	"I am the best, so I just take my current velocity and my previous
	best position into account"
----------------------------------- Parameters/Options
clamping := true/false => whether to use clamping positions or not
randOrder:= true/false => whether to avoid the bias due to the loop
				on particles "for s = 1 to swarm_size ..." or not
rotation := true/false => whether the algorithm is sensitive
				to a rotation of the landscape or not
You may also modify the following ones, although suggested values
are either hard coded or automatically computed:
S := swarm size
K := maximum number of particles _informed_ by a given one
w := first cognitive/confidence coefficient
c := second cognitive/confidence coefficient
 ----------------------------------- Equations
For each particle and each dimension
Equation 1:	v(t+1) = w*v(t) + R(c)*(p(t)-x(t)) + R(c)*(g(t)-x(t))
Equation 2:	x(t+1) = x(t) + v(t+1)
where
v(t) := velocity at time t
x(t) := position at time t
p(t) := best previous position of the particle
g(t) := best position amongst the best previous positions
		of the informants of the particle
R(c) := a number coming from a random distribution, which depends on c
In this standard, the distribution is uniform on [0,c]
Note 1:
When the particle has no informant better than itself,
it implies p(t) = g(t)
Therefore, Equation 1 gets modified to:
v(t+1) = w*v(t) + R(c)*(p(t)-x(t))
Note 2:
When the "non sensitivity to rotation" option is activated
(p(t)-x(t)) (and (g(t)-x(t))) are replaced by rotated vectors,
so that the final DNPP (Distribution of the Next Possible Positions)
is not dependent on the system of co-ordinates.
 ----------------------------------- Information links topology
A lot of work has been done about this topic. The main result is this:
There is no "best" topology. Hence the random approach used here.
 ----------------------------------- Initialisation
Initial positions are chosen at random inside the search space
(which is supposed to be a hyperparallelepiped, and often even
a hypercube), according to a uniform distribution.
This is not the best way, but the one used in the original PSO.
Each initial velocity is simply defined as the half-difference of two
random positions. It is simple, and needs no additional parameter.
However, again, it is not the best approach. The resulting distribution
is not even uniform, as is the case for any method that uses a
uniform distribution independently for each component.
The mathematically correct approach needs to use a uniform
distribution inside a hypersphere. It is not very difficult,
and was indeed used in some PSO versions.  However, it is quite
different from the original one.
Moreover, it may be meaningless for some heterogeneous problems,
when each dimension has a different "interpretation".
------------------------------------ From SPSO-06 to SPSO-07
The main differences are:
1. option "non sensitivity to rotation of the landscape"
	Note: although theoretically interesting, this option is quite
		computer time consuming, and the improvement in result may
		only be marginal.
2. option "random permutation of the particles before each iteration"
	Note: same remark. Time consuming, no clear improvement
3. option "clamping position or not"
	Note: in a few rare cases, not clamping positions may induce an
	infinite run, if the stop criterion is the maximum number of
	evaluations

4. probability p of a particular particle being an informant of another
	particle. In SPSO-06 it was implicit (by building the random infonetwork)
	Here, the default value is directly computed as a function of (S,K),
	so that the infonetwork is exactly the same as in SPSO-06.
	However, now it can be "manipulated" ( i.e. any value can be assigned)

5. The search space can be quantised (however this algorithm is _not_
   for combinatorial problems)
Also, the code is far more modular. It means it is slower, but easier
to translate into another language, and easier to modify.
 ----------------------------------- Use
 Define the problem (you may add your own one in problemDef() and perf())
 Choose your options
 Run and enjoy!

 */


#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>

#define	D_max 114		// Max number of dimensions of the search space
#define R_max 500		// Max number of runs
#define	S_max 910		// Max swarm size
#define fMax 6				// Max number of constraints +1
#define zero  0			// 1.0e-30 // To avoid numerical instabilities

#define ulong unsigned long // To generate pseudo-random numbers with KISS
#define RAND_MAX_KISS ((unsigned long) 4294967295)

#define ERROR(a) {printf("\nERROR: %s\n", a);exit(-1);}

 // Structures
struct quantum
{
	double q[D_max];
	int size;
};

struct SS
{
	int D;
	double max[D_max];
	double maxInit[D_max];
	double min[D_max];
	double minInit[D_max];
	struct quantum q;		// Quantisation step size. 0 => continuous problem
};

struct param
{
	double c;		// Confidence coefficient
	int clamping;	// Position clamping or not
	int K;			// Max number of particles informed by a given one
	double p;		// Probability threshold for random topology	
	// (is actually computed as p(S,K) )
	int randOrder;	// Random choice of particles or not
	int rand; // 0 => use KISS. Any other value: use the standard C RNG
	int initLink; // How to re-init links
	int rotation;	// Sensitive to rotation or not
	int S;			// Swarm size
	int stop;		// Flag for stop criterion
	double w;		// Confidence coefficient
};

struct fitness
{
	int size;
	double f[fMax];
};

struct position
{
	double f;
	int improved;
	int size;
	double x[D_max];
};

struct velocity
{
	int size;
	double v[D_max];
};
struct problem
{
	double epsilon; 	// Admissible error
	int evalMax; 		// Maximum number of fitness evaluations
	int function; 		// Function code
	double objective; 	// Objective value
	// Solution position (if known, just for tests)	
	struct position solution;
	struct SS SS;		// Search space
};

struct swarm
{
	int best; 					// rank of the best particle
	struct position P[S_max];	// Previous best positions found by each particle
	int S; 						// Swarm size 
	struct velocity V[S_max];	// Velocities
	struct position X[S_max];	// Positions 
};

struct result
{
	double nEval; 		// Number of evaluations  
	struct swarm SW;	// Final swarm
	double error;		// Numerical result of the run
};
struct matrix 	// Useful for "non rotation sensitive" option
{
	int size;
	double v[D_max][D_max];
};

// Sub-programs
ulong	rand_kiss_07(); // For the pseudo-random number generator KISS
void	seed_rand_kiss_07(ulong seed);

double alea_07(double a, double b);
int alea_integer_07(int a, int b);
double alea_normal_07(double mean, double stdev);
struct fitness constraint_07(struct position x, int functCode, double eps);
double distanceL_07(struct position x1, struct position x2, double L);

double lennard_jones_07(struct position x);

struct velocity aleaVector_07(int D, double coeff);

struct matrix* matrixProduct_07(struct matrix* M1, struct matrix* M2);
struct matrix* matrixRotation_07(struct velocity* V);
struct velocity*	matrixVectProduct_07(struct matrix* M, struct velocity* V);

double normL_07(struct velocity v, double L);
double perf_07(struct position x, int function, struct SS SS, double objective);	// Fitness evaluation
struct position quantis_07(struct position x, struct SS SS);
struct problem problemDef_07(int functionCode);
struct result* PSO_07(struct param* param, struct problem* problem);
int sign_07(double x);

// Global variables
long double E;			// exp(1). Useful for some test functions
long double pi;			// Useful for some test functions
int randOption;
struct matrix reflex1;
long double sqrtD;

double shift;

// For Network problem
int bcsNb;
int btsNb;

// File(s);
FILE* f_run;
FILE* f_synth;

// =================mainfunc========================
int SPSO2007()
{

	struct position* bestBest; // Best position over all runs
	bestBest = (struct position*)malloc(sizeof(struct position));
	int d;			// Current dimension
	double error;			// Current error
	double errorMean;		// Average error
	double errorMin;		// Best result over all runs
	double errorMeanBest[R_max];
	double evalMean;		// Mean number of evaluations
	int functionCode;
	int i, j;
	int nFailure;		// Number of unsuccessful runs
	double logProgressMean = 0;
	struct param param;
	struct problem* pb;
	pb = (struct problem*)malloc(sizeof(struct problem));
	int run, runMax;
	struct result* result;
	result = (struct result*)malloc(sizeof(struct result));
	double successRate;
	double variance;

	f_run = fopen("f_run.txt", "w");
	f_synth = fopen("f_synth.txt", "w");

	E = exp((long double)1);
	pi = acos((long double)-1);


	// ----------------------------------------------- PROBLEM
	functionCode = 102;
	/* (see problemDef( ) for precise definitions)
	 0 Parabola (Sphere)
	 1 Griewank
	 2 Rosenbrock (Banana)
	 3 Rastrigin
	 4 Tripod (dimension 2)
	 5 Ackley
	 6 Schwefel
	 7 Schwefel 1.2
	 8 Schwefel 2.22
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
	 21 Compression spring
	 22 g03
	 23 Penalized

	  CEC 2005 benchmark  (no more than 30D. See cec2005data.c)
	 100 F1 (shifted Parabola/Sphere)
	 102 F6 (shifted Rosenbrock)
	 103 F9 (shifted Rastrigin)
	 104 F2 Schwefel
	 105 F7 Griewank  (NOT rotated)
	 106 F8 Ackley  (NOT rotated)

	 99 Test

							  */

	runMax = 30; // Numbers of runs
	if (runMax > R_max) runMax = R_max;


	// -----------------------------------------------------
	// PARAMETERS
	// * means "suggested value"		

	param.clamping = 1;
	// 0 => no clamping AND no evaluation. WARNING: the program
	// 				may NEVER stop (in particular with option move 20 (jumps)) 1
	// *1 => classical. Set to bounds, and velocity to zero

	param.initLink = 0; // 0 => re-init links after each unsuccessful iteration
	// 1 => re-init links after each successful iteration

	param.rand = 0; // 0 => Use KISS as random number generator. 
	// Any other value => use the standard C one

	param.randOrder = 1; // 0 => at each iteration, particles are modified
	//     always according to the same order 0..S-1
	//*1 => at each iteration, particles numbers are
	//		randomly permutated
	param.rotation = 0;
	// WARNING. Experimental code, completely valid only for dimension 2
	// 0 =>  sensitive to rotation of the system of coordinates
	// 1 => non sensitive (except side effects), 
	// 			by using a rotated hypercube for the probability distribution
	//			WARNING. Quite time consuming!

	param.stop = 0;	// Stop criterion
	// 0 => error < pb.epsilon
	// 1 => eval >= pb.evalMax		
	// 2 => ||x-solution|| < pb.epsilon

	// -------------------------------------------------------
	// Some information
	printf("\n Function %i ", functionCode);
	printf("\n (clamping, randOrder, rotation, stop_criterion) = (%i, %i, %i, %i)",
		param.clamping, param.randOrder, param.rotation, param.stop);
	if (param.rand == 0) printf("\n WARNING, I am using the RNG KISS");

	// =========================================================== 
	// RUNs

	// Initialize some objects
	*pb = problemDef_07(functionCode);

	// You may "manipulate" S, p, w and c
	// but here are the suggested values
	param.S = (int)(10 + 2 * sqrt(pb->SS.D));	// Swarm size
	//param.S=40;	
	if (param.S > S_max) param.S = S_max;

	printf("\n Swarm size %i", param.S);

	param.K = 3;
	param.p = 1 - pow(1 - (double)1 / (param.S), param.K);
	// (to simulate the global best PSO, set param.p=1)
	//param.p=1;

	// According to Clerc's Stagnation Analysis
	param.w = 1 / (2 * log((double)2)); // 0.721
	param.c = 0.5 + log((double)2); // 1.193

	// According to Poli's Sampling Distribution of PSOs analysis	
	//param.w = ??; // in [0,1[
	//param.c =  
	//    smaller than 12*(param.w*param.w-1)/(5*param.w -7); 


	printf("\n c = %f,  w = %f", param.c, param.w);
	randOption = param.rand; // Global variable
	//---------------
	sqrtD = sqrt((long double)pb->SS.D);

	// Define just once the first reflexion matrix
	if (param.rotation > 0)
	{
		reflex1.size = pb->SS.D;
		for (i = 0; i < pb->SS.D; i++)
		{
			for (j = 0; j < pb->SS.D; j++)
			{
				reflex1.v[i][j] = -2.0 / pb->SS.D;
			}
		}

		for (d = 0; d < pb->SS.D; d++)
		{
			reflex1.v[d][d] = 1 + reflex1.v[d][d];
		}
	}

	errorMean = 0;
	evalMean = 0;
	nFailure = 0;
	//------------------------------------- RUNS
	seed_rand_kiss_07(1294404794);	// For reproducible results, if using KISS
	for (run = 0; run < runMax; run++)
	{
		//srand (clock () / 100);	// May improve pseudo-randomness            
		// 
		result = PSO_07(&param, pb);
		error = result->error;

		if (error > pb->epsilon) // Failure
		{
			nFailure = nFailure + 1;
		}

		// Memorize the best (useful if more than one run)
		if (run == 0) *bestBest = result->SW.P[result->SW.best];
		else
			if (error < bestBest->f) *bestBest = result->SW.P[result->SW.best];

		// Result display
		printf("\nRun %i. Eval %f. Error %e ", run + 1, result->nEval, error);
		printf("  Success  %.2f%% \n", 100 * (1 - (double)nFailure / (run + 1)));
		//	for (d=0;d<pb.SS.D;d++) printf(" %f",bestBest.x[d]);

			// Save result
		fprintf(f_run, "\n%i %.0f %e ", run + 1, result->nEval, error);
		for (d = 0; d < pb->SS.D; d++) fprintf(f_run, " %f", bestBest->x[d]);

		// Compute/store some statistical information
		if (run == 0)
			errorMin = error;
		else if (error < errorMin)
			errorMin = error;
		evalMean = evalMean + result->nEval;
		errorMean = errorMean + error;
		errorMeanBest[run] = error;
		logProgressMean = logProgressMean - log(error);
	}		// End loop on "run"

	// ---------------------END 
	// Display some statistical information
	evalMean = evalMean / (double)runMax;
	errorMean = errorMean / (double)runMax;
	logProgressMean = logProgressMean / (double)runMax;

	printf("\n Eval. (mean)= %f", evalMean);
	printf("\n Error (mean) = %1.5e", errorMean);
	// Variance
	variance = 0;

	for (run = 0; run < runMax; run++)
		variance = variance + pow(errorMeanBest[run] - errorMean, 2);

	variance = sqrt(variance / runMax);
	printf("\n Std. dev. %f", variance);
	printf("\n Log_progress (mean) = %f", logProgressMean);
	// Success rate and minimum value
	printf("\n Failure(s) %i", nFailure);
	successRate = 100 * (1 - nFailure / (double)runMax);
	printf("\n Success rate = %.2f%%", successRate);

	//if (run > 1)
	{
		printf("\n Best min value = %1.5e", errorMin);
		printf("\nPosition of the optimum: ");
		for (d = 0; d < pb->SS.D; d++) printf(" %f", bestBest->x[d]);
	}

	// Save	
	fprintf(f_synth, " %1.5e %1.5e %.0f%% %f   ",
		errorMean, variance, successRate, evalMean);
	for (d = 0; d < pb->SS.D; d++) fprintf(f_synth, " %f", bestBest->x[d]);
	fprintf(f_synth, "\n");

	return 0; // End of main program
}

// ===============================================================
// PSO
struct result* PSO_07(struct param* param, struct problem* pb)
{
	struct velocity* aleaV;
	aleaV = (struct velocity*)malloc(sizeof(struct velocity));
	int d;
	double error = 0;
	struct result* R;
	R = (struct result*)malloc(sizeof(struct result));
	double errorPrev;
	struct velocity *expt1, *expt2;
	expt1 = (struct velocity*)malloc(sizeof(struct velocity));
	expt2 = (struct velocity*)malloc(sizeof(struct velocity));
	int g;
	struct velocity* GX;
	GX = (struct velocity*)malloc(sizeof(struct velocity));
	int index[S_max], indexTemp[S_max];
	int initLinks;	// Flag to (re)init or not the information links
	int iter; 		// Iteration number (time step)
	int iterBegin;
	int length;
	int m;
	int noEval;
	double normPX, normGX;
	int noStop;
	int outside;
	double p;
	struct velocity* PX;
	PX = (struct velocity*)malloc(sizeof(struct velocity));
	double r[D_max];
	int rank;
	struct matrix* RotatePX;
	RotatePX = (struct matrix*)malloc(sizeof(struct matrix));
	struct matrix* RotateGX;
	RotateGX = (struct matrix*)malloc(sizeof(struct matrix));
	int s0, s, s1;
	int t;
	double zz;

	aleaV->size = pb->SS.D;
	RotatePX->size = pb->SS.D;
	RotateGX->size = pb->SS.D;
	// -----------------------------------------------------
	// INITIALISATION
	p = param->p; // Probability threshold for random topology
	R->SW.S = param->S; // Size of the current swarm

	// Position and velocity
	for (s = 0; s < R->SW.S; s++)
	{
		R->SW.X[s].size = pb->SS.D;
		R->SW.V[s].size = pb->SS.D;

		for (d = 0; d < pb->SS.D; d++)
		{
			R->SW.X[s].x[d] = alea_07(pb->SS.minInit[d], pb->SS.maxInit[d]);
		}

		for (d = 0; d < pb->SS.D; d++)
		{
			R->SW.V[s].v[d] =
				(alea_07(pb->SS.min[d], pb->SS.max[d]) - R->SW.X[s].x[d]) / 2;
		}
		// Take quantisation into account
		R->SW.X[s] = quantis_07(R->SW.X[s], pb->SS);
		//printf("\n %1.20e",R.SW.X[s].x[0]);	
	}
	//ERROR("stop SPSO");	

		// First evaluations
	for (s = 0; s < R->SW.S; s++)
	{
		R->SW.X[s].f =
			perf_07(R->SW.X[s], pb->function, pb->SS, pb->objective);

		R->SW.P[s] = R->SW.X[s];	// Best position = current one
		R->SW.P[s].improved = 0;	// No improvement
		//printf("\n %1.20e",R.SW.X[s].f);	
	}

	// If the number max of evaluations is smaller than 
	// the swarm size, just keep evalMax particles, and finish
	if (R->SW.S > pb->evalMax) R->SW.S = pb->evalMax;
	R->nEval = R->SW.S;

	// Find the best
	R->SW.best = 0;
	switch (param->stop)
	{
	default:
		errorPrev = R->SW.P[R->SW.best].f; // "distance" to the wanted f value (objective)
		break;

	case 2:
		errorPrev = distanceL_07(R->SW.P[R->SW.best], pb->solution, 2); // Distance to the wanted solution
		break;
	}

	for (s = 1; s < R->SW.S; s++)
	{
		switch (param->stop)
		{
		default:
			zz = R->SW.P[s].f;
			if (zz < errorPrev)
			{
				R->SW.best = s;
				errorPrev = zz;
			}
			break;

		case 2:
			zz = distanceL_07(R->SW.P[R->SW.best], pb->solution, 2);
			if (zz < errorPrev)
			{
				R->SW.best = s;
				errorPrev = zz;
			}
			break;
		}
	}
	// Display the best
	printf(" Best value after init. %1.20e ", errorPrev);

	//	printf( "\n Position :\n" );
	//	for ( d = 0; d < SS.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );

	initLinks = 1;		// So that information links will beinitialized
	// Note: It is also a flag saying "No improvement"
	noStop = 0;
	error = errorPrev;
	// ---------------------------------------------- ITERATIONS
	iter = 0; iterBegin = 0;

	int** LINKS;	// Information links
	LINKS = (int**)malloc(R->SW.S * sizeof(int*));
	for (s = 0; s < R->SW.S; s++)
	{
		LINKS[s] = (int*)malloc(R->SW.S * sizeof(int));
	}


	while (noStop == 0)
	{
		iter = iter + 1;

		if (initLinks == 1)	// Random topology
		{
			// Who informs who, at random
			for (s = 0; s < R->SW.S; s++)
			{
				for (m = 0; m < R->SW.S; m++)
				{
					if (alea_07(0, 1) < p) LINKS[m][s] = 1;	// Probabilistic method
					else LINKS[m][s] = 0;
				}
			}
			/*
			 // Ring topology  (Just for test)
			 for (s = 0; s < R.SW.S; s++)
			 {
				 for (m = 0; m < R.SW.S; m++)
				 {
					 LINKS[m][s] = 0;
				 }
			 }
			 for (s = 0; s < R.SW.S-1; s++)
			 {
				 for (m = s+1; m < R.SW.S; m++)
				 {
					 LINKS[m][s] = 1;
				 }
			 }
			 LINKS[ 0 ][R.SW.S-1]=1;
			 */
			 // Each particle informs itself
			for (m = 0; m < R->SW.S; m++)
			{
				LINKS[m][m] = 1;
			}
		}

		// The swarm MOVES
		//printf("\nIteration %i",iter);
		for (s = 0; s < R->SW.S; s++)  index[s] = s;

		switch (param->randOrder)
		{
		default:
			break;

		case 1: //Define a random permutation
			length = R->SW.S;
			for (s = 0; s < length; s++) indexTemp[s] = index[s];

			for (s = 0; s < R->SW.S; s++)
			{
				rank = alea_integer_07(0, length - 1);
				index[s] = indexTemp[rank];
				if (rank < length - 1)	// Compact
				{
					for (t = rank; t < length; t++)
						indexTemp[t] = indexTemp[t + 1];
				}
				length = length - 1;
			}
			break;
		}
		for (s0 = 0; s0 < R->SW.S; s0++)	// For each particle ...
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
			PX->size = pb->SS.D;

			if (g != s)
			{
				for (d = 0; d < pb->SS.D; d++) // g-x
				{
					GX->v[d] = R->SW.P[g].x[d] - R->SW.X[s].x[d];
				}
				GX->size = pb->SS.D;
			}

			// Option "non sentivity to rotation"				
			if (param->rotation > 0)
			{
				normPX = normL_07(*PX, 2);
				if (g != s) normGX = normL_07(*GX, 2);
				if (normPX > 0)
				{
					*RotatePX = *matrixRotation_07(PX);
				}

				if (g != s && normGX > 0)
				{
					*RotateGX = *matrixRotation_07(GX);
				}
			}

			// Exploitation tendencies
			switch (param->rotation)
			{
			default:
				for (d = 0; d < pb->SS.D; d++)
				{
					r[d] = alea_07(0, param->c);
					R->SW.V[s].v[d] = R->SW.V[s].v[d] +
						+r[d] * PX->v[d];
				}


				if (g != s)
				{
					for (d = 0; d < pb->SS.D; d++)
					{
						r[d] = alea_07(0, param->c);
						//r[d]=param.c-r[d]; //********** TEST
						R->SW.V[s].v[d] = R->SW.V[s].v[d]
							+ r[d] * GX->v[d];
					}
				}

				break;

			case 1:
				// First exploitation tendency
				if (normPX > 0)
				{
					zz = param->c * normPX / sqrtD;
					*aleaV = aleaVector_07(pb->SS.D, zz);
					expt1 = matrixVectProduct_07(RotatePX, aleaV);

					for (d = 0; d < pb->SS.D; d++)
					{
						R->SW.V[s].v[d] = R->SW.V[s].v[d] + expt1->v[d];
					}
				}

				// Second exploitation tendency
				if (g != s && normGX > 0)
				{
					zz = param->c * normGX / sqrtD;
					*aleaV = aleaVector_07(pb->SS.D, zz);
					expt2 = matrixVectProduct_07(RotateGX, aleaV);

					for (d = 0; d < pb->SS.D; d++)
					{
						R->SW.V[s].v[d] = R->SW.V[s].v[d] + expt2->v[d];
					}
				}
				break;
			}

			// Update the position
			for (d = 0; d < pb->SS.D; d++)
			{
				R->SW.X[s].x[d] = R->SW.X[s].x[d] + R->SW.V[s].v[d];
			}

			if (R->nEval >= pb->evalMax)
			{
				//error= fabs(error - pb.objective);
				goto end;
			}
			// --------------------------
			noEval = 1;

			// Quantisation
			R->SW.X[s] = quantis_07(R->SW.X[s], pb->SS);

			switch (param->clamping)
			{
			case 0:	// No clamping AND no evaluation
				outside = 0;

				for (d = 0; d < pb->SS.D; d++)
				{
					if (R->SW.X[s].x[d] < pb->SS.min[d] || R->SW.X[s].x[d] > pb->SS.max[d])
						outside++;
				}

				if (outside == 0)	// If inside, the position is evaluated
				{
					R->SW.X[s].f =
						perf_07(R->SW.X[s], pb->function, pb->SS, pb->objective);
					R->nEval = R->nEval + 1;
				}
				break;

			case 1:	// Set to the bounds, and v to zero
				for (d = 0; d < pb->SS.D; d++)
				{
					if (R->SW.X[s].x[d] < pb->SS.min[d])
					{
						R->SW.X[s].x[d] = pb->SS.min[d];
						R->SW.V[s].v[d] = 0;
					}

					if (R->SW.X[s].x[d] > pb->SS.max[d])
					{
						R->SW.X[s].x[d] = pb->SS.max[d];
						R->SW.V[s].v[d] = 0;
					}
				}

				R->SW.X[s].f = perf_07(R->SW.X[s], pb->function, pb->SS, pb->objective);
				R->nEval = R->nEval + 1;
				break;
			}

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
		}			// End of "for (s0=0 ...  "	
		// Check if finished
		switch (param->stop)
		{
		default:
			error = R->SW.P[R->SW.best].f;
			break;

		case 2:
			error = distanceL_07(R->SW.P[R->SW.best], pb->solution, 2);
			break;
		}
		//error= fabs(error - pb.epsilon);

		if (error < errorPrev)	// Improvement
		{
			initLinks = 0;
		}
		else			// No improvement
		{
			initLinks = 1;	// Information links will be	reinitialized	
		}

		if (param->initLink == 1) initLinks = 1 - initLinks;

		errorPrev = error;
	end:

		switch (param->stop)
		{
		case 0:
		case 2:
			if (error > pb->epsilon && R->nEval < pb->evalMax)
			{
				noStop = 0;	// Won't stop
			}
			else
			{
				noStop = 1;	// Will stop
			}
			break;

		case 1:
			if (R->nEval < pb->evalMax)
				noStop = 0;	// Won't stop
			else
				noStop = 1;	// Will stop
			break;
		}

	} // End of "while nostop ...

	 printf( "\n and the winner is ... %i", R->SW.best );			
	 fprintf(f_synth, "\nEND" );

	R->error = error;
	return R;
}


// ===========================================================
double alea_07(double a, double b)
{				// random number (uniform distribution) in  [a b]
	// randOption is a global parameter
	double r;
	if (randOption == 0)
		r = a + (double)rand_kiss_07() * (b - a) / RAND_MAX_KISS;
	else
		r = a + (double)rand() * (b - a) / RAND_MAX;
	//printf("\nalea 944 r %f",r);
	return r;
}

// ===========================================================
int alea_integer_07(int a, int b)
{				// Integer random number in [a b]
	int ir;
	double r;

	r = alea_07(0, 1);
	ir = (int)(a + r * (b + 1 - a));

	if (ir > b)	ir = b;

	return ir;
}

// ===========================================================
double alea_normal_07(double mean, double std_dev)
{
	/*
	 Use the polar form of the Box-Muller transformation to obtain a pseudo
	 random number from a Gaussian distribution
	 */
	double x1, x2, w, y1;
	// double y2;

	do
	{
		x1 = 2.0 * alea_07(0, 1) - 1.0;
		x2 = 2.0 * alea_07(0, 1) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);

	w = sqrt(-2.0 * log(w) / w);
	y1 = x1 * w;
	// y2 = x2 * w;
	if (alea_07(0, 1) < 0.5) y1 = -y1;
	y1 = y1 * std_dev + mean;
	return y1;
}


// =============================================================
struct velocity aleaVector_07(int D, double coeff)
{
	struct velocity V;
	int d;
	int i;
	int K = 2; // 1 => uniform distribution in a hypercube
	// 2 => "triangle" distribution
	double rnd;

	V.size = D;

	for (d = 0; d < D; d++)
	{
		rnd = 0;
		for (i = 1; i <= K; i++) rnd = rnd + alea_07(0, 1);
		V.v[d] = rnd * coeff / K;
	}

	return V;
}

// ===========================================================
double normL_07(struct velocity v, double L)
{   // L-norm of a vector
	int d;
	double n;

	n = 0;

	for (d = 0; d < v.size; d++)
		n = n + pow(fabs(v.v[d]), L);

	n = pow(n, 1 / L);
	return n;
}

// ===========================================================
double distanceL_07(struct position x1, struct position x2, double L)
{  // Distance between two positions
	// L = 2 => Euclidean	
	int d;
	double n;

	n = 0;

	for (d = 0; d < x1.size; d++)
		n = n + pow(fabs(x1.x[d] - x2.x[d]), L);

	n = pow(n, 1 / L);
	return n;
}

//============================================================
struct matrix* matrixProduct_07(struct matrix* M1, struct matrix* M2)
{
	// Two square matrices of same size
	struct matrix* Product;
	Product = (struct matrix*)malloc(sizeof(struct matrix));
	int D;
	int i, j, k;
	double sum;
	D = M1->size;
	for (i = 0; i < D; i++)
	{
		for (j = 0; j < D; j++)
		{
			sum = 0;
			for (k = 0; k < D; k++)
			{
				sum = sum + M1->v[i][k] * M2->v[k][j];
			}
			Product->v[i][j] = sum;
		}
	}
	Product->size = D;
	return Product;
}
//=========================================================
struct matrix* matrixRotation_07(struct velocity* V)
{
	/*
	 Define the matrice of the rotation V' => V
	 where V'=(1,1,...1)*normV/sqrt(D)  (i.e. norm(V') = norm(V) )

	 */
	struct velocity* B;
	B = (struct velocity*)malloc(sizeof(struct velocity));
	int i, j, d, D;
	double normB, normV, normV2;
	//struct matrix reflex1; // Global variable
	struct matrix *reflex2;
	reflex2 = (struct matrix*)malloc(sizeof(struct matrix));
	struct matrix *rotateV;
	rotateV = (struct matrix*)malloc(sizeof(struct matrix));
	double temp;

	D = V->size;
	normV = normL_07(*V, 2); normV2 = normV * normV;
	reflex2->size = D;

	// Reflection relatively to the vector V'=(1,1, ...1)/sqrt(D)	
	// norm(V')=1
	// Has been computed just once  (global matrix reflex1)	

	//Define the "bisectrix" B of (V',V) as an unit vector
	B->size = D;
	temp = normV / sqrtD;

	for (d = 0; d < D; d++)
	{
		B->v[d] = V->v[d] + temp;
	}
	normB = normL_07(*B, 2);

	if (normB > 0)
	{
		for (d = 0; d < D; d++)
		{
			B->v[d] = B->v[d] / normB;
		}
	}

	// Reflection relatively to B
	for (i = 0; i < D; i++)
	{
		for (j = 0; j < D; j++)
		{
			reflex2->v[i][j] = -2 * B->v[i] * B->v[j];
		}
	}

	for (d = 0; d < D; d++)
	{
		reflex2->v[d][d] = 1 + reflex2->v[d][d];
	}

	// Multiply the two reflections
	// => rotation				
	rotateV = matrixProduct_07(reflex2, &reflex1);
	return rotateV;

}
//==========================================================
struct velocity* matrixVectProduct_07(struct matrix* M, struct velocity* V)
{
	struct velocity* Vp;
	Vp = (struct velocity*)malloc(sizeof(struct velocity));
	int d, j;
	int Dim;
	double sum;
	Dim = V->size;
	for (d = 0; d < Dim; d++)
	{
		sum = 0;
		for (j = 0; j < Dim; j++)
		{
			sum = sum + M->v[d][j] * V->v[j];
		}
		Vp->v[d] = sum;
	}
	Vp->size = Dim;
	return Vp;
}
// ===========================================================
int sign_07(double x)
{
	if (x == 0)	return 0;
	if (x < 0)	return -1;
	return 1;
}

// ===========================================================
struct position quantis_07(struct position x, struct SS SS)
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
			//qd = qd * (SS.max[d] - SS.min[d]) / 2;	      
			quantx.x[d] = qd * floor(0.5 + x.x[d] / qd);
		}
	}
	return quantx;
}
//================================================== KISS
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

static ulong kiss_x = 1;
static ulong kiss_y = 2;
static ulong kiss_z = 4;
static ulong kiss_w = 8;
static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;


void seed_rand_kiss_07(ulong seed)
{
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
}

ulong rand_kiss_07()
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
// ===========================================================
double perf_07(struct position x, int function, struct SS SS, double objective)
{
	// Evaluate the fitness value for the particle of rank s 
	double beta;
	double c;
	int d;
	double DD;
	int  k;
	int n;
	struct fitness ff;
	double f, p, xd, x1, x2, x3, x4;
	double s11, s12, s21, s22;
	double sum1, sum2;
	double t0, tt, t1;
	double u;
	struct position xs;
	//----------------#include "cec2005data.c"----------------
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

	// Shifted Schwefel (F2 CEC 2005)
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
	//--------------------------------------------------------

	/*
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
	*/
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

	xs = x;
	switch (function)
	{
		//---------------------#include "cec2005.c"---------------
	case 100: // Parabola (Sphere) CEC 2005
		f = -450;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d] - offset_0[d];
			f = f + xd * xd;
		}
		break;

	case 102:  // Rosenbrock 

		for (d = 0; d < xs.size; d++)
		{
			xs.x[d] = xs.x[d] - offset_2[d];
		}

		f = 390;

		for (d = 1; d < xs.size; d++)
		{
			tt = xs.x[d - 1] - 1;
			f = f + tt * tt;

			tt = xs.x[d - 1] * xs.x[d - 1] - xs.x[d];
			f = f + 100 * tt * tt;
		}

		break;

	case 103: // Rastrigin 
		for (d = 0; d < xs.size; d++)
		{
			xs.x[d] = xs.x[d] - offset_3[d];
		}
		f = -330;
		k = 10;

		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d];
			f = f + xd * xd - k * cos(2 * pi * xd);
		}
		f = f + xs.size * k;
		break;

	case 104: // Schwefel (F2)
		for (d = 0; d < xs.size; d++)
		{
			xs.x[d] = xs.x[d] - offset_4[d];
		}

		f = -450;
		for (d = 0; d < xs.size; d++)
		{
			sum2 = 0.0;
			for (k = 0; k <= d; k++)
			{
				sum2 += xs.x[k];
			}
			f += sum2 * sum2;
		}
		break;


	case 105: // Griewank. WARNING: in the CEC 2005 benchmark it is rotated
		sum1 = 0.0;
		sum2 = 1.0;
		f = -180;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d] - offset_5[d];
			sum1 += xd * xd;
			sum2 *= cos(xd / sqrt(1.0 + d));
		}
		f = f + 1.0 + sum1 / 4000.0 - sum2;
		break;

	case 106: // Ackley 
		f = -140;

		sum1 = 0.0;
		sum2 = 0.0;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d] - offset_6[d];
			sum1 += xd * xd;
			sum2 += cos(2.0 * pi * xd);
		}
		sum1 = -0.2 * sqrt(sum1 / xs.size);
		sum2 /= xs.size;
		f = f + 20.0 + E - 20.0 * exp(sum1) - exp(sum2);
		break;
		//--------------------------------------------------------

	case 0:		// Parabola (Sphere)
		f = 0;

		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d];
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
		s11 = (1.0 - sign_07(x1)) / 2;
		s12 = (1.0 + sign_07(x1)) / 2;
		s21 = (1.0 - sign_07(x2)) / 2;
		s22 = (1.0 + sign_07(x2)) / 2;

		//f = s21 * (fabs (x1) - x2); // Solution on (0,0)
		f = s21 * (fabs(x1) + fabs(x2 + 50)); // Solution on (0,-50)  
		f = f + s22 * (s11 * (1 + fabs(x1 + 50) +
			fabs(x2 - 50)) + s12 * (2 +
				fabs(x1 - 50) +
				fabs(x2 - 50)));
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
		f = lennard_jones_07(xs);
		break;

	case 18: // Gear train
		f = pow(1. / 6.931 - x.x[0] * x.x[1] / (x.x[2] * x.x[3]), 2);
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
		// Constraints
		ff = constraint_07(xs, function, 0);

		if (ff.f[1] > 0) { c = 1 + ff.f[1]; f = f * c * c * c; }
		if (ff.f[2] > 0) { c = 1 + ff.f[1]; f = f * c * c * c; }
		if (ff.f[3] > 0) { c = 1 + ff.f[3]; f = f * c * c * c; }
		if (ff.f[4] > 0) { c = 1 + pow(10, 10) * ff.f[4]; f = f * c * c * c; }
		if (ff.f[5] > 0) { c = 1 + pow(10, 10) * ff.f[5]; f = f * c * c * c; }
		break;

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
	case 99: // Test

		f = pow(fabs(xs.x[0]) - 5, 2) + pow(fabs(xs.x[1]) - 5, 2);

		break;

		f = 0;
		for (d = 0; d < x.size; d++)
		{
			f = f + pow((int)(xs.x[d] + 0.5), 2);
		}
		break;
		// Periodic
		// On [-10,10]^2
		// 49 local minima all with minimum value 1 
		// and one global minimum located at (0,0) with f = 0.9.

		x1 = xs.x[0]; x2 = xs.x[1];
		f = 1 + pow(sin(x1), 2) + pow(sin(x2), 2) - 0.1 * exp(-x1 * x1 - x2 * x2);
		break;

		// Modified Rosenbrock
		// On [-5,5]^2
		// Two global minima (0.3412, 0.1164), and (1, 1), on which f=0

		x1 = xs.x[0]; x2 = xs.x[1];

		f = 100 * pow(x2 - x1 * x1, 2) + pow(6.4 * pow(x2 - 0.5, 2) - x1 - 0.6, 2);

		break;
		// Wood's  
		// on [-10,10]^4
		// Solution (1,1,1,1)
		x1 = xs.x[0]; x2 = xs.x[1];
		x3 = xs.x[2]; x4 = xs.x[3];

		f = 100 * pow(x2 - x1 * x1, 2) + pow(1 - x1, 2) + 90 * pow(x4 - x3 * x3, 2)
			+ pow(1 - x3, 2) + 10.1 * (pow(x2 - 1, 1) + pow(x4 - 1, 2))
			+ 19.8 * (x2 - 1) * (x4 - 1);

		break;
		x1 = xs.x[0]; x2 = xs.x[1];
		sum1 = x1 * x1 + x2 * x2;
		f = 0.5 + (pow(sin(sqrt(sum1)), 2) - 0.5) / (1 + 0.001 * sum1 * sum1);
		break;



		f = 0;

		for (d = 0; d < x.size; d++)
		{
			sum1 = 0;
			for (k = 0; k < d + 1; k++)
			{
				sum1 = sum1 + xs.x[d];
			}

			f = f + sum1 * sum1;
		}

		break;

		f = 1.e6 * xs.x[0] * xs.x[0];

		for (d = 1; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f + xd * xd;
		}

		break;

		// 2D Peaks function
		x1 = xs.x[0];
		x2 = xs.x[1];

		f = 3 * (1 - x1) * (1 - x1) * exp(-x1 * x1 - (x2 + 1) * (x2 + 1))
			- 10 * (x1 / 5 - pow(x1, 3) - pow(x2, 5)) * exp(-x1 * x1 - x2 * x2)
			- (1. / 3) * exp(-(x1 + 1) * (x1 + 1) - x2 * x2);

		break;

		// Quartic
		f = 0;
		for (d = 0; d < x.size; d++)
		{
			xd = xs.x[d];
			f = f + (d + 1) * pow(xd, 4) + alea_normal_07(0, 1);
		}

		break;

		x1 = xs.x[0]; x2 = xs.x[1];
		f = (1 - x1) * (1 - x1) * exp(-x1 * x1 - (x2 + 1) * (x2 + 1)) - (x1 - x1 * x1 * x1 - pow(x2, 5)) * exp(-x1 * x1 - x2 * x2);
		f = -f; // To minimise
		break;

		xd = xs.x[0];
		f = xd * (xd + 1) * cos(xd);
		break;

	default:
		break;
	}

	return fabs(f - objective);
}

//==========================================================
struct fitness constraint_07(struct position x, int functCode, double epsConstr)
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

	}

	return ff;
}

struct problem problemDef_07(int functionCode)
{
	int d;
	struct problem pb;

	int nAtoms; // For Lennard-Jones problem
	static double lennard_jones[14] =
	{ -1, -3, -6, -9.103852, -12.71, -16.505384,-19.821489,-24.113360,-28.422532,
	-32.77,-37.97,-44.33,-47.84,-52.32 };


	pb.function = functionCode;
	pb.epsilon = 0.00000;	// Acceptable error (default). May be modified below
	pb.objective = 0;       // Objective value (default). May be modified below

	// Define the solution point, for test
	// NEEDED when param.stop = 2 
	// i.e. when stop criterion is distance_to_solution < epsilon
	for (d = 0; d < 30; d++)
	{
		pb.solution.x[d] = 0;
	}

	// ------------------ Search space
	switch (pb.function)
	{

	case 0:			// Parabola
		pb.SS.D = 30;//  Dimension							

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -5.12; // -100
			pb.SS.max[d] = 5.12;	// 100
			pb.SS.q.q[d] = 0;	// Relative quantisation, in [0,1].   
		}

		pb.evalMax = 20000;// Max number of evaluations for each run
		pb.epsilon = 0.9; // 1e-3;	
		pb.objective = 0;

		// For test purpose, the initialisation space may be different from
		// the search space. If so, just modify the code below

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d]; // May be a different value
			pb.SS.minInit[d] = pb.SS.min[d]; // May be a different value
		}

		break;
		// ------------------cec2005pb--------------------
	case 100: // CEC 2005 F1
		pb.SS.D = 30;//30; 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;

		}
		pb.evalMax = pb.SS.D * 10000;
		pb.epsilon = 0.000001;	//Acceptable error
		pb.objective = -450;       // Objective value

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 102:		// Rosenbrock. CEC 2005 F6
		pb.SS.D = 30;

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100; pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;

		}

		pb.evalMax = pb.SS.D * 10000;
		pb.epsilon = 0.01;	//0.01 Acceptable error
		pb.objective = 390;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 103:// CEC 2005 F9, Rastrigin
		pb.SS.D = 30;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -5.12;
			pb.SS.max[d] = 5.12;
			pb.SS.q.q[d] = 0;

		}
		pb.epsilon = 0.01; // 0.01;	// Acceptable error
		pb.objective = -330;       // Objective value
		pb.evalMax = pb.SS.D * 10000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 104:// CEC 2005 F2  Schwefel
		pb.SS.D = 10;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;

		}
		pb.epsilon = 0.00001;	// Acceptable error
		pb.objective = -450;       // Objective value
		pb.evalMax = pb.SS.D * 10000;

		pb.epsilon = 0.0001; pb.evalMax = 100000;


		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;


	case 105:// CEC 2005 F7  Griewank (NON rotated)
		pb.SS.D = 10;	 // 10 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -600;
			pb.SS.max[d] = 600;
			pb.SS.q.q[d] = 0;

		}
		pb.epsilon = 0.01;	//Acceptable error
		pb.objective = -180;       // Objective value
		pb.evalMax = pb.SS.D * 10000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;


	case 106:// CEC 2005 F8 Ackley (NON rotated)

		pb.SS.D = 30; // 10;	 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -32;
			pb.SS.max[d] = 32;
			pb.SS.q.q[d] = 0;
		}
		pb.epsilon = 0.0001;	// Acceptable error
		pb.objective = -140;       // Objective value
		pb.evalMax = pb.SS.D * 10000;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;
		// -----------------------------------------------


	case 1:		// Griewank
		pb.SS.D = 30; //30;	

		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -600;
			pb.SS.max[d] = 600;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 30000;
		pb.epsilon = 0.00; //0.001;
		pb.objective = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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
		pb.epsilon = 0.9;
		pb.evalMax = 20000; //2.e6;  // 40000 
		pb.objective = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		pb.evalMax = 20000; //3200; 
		pb.epsilon = 0.9;
		pb.objective = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.minInit[d] = pb.SS.min[d];
			pb.SS.maxInit[d] = pb.SS.max[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 5: // Ackley
		pb.SS.D = 30;
		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -30; // 32
			pb.SS.max[d] = 30;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 30000;
		pb.epsilon = 0.01;
		pb.objective = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 8: // Schwefel 2.22
		pb.SS.D = 30;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 100000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

		break;

	case 11: // Network
		btsNb = 19; bcsNb = 2;
		pb.SS.D = bcsNb * btsNb + 2 * bcsNb;
		pb.objective = 0;
		for (d = 0; d < bcsNb * btsNb; d++) // Binary representation. 1 means: there is a link
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 1;
		}

		for (d = bcsNb * btsNb; d < pb.SS.D; d++) // 2D space for the BSC positions
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 20; //15;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 5000;
		pb.objective = 0;
		pb.epsilon = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 12: // Schwefel
		pb.SS.D = 30;
		pb.objective = -418.98288727243369 * pb.SS.D;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -500;
			pb.SS.max[d] = 500;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 60000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 17: // Lennard-Jones
		nAtoms = 6; // in {2, ..., 15}
		pb.SS.D = 3 * nAtoms;
		pb.objective = lennard_jones[nAtoms - 2];
		pb.evalMax = 5000 + 3000 * nAtoms * (nAtoms - 1); // Empirical rule
		pb.epsilon = 1.e-6;


		//pb.SS.D=3*21; pb.objective=-81.684;	
		//pb.SS.D=3*27; pb.objective=-112.87358;
		//pb.SS.D=3*38; pb.objective=-173.928427;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -2;
			pb.SS.max[d] = 2;
			pb.SS.q.q[d] = 0;
		}

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 18: // Gear train
		// solution (16,19,43,49) and equivalent (like (19,16,49,43)
		// Success rate 9%
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 21: // Compression spring
		pb.SS.D = 3;

		pb.SS.min[0] = 1; pb.SS.max[0] = 70; pb.SS.q.q[0] = 1;
		pb.SS.min[1] = 0.6; pb.SS.max[1] = 3; pb.SS.q.q[1] = 0;
		pb.SS.min[2] = 0.207; pb.SS.max[2] = 0.5; pb.SS.q.q[2] = 0.001;

		pb.evalMax = 20000;
		pb.epsilon = 1.e-10;
		pb.objective = 2.6254214578;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	case 22: // g03
		pb.SS.D = 10;
		pb.objective = 0;
		pb.epsilon = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = 0;
			pb.SS.max[d] = 1;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 50000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
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

		pb.evalMax = 50000; // 250000;	 
		pb.epsilon = 0;
		pb.objective = 0;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;


	case 99: // Test

		pb.SS.D = 2;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 30000;
		pb.objective = 0;
		pb.epsilon = 0.00;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

		break;


		pb.SS.D = 10;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -100;
			pb.SS.max[d] = 100;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 100000;
		pb.objective = 0;
		pb.epsilon = 0.00;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

		break;
		pb.SS.D = 2;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 50000;
		pb.objective = 0.9;
		pb.epsilon = 0.001;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;
		pb.SS.D = 4;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 50000;
		pb.objective = 0.0;
		pb.epsilon = 0.001;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}


		break;


		// 2D Peaks function
		pb.SS.D = 2;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -3;
			pb.SS.max[d] = 3;
			pb.SS.q.q[d] = 0;
		}

		pb.evalMax = 50000;
		pb.objective = -6.551133;
		pb.epsilon = 0.001;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;
		// Quartic
		pb.SS.D = 50;
		pb.objective = 0;
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}
		pb.evalMax = 25000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

		break;


		pb.SS.D = 2;	// Dimension
		pb.objective = -2;

		pb.SS.min[0] = -2;
		pb.SS.max[0] = 2;
		pb.SS.q.q[0] = 0;
		pb.SS.min[1] = -3;
		pb.SS.max[1] = 3;
		pb.SS.q.q[1] = 0;

		pb.evalMax = 10000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}

		break;
		pb.SS.D = 1;	// Dimension
		// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] = -10;
			pb.SS.max[d] = 10;
			pb.SS.q.q[d] = 0;
		}
		pb.objective = -1000; // Just a sure too small value for the above search space
		pb.evalMax = 1000;

		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.maxInit[d] = pb.SS.max[d];
			pb.SS.minInit[d] = pb.SS.min[d];
		}
		break;

	default:
		break;
	}

	pb.SS.q.size = pb.SS.D;
	return pb;
}

// -----------------------lennard_jones---------------------
double lennard_jones_07(struct position x)
{
	/*
		This is for black-box optimisation. Therefore, we are not supposed to know
		that there are some symmetries. That is why the dimension of the problem is
		3*nb_of_atoms, as it could be 3*nb_of_atoms-6
	*/
	int d;
	int dim = 3;
	double dist;
	double f;
	int i, j;
	int nPoints = x.size / dim;
	struct position x1, x2;
	double zz;

	x1.size = dim; x2.size = dim;

	f = 0;
	for (i = 0; i < nPoints - 1; i++)
	{
		for (d = 0; d < dim; d++)  x1.x[d] = x.x[3 * i + d];
		for (j = i + 1; j < nPoints; j++)
		{
			for (d = 0; d < dim; d++)  x2.x[d] = x.x[3 * j + d];

			dist = distanceL_07(x1, x2, 2);
			zz = pow(dist, -6);
			f = f + zz * (zz - 1);
		}
	}
	f = 4 * f;
	return f;
}
