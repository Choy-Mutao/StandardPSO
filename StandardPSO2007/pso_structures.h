#ifndef PSO_STRUCTURES_H
#define  PSO_STRUCTURES_H

#define D_max 114	// Max number of dimensions of the search space
#define S_max 910	// Max sarm size
#define fMax 6		// Max number of constraints +1

// Structures
typedef struct Quantum
{
	double q[D_max];
	int size;
} Quantum;

typedef struct SS
{
	int D;
	double max[D_max];
	double maxInit[D_max];
	double min[D_max];
	double minInit[D_max];
	Quantum q;	// Quantisation step size .0 => continuous problem
} SS;

typedef struct Param
{
	double c;		// Confidence coefficient
	double w;		// Confidence coefficient

	int clamping;	// Position clamping or not
	int K;			// Max number of particles informed by a given one

	double p;		// Probability threshold for random topology
	// (is actually computed as p(S,K))

	int randOrder;	// Random choice of particles or not
	int rand;		// 0 => use Kiss. Any other value: use the standard C RNG
	int initLink;	// How to re-init links
	int rotation;	// Swarm size
	int S;			// Swarm size
	int stop;		// Flag for stop criterion
} Param;

typedef struct Fitness
{
	int size;
	double f[fMax];
} Fitness;

typedef struct Position
{
	double f;
	int improved;
	int size;
	double x[D_max];
} Position;

typedef struct Velocity
{
	int size;
	double v[D_max];
} Velocity;

typedef struct Problem
{
	double epsilon;		// Admissible error
	int evalMax;		// Maximum number of fitness evaluations
	int function;		// Function code

	double objective;	// Objective value
	// Solution Position(if konwn, just for tests)
	Position solution;
	SS SS;		// Search space
} Problem;

typedef struct Swarm
{
	int best;					// rank of the best particle
	Position P[S_max];	// Previous best positions found by each particle
	int S;						// Swarm size
	Velocity	V[S_max];	// Velocities
	Position X[S_max];	// Positions
} Swarm;

typedef struct Result
{
	double nEval;		// Number of evaluations
	Swarm SW;	// Final swarm
	double error;		// Numberical result of the run
} Result;

typedef struct Matrix			// Useful for "non rotation sensitive" option
{
	int size;
	double v[D_max][D_max];
} Matrix;

#endif // !PSO_STRUCTURES_H
