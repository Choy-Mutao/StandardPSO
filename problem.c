// 问题的边界条件的定义

#include "pso_structures.h"

//===================================================
Problem problemDef(int functionCode)
{
	Problem pb;
	pb.function = functionCode;
	pb.epsilon = 0.00000;			// Acceptable error (default). May be modified below
	pb.objective = 0;				// Objective value (default). May be modified below

	pb.SS.q.size = pb.SS.D;
	switch (functionCode)
	{
	case 0:
		pb = Parabola();
		break;
	default:
		break;
	}

	return pb;
};

Problem Parabola()
{
	int d = 0;
	Problem pb;

	pb.SS.D = 30;
	for (d = 0; d < pb.SS.D; d++)
	{
		pb.SS.min[d] = -5.12;		// -100
		pb.SS.max[d] = 5.12;		// 100
		pb.SS.q.q[d] = 0;			// Relative quantisation, in [0,1]
	}

	pb.evalMax = 20000;				// Max number of evaluations for each run
	pb.epsilon = 0.9;				// 1e-3
	pb.objective = 0;

	// For test purpose, the initialisation space may be different from the search space. If so, just modify the code below

	for (d = 0; d < pb.SS.D; d++)
	{
		pb.SS.maxInit[d] = pb.SS.max[d]; // May be a different value
		pb.SS.minInit[d] = pb.SS.min[d]; // May be a different value
	}
	
	return pb;
}

Problem Griewank()
{
	int d = 0;
	Problem pb;
	pb.SS.D = 30;

	// Boundariew
	for (d = 0; d < pb.SS.D; d++)
	{
		pb.SS.min[d] = -600;
		pb.SS.max[d] = 600;
		pb.SS.q.q[d] = 0;
	}

	pb.evalMax = 30000;
	pb.epsilon = 0.00; // 0.001
	pb.objective = 0;

	for ( d = 0; d < pb.SS.D; d++)
	{
		pb.SS.maxInit[d] = pb.SS.max[d];
		pb.SS.minInit[d] = pb.SS.min[d];
	}

	return pb;
}