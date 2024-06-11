//#include "stdio.h"
//
//#include "pso_structures.h"
//#include "sub_programs.h"
//
////#include "perf.c"
//
//// ===============================================================
//// PSO
//
//Result PSO(Param param, Problem pb, int randOption)
//{
//	Velocity aleaV;
//	int d;
//	double error;
//	double errorPrev;
//	Velocity expt1, expt2;
//	int g;
//	Velocity GX;
//	int index[S_max], indexTemp[S_max];
//	int initLinks;						// Flag to (re)init or not the information links
//	int iter;							// Iteration number (time step)
//	int iterBegin;
//	int length;
//	int LINKS[S_max][S_max];			// Information links
//	int m;
//	int noEval;
//	double normPX, normGX;
//	int noStop;
//	int outside;
//	double p;
//	Velocity PX;
//	Result R;
//	double r[D_max];
//	int rank;
//	Matrix RotatePX;
//	Matrix RotateGX;
//	int s0, s, s1;
//	int t;
//	double zz;
//
//	aleaV.size = pb.SS.D;
//	RotatePX.size = pb.SS.D;
//	RotateGX.size = pb.SS.D;
//
//	long double sqrtD = sqrt((long double)pb.SS.D);
//
//	// -----------------------------------------------------
//	// INITIALISATION
//	p = param.p;						// Probability threshold for random topology
//	R.SW.S = param.S;					// Size of the current swarm
//
//	// Position and velocity
//	for (s = 0; s < R.SW.S; s++)
//	{
//		R.SW.X[s].size = pb.SS.D;
//		R.SW.V[s].size = pb.SS.D;
//
//		for (d = 0; d < pb.SS.D; d++)
//		{
//			R.SW.X[s].x[d] = alea(pb.SS.minInit[d], pb.SS.maxInit[d], randOption); // Initial random Position
//		}
//
//		for (d = 0; d < pb.SS.D; d++)
//		{
//			R.SW.V[s].v[d] = (alea(pb.SS.min[d], pb.SS.max[d], randOption) - R.SW.X[s].x[d]) / 2; // Velocity calculate function
//		}
//
//		// Take quantisation into account
//		R.SW.X[s] = quantis(R.SW.X[s], pb.SS);
//		// printf("\n %1.20e, R.SW.X[s].x[0]");
//	}
//
//	// First evaluations
//	for (s = 0; s < R.SW.S; s++)
//	{
//		R.SW.X[s].f = perf(R.SW.X[s], pb.function, pb.SS, pb.objective);
//		R.SW.P[s] = R.SW.X[s];												// Best Position = current one
//		R.SW.P[s].improved = 0;												// No improvement
//		//printf("\n %1.20e",R.SW.X[s].f);	
//	}
//
//	// If the number max of evaluations is smaller than the swarm size,
//	// just keep evalMax particles, and finish
//	if (R.SW.S > pb.evalMax) R.SW.S = pb.evalMax;
//	R.nEval = R.SW.S;
//
//	// Find the best
//	R.SW.best = 0;
//	switch (param.stop)
//	{
//	case 2:
//		errorPrev = distanceL(R.SW.P[R.SW.best], pb.solution, 2);			// Distance to the wanted solution
//		break;
//	default:
//		errorPrev = R.SW.P[R.SW.best].f;									// "distance" to the wanted f value (objective)
//		break;
//	}
//
//	for (s = 1; s < R.SW.S; s++)
//	{
//		switch (param.stop)
//		{
//		case 2:
//			zz = distanceL(R.SW.P[R.SW.best], pb.solution, 2);
//			if (zz < errorPrev)
//			{
//				R.SW.best = s;
//				errorPrev = zz;
//			}
//			break;
//		default:
//			zz = R.SW.P[s].f;
//			if (zz < errorPrev)
//			{
//				R.SW.best = s;
//				errorPrev = zz;
//			}
//			break;
//		}
//	}
//
//	// Display the best
//	printf("Best value after init. %1.20e", errorPrev);
//
//	// printf("\n Position: \n")
//	// for(d=0; d<SS.D; d++) printf("%f", R.SW.P[R.SW.best].x[d]);
//
//	initLinks = 1;															// So that information links will beinitialized
//	// Note: It is also a flag saying "No improvement"
//	noStop = 0;
//	error = errorPrev;
//
//	// ---------------------------------------------- ITERATIONS
//	iter = 0; iterBegin = 0;
//	while (noStop == 0)
//	{
//		iter = iter + 1;
//
//		if (initLinks == 1)													// Random topology
//		{
//			// Who informs who, at random
//			for (s = 0; s < R.SW.S; s++)
//			{
//				for (s = 0; s < R.SW.S; s++)
//				{
//					if (alea(0, 1, randOption) < p) LINKS[m][s] = 1;					// Probabilistic method
//					else LINKS[m][s] = 0;
//				}
//			}
//			/*
//			 // Ring topology  (Just for test)
//			 for (s = 0; s < R.SW.S; s++)
//			 {
//				 for (m = 0; m < R.SW.S; m++)
//				 {
//					 LINKS[m][s] = 0;
//				 }
//			 }
//			 for (s = 0; s < R.SW.S-1; s++)
//			 {
//				 for (m = s+1; m < R.SW.S; m++)
//				 {
//					 LINKS[m][s] = 1;
//				 }
//			 }
//			 LINKS[ 0 ][R.SW.S-1]=1;
//			 */
//			 // Each particle informs itself
//			for (m = 0; m < R.SW.S; m++)
//			{
//				LINKS[m][m] = 1;
//			}
//		}
//
//		// The swarm MOVES
//		// printf("\n Iteration %i", iter);
//		for (s = 0; s < R.SW.S; s++) index[s] = s;
//
//		switch (param.randOrder)
//		{
//		case 1:	// Define a random permutation
//			length = R.SW.S;
//			for (s = 0; s < length; s++)
//				indexTemp[s] = index[s];
//
//			for (s = 0; s < R.SW.S; s++)
//			{
//				rank = alea_integer(0, length - 1, randOption);
//				index[s] = indexTemp[rank];
//				if (rank < length - 1)	// Compact
//				{
//					for (t = rank; t < length; t++)
//						indexTemp[t] = indexTemp[t + 1];
//				}
//				length = length - 1;
//			}
//			break;
//		default:
//			break;
//		}
//
//		for (s0 = 0; s0 < R.SW.S; s0++)	// For each particle ...
//		{
//			s = index[s0];
//			// ... find the first informant
//			s1 = 0;
//			while (LINKS[s1][s] == 0)
//				s1++;
//			if (s1 >= R.SW.S) s1 = s;
//
//			// Find the best informant
//			g = s1;
//			for (m = s1; m < R.SW.S; m++)
//			{
//				if (LINKS[m][s] == 1 && R.SW.P[m].f < R.SW.P[g].f)
//					g = m;
//			}
//
//			// ... compute the new velocity, and move
//
//			// Exploration tendency
//			for (d = 0; d < pb.SS.D; d++)
//			{
//				//R.SW.V[s].v[d] = param.w * R.SW.V[s].v[d];
//				R.SW.V[s].v[d] *= param.w ;
//			}
//
//			// Prepare Exploitation tendency p-x
//			for (d = 0; d < pb.SS.D; d++)
//			{
//				PX.v[d] = R.SW.P[s].x[d] - R.SW.X[s].x[d];
//			}
//			PX.size = pb.SS.D;
//
//			if (g != s)
//			{
//				for (d = 0; d < pb.SS.D; d++)	// g-x
//				{
//					GX.v[d] = R.SW.P[g].x[d] - R.SW.X[s].x[d];
//				}
//				GX.size = pb.SS.D;
//			}
//
//			// Option "non sentivity to rotation"			
//			if (param.rotation > 0)
//			{
//				normPX = normL(PX, 2);
//				if (g != s) normGX = normL(GX, 2);
//				if (normPX > 0)
//				{
//					RotatePX = matrixRotation(PX);
//				}
//
//				if (g != s && normGX > 0)
//				{
//					RotateGX = matrixRotation(GX);
//				}
//			}
//
//			// Exploitation tendencies
//			switch (param.rotation)
//			{
//			default:
//				for (d = 0; d < pb.SS.D; d++)
//				{
//					r[d] = alea(0, param.c, randOption);
//					R.SW.V[s].v[d] = R.SW.V[s].v[d] +
//						+r[d] * PX.v[d];
//				}
//
//
//				if (g != s)
//				{
//					for (d = 0; d < pb.SS.D; d++)
//					{
//						r[d] = alea(0, param.c, randOption);
//						//r[d]=param.c-r[d]; //********** TEST
//						R.SW.V[s].v[d] = R.SW.V[s].v[d]
//							+ r[d] * GX.v[d];
//					}
//				}
//
//				break;
//
//			case 1:
//				// First exploitation tendency
//				if (normPX > 0)
//				{
//					zz = param.c * normPX / sqrtD;
//					aleaV = aleaVector(pb.SS.D, zz, randOption);
//					expt1 = matrixVectProduct(RotatePX, aleaV);
//
//					for (d = 0; d < pb.SS.D; d++)
//					{
//						R.SW.V[s].v[d] = R.SW.V[s].v[d] + expt1.v[d];
//					}
//				}
//
//				// Second exploitation tendency
//				if (g != s && normGX > 0)
//				{
//					zz = param.c * normGX / sqrtD;
//					aleaV = aleaVector(pb.SS.D, zz, randOption);
//					expt2 = matrixVectProduct(RotateGX, aleaV);
//
//					for (d = 0; d < pb.SS.D; d++)
//					{
//						R.SW.V[s].v[d] = R.SW.V[s].v[d] + expt2.v[d];
//					}
//				}
//				break;
//			}
//
//			// Update the Position
//			for (d = 0; d < pb.SS.D; d++)
//			{
//				R.SW.X[s].x[d] = R.SW.X[s].x[d] + R.SW.V[s].v[d];
//			}
//
//			if (R.nEval >= pb.evalMax)
//			{
//				//error= fabs(error - pb.objective);
//				goto end;
//			}
//			// --------------------------
//			noEval = 1;
//
//			// Quantisation
//			R.SW.X[s] = quantis(R.SW.X[s], pb.SS);
//
//			switch (param.clamping)
//			{
//			case 0:	// No clamping AND no evaluation
//				outside = 0;
//
//				for (d = 0; d < pb.SS.D; d++)
//				{
//					if (R.SW.X[s].x[d] < pb.SS.min[d] || R.SW.X[s].x[d] > pb.SS.max[d])
//						outside++;
//				}
//
//				if (outside == 0)	// If inside, the Position is evaluated
//				{
//					R.SW.X[s].f =
//						perf(R.SW.X[s], pb.function, pb.SS, pb.objective);
//					R.nEval = R.nEval + 1;
//				}
//				break;
//
//			case 1:	// Set to the bounds, and v to zero
//				for (d = 0; d < pb.SS.D; d++)
//				{
//					if (R.SW.X[s].x[d] < pb.SS.min[d])
//					{
//						R.SW.X[s].x[d] = pb.SS.min[d];
//						R.SW.V[s].v[d] = 0;
//					}
//
//					if (R.SW.X[s].x[d] > pb.SS.max[d])
//					{
//						R.SW.X[s].x[d] = pb.SS.max[d];
//						R.SW.V[s].v[d] = 0;
//					}
//				}
//
//				R.SW.X[s].f = perf(R.SW.X[s], pb.function, pb.SS, pb.objective);
//				R.nEval = R.nEval + 1;
//				break;
//			}
//
//			// ... update the best previous Position
//			if (R.SW.X[s].f < R.SW.P[s].f)	// Improvement
//			{
//				R.SW.P[s] = R.SW.X[s];
//
//				// ... update the best of the bests
//				if (R.SW.P[s].f < R.SW.P[R.SW.best].f)
//				{
//					R.SW.best = s;
//				}
//			}
//		}	// End of "for (s0 = 0)..."
//
//		// Check if finished
//		switch (param.stop)
//		{
//		default:
//			error = R.SW.P[R.SW.best].f;
//			break;
//
//		case 2:
//			error = distanceL(R.SW.P[R.SW.best], pb.solution, 2);
//			break;
//		}
//		//error= fabs(error - pb.epsilon);
//
//		if (error < errorPrev)	// Improvement
//		{
//			initLinks = 0;
//		}
//		else			// No improvement
//		{
//			initLinks = 1;	// Information links will be	reinitialized	
//		}
//
//		if (param.initLink == 1) initLinks = 1 - initLinks;
//
//		errorPrev = error;
//	end:
//
//		switch (param.stop)
//		{
//		case 0:
//		case 2:
//			if (error > pb.epsilon && R.nEval < pb.evalMax)
//			{
//				noStop = 0;	// Won't stop
//			}
//			else
//			{
//				noStop = 1;	// Will stop
//			}
//			break;
//
//		case 1:
//			if (R.nEval < pb.evalMax)
//				noStop = 0;	// Won't stop
//			else
//				noStop = 1;	// Will stop
//			break;
//		}
//	} // End of "whild no stop ..."
//
//	// printf( "\n and the winner is ... %i", R.SW.best );			
//	// fprintf( f_stag, "\nEND" );
//	R.error = error;
//	return R;
//};
