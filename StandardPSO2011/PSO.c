struct result PSO (struct param param, struct problem pb) 
{  
	int d; 
	double error;   
	double errorPrev;
	int g;  
	struct position Gr;
	struct vector GX={0};   
	int index[S_max];     
	int initLinks;	// Flag to (re)init or not the information links
	int iter; 		// Iteration number (time step)
	int iterBegin;
	int LINKS[S_max][S_max];	// Information links
	int m; 
	int noStop;
	double out;
	double p;
	struct vector PX;	
	struct vectorList qRand;
	struct result R;
	double rad;
	int randCase;
	int s0, s,s1; 
	struct vector V1,V2; 
	double w1,w2,w3;
	double xMax,xMin;
	//struct position XPrev;
	double zz;

	PX.size=pb.SS.D; 
	GX.size=pb.SS.D; 
	V1.size=pb.SS.D;
	V2.size=pb.SS.D;
	Gr.size=pb.SS.D;
	// -----------------------------------------------------
	// INITIALISATION
	p=param.p; // Probability threshold for random topology
	R.SW.S = param.S; // Size of the current swarm
randCase=param.BW[2]; 
if(randCase>300)  randCase=3;
if(randCase<-200) randCase=-2;

	// Position and velocity
	for (s = 0; s < R.SW.S; s++)   
	{
		R.SW.X[s].size = pb.SS.D;
		R.SW.V[s].size = pb.SS.D;
	}

	if(pb.SS.normalise>0) {xMin=0; xMax=pb.SS.normalise;} // [0,normalise]^D

switch(randCase)
{
	default:
	for (s = 0; s < R.SW.S; s++)   
	{
		
		for (d = 0; d < pb.SS.D; d++)  
		{  
			if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}
			R.SW.X[s].x[d] = alea (xMin,xMax,randCase);		
			R.SW.V[s].v[d] = alea( xMin-R.SW.X[s].x[d], xMax-R.SW.X[s].x[d],randCase ); 
			// So that  xMin<= x+v <= xMax
		}
	}
	break;
	
	case 1:
	case 2:
	qRand=quasiRand(pb.SS.D,R.SW.S,randCase);
	
		for (s = 0; s < R.SW.S; s++)   
	{
		for (d = 0; d < pb.SS.D; d++)  
		{  
			if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}
			R.SW.X[s].x[d] = xMin+(xMax-xMin)*qRand.V[s].v[d];		
			R.SW.V[s].v[d] = alea( xMin-R.SW.X[s].x[d], xMax-R.SW.X[s].x[d],randCase ); 
		}
	}	
	break;
}

	// Take quantisation into account
	if(pb.SS.quantisation==1)
	{
		for (s = 0; s < R.SW.S; s++)
		{  
			R.SW.X[s] = quantis (R.SW.X[s], pb.SS);
		}
	}

	// First evaluations
	errMax=0;
	errMin=infinity;
	for (s = 0; s < R.SW.S; s++) 
	{	
		R.SW.X[s].f = perf (R.SW.X[s], pb.function,pb.SS,pb.objective);
		R.SW.P[s] = R.SW.X[s];	// Best position = current one
	}

//
	// If the number max of evaluations is smaller than 
	// the swarm size, just keep evalMax particles, and finish
	if (R.SW.S>pb.evalMax) R.SW.S=pb.evalMax;	
	R.nEval = R.SW.S;

	// Find the best
	R.SW.best = 0;

	errorPrev =R.SW.P[R.SW.best].f; // "distance" to the wanted f value (objective)

	for (s = 1; s < R.SW.S; s++)     
	{
		zz=R.SW.P[s].f;
		if (zz < errorPrev)
		{
			R.SW.best = s;
			errorPrev=zz;
		} 
	}

	// Display the best
//	printf( "\nBest value after init. %1.20e ", errorPrev ); printf(" ");
//		printf( "\n Position :\n" );
//		for ( d = 0; d < pb.SS.D; d++ ) printf( " %f", R.SW.P[R.SW.best].x[d] );

	// ---------------------------------------------- ITERATIONS	
	noStop = 0;	
	error=errorPrev;	
	iter=0; iterBegin=0;
	initLinks = 1;		// So that information links will be initialized

	// Each particle informs itself
	for (m = 0; m < R.SW.S; m++) LINKS[m][m] = 1;	

	while (noStop == 0) 
	{	
		iter=iter+1;
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

	switch(param.BW[3])
	{
		default:
		case 0:
			for (s = 0; s < R.SW.S; s++)  index[s]=s;  // No random permutation
		break;
		
		case 1:
	aleaIndex(index, R.SW.S,randCase); // Random numbering of the particles
	break;

	}	

		if (initLinks==1)	// Modify topology
		{
			switch(param.topology)
			{
				default: // case 0. As in SPSO 2007
					// Who informs who, at random
					for (s = 0; s < R.SW.S; s++)
				{				
					for (m = 0; m < R.SW.S; m++)
					{		    
						if(m==s) continue;		

						if (alea (0, 1,randCase)<p)
						{ 
							LINKS[m][s] = 1;	// Probabilistic method
						}
						else LINKS[m][s] = 0;
					}	
				}				 
				break;

				case 1: // Random ring
					// Init to zero (no link)
					for (s = 0; s < R.SW.S; s++) 
				{	
					for (m = 0; m < R.SW.S; m++)    LINKS[m][s] = 0;
				}

				// Information links (bidirectional ring)
				for (s = 0; s < R.SW.S-1; s++)
				{	
					LINKS[index[s]][index[s+1]] = 1;
					LINKS[index[s+1]][index[s]] = 1;		
				}

				LINKS[index[0]][index[R.SW.S-1]] = 1;
				LINKS[index[R.SW.S-1]][index[0]] = 1;

				// Each particle informs itself
				for (m = 0; m < R.SW.S; m++) LINKS[m][m] = 1;	 
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
		for (s0 = 0; s0 < R.SW.S; s0++)	// For each particle ...
		{		
			s=index[s0];
			// ... find the first informant
			s1 = 0;    
			while (LINKS[s1][s] == 0)	s1++;		
		
			if (s1 >= R.SW.S)	s1 = s;

			// Find the best informant			
			g = s1;	
			for (m = s1; m < R.SW.S; m++) 
			{	    
				if (LINKS[m][s] == 1 && R.SW.P[m].f < R.SW.P[g].f)
					g = m;
			}	

			//.. compute the new velocity, and move

			// Exploration tendency
			for (d = 0; d < pb.SS.D; d++)
			{
				R.SW.V[s].v[d]=param.w *R.SW.V[s].v[d];
			}

			// Prepare Exploitation tendency  p-x
			for (d = 0; d < pb.SS.D; d++)
			{
				PX.v[d]= R.SW.P[s].x[d] - R.SW.X[s].x[d];
			}		
		
			if(g!=s) // If the particle is not its own local best, prepare g-x
			{
				for (d = 0; d < pb.SS.D; d++) 
				{
					GX.v[d]= R.SW.P[g].x[d] - R.SW.X[s].x[d];
				}
			}
			else // this is the best particle. Define another random "previous best" 
				// May be used or not, though (see below)
			{
//	if(param.BW[1]==1) // WARNING: un-comment this line largely modify the performances
                      // of list-based option BW[2]=4. Not the same lists are valid
	{			
				s1:					
				s1=alea_integer(0,param.S-1,randCase); 

				if(s1==s) goto s1;	
				// *** WARNING, may be infinite
	}
			}

			// Gravity centre Gr
						w1=1; w2=1; w3=1; 
			switch(param.BW[1])
				{
				default:if(g==s) w3=0; break; // Pure standard
				case 1: break; // Will use a specific method (see below)
				case 2: break;
				case 3: if(g==s) {w1=0; w3=0;} break;
				case 4: // TEST
				w1=R.SW.P[g].f; w2=R.SW.P[s].f;
				if(g!=s)
				{
	 					w3=R.SW.X[s].f; 
				}
				else
				{	
				w3=0;
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

			zz=w1+w2+w3;
			w1=w1/zz; w2=w2/zz; w3=w3/zz;

			for (d = 0; d < pb.SS.D; d++) 
			{
				Gr.x[d]=w1*R.SW.X[s].x[d] + w2*(R.SW.X[s].x[d] + param.c*PX.v[d]);
	
				if(g!=s) 
				{
					Gr.x[d]=Gr.x[d]+ w3*(R.SW.X[s].x[d] +param.c*GX.v[d]) ;	
				}
				else // Here P=G
				switch(param.BW[1])
				{
					default: break; //0. If "pure standard", do nothing.
					case 1: //G is a random informant (see above). 
						// However, it is used only for high dimension
						// and not even always (probabilistic choice)
						if(alea(0,1,randCase) <(double)param.S/(pb.SS.D*param.K))	// Approximation
						{ // Use P twice
							Gr.x[d]=Gr.x[d]+w3*(R.SW.X[s].x[d] + param.c*PX.v[d]);  
						}
						else
						{	// Random informant
							Gr.x[d]=Gr.x[d]+ w3*(R.SW.X[s].x[d] 
							    +param.c*(R.SW.P[s1].x[d] - R.SW.X[s].x[d]));
						}
						break;
						case 2: // More conservative
						Gr.x[d]=R.SW.P[g].x[d];
						break;
						
				}
				V1.v[d]= Gr.x[d]-R.SW.X[s].x[d]; // Vector X-G
			}

			// Random point around
			switch(param.BW[1])
				{
				default:
			rad=distanceL(R.SW.X[s],Gr,2);
			break;
			case 3:
			rad=(param.c-1)*distanceL(R.SW.X[s],R.SW.P[s],2);
			break;
			}

			V2=alea_sphere(pb.SS.D,rad,param.distrib, param.mean,param.sigma,randCase); 

			// New "velocity"
			for (d = 0; d < pb.SS.D; d++) 
			{
				R.SW.V[s].v[d]=R.SW.V[s].v[d]+V1.v[d] + V2.v[d]; // New "velocity"
			}

			// New position

			for (d = 0; d < pb.SS.D; d++) 
			{	
				R.SW.X[s].x[d] = R.SW.X[s].x[d] + R.SW.V[s].v[d];			
			}

			if (R.nEval >= pb.evalMax)  goto end;

			// --------------------------

			if(pb.SS.quantisation==1)
				R.SW.X[s] = quantis (R.SW.X[s], pb.SS);

			// Confinement			
			out=0;
			switch(param.confin)
			{
				default: 
					for (d = 0; d < pb.SS.D; d++)
				{	
					if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}

					if (R.SW.X[s].x[d] < xMin)
					{	
						R.SW.X[s].x[d] = xMin;
						R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
						out=1;
					}
					else
					{
						if (R.SW.X[s].x[d] > xMax)
						{			
							R.SW.X[s].x[d] = xMax;
							R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
							out=1;
						}
					}				
				}	
					break;
					
				case 1: // No confinement and no evaluation if outside ("let if fly")
					for (d = 0; d < pb.SS.D; d++)
				{	
					if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}

					if (R.SW.X[s].x[d] < xMin || R.SW.X[s].x[d] > xMax) out=1;		
				}	
				break;
			}

			if(pb.SS.quantisation==1 && out>0)
			{
				R.SW.X[s] = quantis (R.SW.X[s], pb.SS);
			}		
			// If the position is inside
			if(param.confin==0 || (param.confin==1 && out<zero))
			{
				// Evaluation
				R.SW.X[s].f =perf(R.SW.X[s],pb.function, pb.SS,pb.objective);
				R.nEval = R.nEval + 1;				
				// ... update the best previous position		
				if (R.SW.X[s].f < R.SW.P[s].f)	// Improvement
				{															
						R.SW.P[s] = R.SW.X[s];

						// ... update the best of the bests
						if (R.SW.P[s].f < R.SW.P[R.SW.best].f)
						{		
							R.SW.best = s;			
						}			
				}		
			}
			
			if(param.trace>0)
			{
			// Keep trace of every position, for further tests
			fprintf(f_trace,"%i %f ",s,R.SW.X[s].f);
				for (d = 0; d < pb.SS.D; d++)
				{
				  fprintf(f_trace,"%f ",R.SW.X[s].x[d]);
				}
			fprintf(f_trace,"\n");
			}	
		}			// End of "for (s0=0 ...  "	

		// Check if finished
		error = R.SW.P[R.SW.best].f;

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

			if (error > pb.epsilon && R.nEval < pb.evalMax)
		{
			noStop = 0;	// Won't stop
		}
			else
		{
			noStop = 1;	// Will stop
		}


	} // End of "while nostop ...

	// printf( "\n and the winner is ... %i", R.SW.best );			
	R.error = error;
	return R;  
}
// ===========================================================
struct position quantis (struct position x, struct SS SS ) 
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
			quantx.x[d] = qd * floor (0.5 + x.x[d] / qd);	    
		}	
	}
	return quantx;    
}


