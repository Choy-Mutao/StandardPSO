static int compareDoubles (void const *a, void const *b)
{
      // Increasing order
       const double *da = (const double *) a;
       const double *db = (const double *) b;
     
       return (*da > *db) -(*da < *db)  ;
}
// ===========================================================
double distanceL (struct position x1, struct position x2,double L) 
{  // Distance between two positions
	// L = 2 => Euclidean	
	int d;     
	double n;

	n = 0;

	for (d = 0; d < x1.size; d++)
		n = n + pow (fabs(x1.x[d] - x2.x[d]), L);

	n = pow (n, 1/L);
	return n;    
}

// ===========================================================
double Gamma(double u)
{
// Only for two particular cases:
// u = integer
// u = integer + 0.5
	int cas;
	double G;
	int i;
	int k;
	
	k=(int)(u+0.5);
	if(fabs(k-(int)u)<0.5) cas=1; else cas=2;
	
	switch(cas)
	{
		case 1: // u integer
		if(k==1) return 1;
			G=1;
			for(i=2;i<k;i++) G=G*i;
			return G;
		
		case 2: // u = k +1/2. We use the duplication formula
		k=k-1;  	
		G=sqrt(pi)*pow(2,1-2*k)*Gamma(2*k)/Gamma(k);
		return G;
		default:
		ERROR("tools 79. In Gamma, u %f is neither inter nor half-integer");
	}
}

// ===========================================================
double max(double a, double b)
{
	if(a>b) return a;
	return b;
}
// ===========================================================
double min(double a, double b)
{
	if(a<b) return a;
	return b;
}
// ===========================================================
int sign (double x) 
{     
	if (x == 0)	return 0;
	if (x < 0)	return -1;    
	return 1;   
}

