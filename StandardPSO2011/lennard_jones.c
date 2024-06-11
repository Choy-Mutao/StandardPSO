double lennard_jones (struct position x) 
{
	/*
		This is for black-box optimisation. Therefore, we are not supposed to know
		that there are some symmetries. That is why the dimension of the problem is
		3*nb_of_atoms, as it could be 3*nb_of_atoms-6
	*/
	int d;
	int dim=3;
	double dist;
	double f;
	int i,j;
	int nPoints=x.size/dim;
	int alpha=6;
	struct position x1, x2;
	double zz;

	x1.size=dim; x2.size=dim;

	f=0;
	for(i=0;i<nPoints-1;i++)
	{
		for(d=0;d<dim;d++)  x1.x[d]=x.x[3*i+d];
		for(j=i+1;j<nPoints;j++)
		{
			for(d=0;d<dim;d++)  x2.x[d]=x.x[3*j+d];
			dist=distanceL(x1,x2,2);
			if(dist<zero) zz=infinity;
			else zz=pow(dist,-alpha); 
			f=f+zz*(zz-1);
		}
	}
	f=4*f;
	return f;
}
