// Repulsion

		min=0;max=1; // WARNING. Hard coded. Must be the same as in problemDef
	
		n=(int)xs.size/Dim; // Number of points	
		             
	u=0.5*(n-1); //1./100; // For the repulsive action of the boundaries
		          // The smaller u, the nearer to the boundaries can be some points


		f=0; // Total repulsive force to minimise
		for (i=0;i<n;i++) // For each charged point ...
		{
		// ... add the repulsive force between the point and the boundaries
			for(d=0;d<Dim;d++)
			{
			f=f+u*(1./pow(xs.x[d+Dim*i]-min,2) + 1./pow(max-xs.x[d+Dim*i],2));
			}
			
			// ... add the repulsive force between the charged points
			for(j=0;j<n;j++)
			{
				if(j==i) continue;
				
				// Distance^2 between point i and point j
				z=0;
				for(d=0;d<Dim;d++) 
				  z=z+pow(xs.x[d+Dim*i]-xs.x[d+Dim*j] ,2);

				// Repulsion
				f=f+1./z;
			}
		}
		
		

