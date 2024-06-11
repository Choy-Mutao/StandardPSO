		case 100: // Parabola (Sphere) CEC 2005		f=-450;		for (d = 0; d < xs.size; d++) 	  {			xd =xs.x[d]-offset_0[d];
			f = f + xd * xd;  		}		break;
			
		case 102:  // Rosenbrock 
				// Solution point on offset_2 => fitness value = 390				for (d = 0; d < xs.size; d++) 		{			xs.x[d]=xs.x[d]-offset_2[d]+1;		}
		f=390;	
 		for (d = 1; d < xs.size; d++)
	  {     	      	      
			tt=xs.x[d-1]-1;
			f =f+tt*tt;
			
			tt=xs.x[d-1]*xs.x[d-1]-xs.x[d];     			f =f+ 100 * tt*tt;	   		}

//f=log(1+f);		break;
		
		case 103: // Rastrigin 		for (d = 0; d < xs.size; d++) 		{			xs.x[d]=xs.x[d]-offset_3[d];		}		f=-330;
	  k = 10;  			for (d = 0; d < xs.size; d++)    
	  {     			xd = xs.x[d];	      			f =f+ xd * xd - k * cos (2 * pi * xd);	    		}	  		f =f+ xs.size * k;  		break;	

		case 104: // Schwefel (F2)
		case 107: // Schwefel + noise (F4)
		for (d = 0; d < xs.size; d++) 		{			xs.x[d]=xs.x[d]-offset_4[d];		}

    f = 0;
    for (d=0; d<xs.size; d++)
    {
        sum2 = 0.0;
        for (k=0; k<=d; k++)
        {
            sum2 += xs.x[k];
        }
		 f += sum2*sum2;		
    }
if(function==107) f=f*(1+0.4*fabs(alea_normal(0,1,10)));
f=f-450;
		break;

		case 105: // Griewank. WARNING: in the CEC 2005 benchmark it is rotated
	   sum1 = 0.0;
	   sum2 = 1.0;
		 f=-180;
	   for (d=0; d<xs.size; d++)
		 {
		    xd=xs.x[d]-offset_5[d];
				sum1 += xd*xd;
        sum2 *= cos(xd/sqrt(1.0+d));
 		 }
    f =f+ 1.0 + sum1/4000.0 - sum2;
		break;

		case 106: // Ackley 		f=-140;

    sum1 = 0.0;
    sum2 = 0.0;
    for (d=0; d<xs.size; d++)
    {
        xd = xs.x[d]-offset_6[d];
				sum1 += xd*xd;
        sum2 += cos(2.0*pi*xd);
    }
    sum1 = -0.2*sqrt(sum1/xs.size);
    sum2 /= xs.size;
    f = f+ 20.0 + E - 20.0*exp(sum1) - exp(sum2);
		break;