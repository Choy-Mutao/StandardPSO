	case 100: // CEC 2005 F1
			pb.SS.D =30;//30; 
			for (d = 0; d < pb.SS.D; d++)
			{   
				pb.SS.min[d] =- -100;
				pb.SS.max[d] =  100;
				pb.SS.q.q[d] = 0;	
  
			}
			pb.evalMax =pb.SS.D*10000;
			pb.epsilon = 0.000001;	//Acceptable error
			pb.objective =-450;       // Objective value
			
				for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
			break;
			
			case 102:		// Rosenbrock. CEC 2005 F6
		pb.SS.D =30;
				
					// Boundaries
		for (d = 0; d < pb.SS.D; d++)
		{	
			pb.SS.min[d] = -100; pb.SS.max[d] = 100; 			
			pb.SS.q.q[d] = 0;	
      
		}
				
		pb.evalMax = pb.SS.D*10000;
		pb.epsilon = 0.01;	//0.01 Acceptable error
		pb.objective =390; 

	for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
		break;
		
			case 103:// CEC 2005 F9, Rastrigin
		pb.SS.D = 30;	 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] =-5.12; 
			pb.SS.max[d] =5.12; 	  
			pb.SS.q.q[d] = 0;	

		 }
		pb.epsilon =0.01; // 0.01;	// Acceptable error
		pb.objective =-330;       // Objective value
		pb.evalMax=pb.SS.D*10000;

	for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
		break;
		
			case 104:// CEC 2005 F2  Schwefel
		pb.SS.D = 10;	 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] =-100; 
			pb.SS.max[d] =100; 	  
			pb.SS.q.q[d] = 0;	

		 }
		pb.epsilon = 0.00001;	// Acceptable error
		pb.objective =-450;       // Objective value
		pb.evalMax=pb.SS.D*10000;
		
		pb.epsilon=0.0001	; pb.evalMax=100000;	


	for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
		break;
				
		
		case 105:// CEC 2005 F7  Griewank (NON rotated)
		pb.SS.D = 10;	 // 10 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] =-600; 
			pb.SS.max[d] =600; 	  
			pb.SS.q.q[d] = 0;	

		 }
		pb.epsilon = 0.01;	//Acceptable error
		pb.objective =-180;       // Objective value
		pb.evalMax=pb.SS.D*10000;
		
			for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
		break;
				
				
		case 106:// CEC 2005 F8 Ackley (NON rotated)
		
			pb.SS.D =30; // 10;	 
		for (d = 0; d < pb.SS.D; d++)
		{
			pb.SS.min[d] =-32; 
			pb.SS.max[d] =32; 	  
			pb.SS.q.q[d] = 0;	
		 }
		pb.epsilon = 0.0001;	// Acceptable error
		pb.objective =-140;       // Objective value
		pb.evalMax=pb.SS.D*10000;
	for (d = 0; d < pb.SS.D; d++)
		{  
			pb.SS.maxInit[d]=pb.SS.max[d];
			pb.SS.minInit[d]=pb.SS.min[d];
		}
		break;
				