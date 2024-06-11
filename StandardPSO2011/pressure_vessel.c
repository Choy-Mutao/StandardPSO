
		case 25: //  Pressure vessel (penalty method)		
			pb.SS.D=4;

		pb.SS.min[0] = 1.125; 
		pb.SS.max[0] = 12.5;
		pb.SS.q.q[0] = 0.0625;
		
		pb.SS.min[1] = 0.625; 
		pb.SS.max[1] = 12.5; 
		pb.SS.q.q[1] =0.0625;

		pb.SS.min[2] = 0;pb.SS.max[2] = 240; 
		pb.SS.q.q[2] = 0;
		pb.SS.min[3] = 0; pb.SS.max[3] = 240; 
		pb.SS.q.q[3] = 0;	 	

pb.objective = 7197.72893; 
pb.SS.quantisation=1;
pb.SS.normalise=1;	 // Very important here.
		pb.evalMax =50000 ; 
		pb.epsilon =0; // 0.00001; 
		break;