
			// Ref New Optim. Tech. in Eng. p 638
			// D = 4  		
			        
		x1=xs.x[0]; // [1.1,12.5] granularity 0.0625
		x2=xs.x[1];// [0.6,12.5] granularity 0.0625
		x3=xs.x[2]; // [0,240]
		x4=xs.x[3];// [0,240]

		f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1611*x1*x1*x4 + 19.84*x1*x1*x3;
//		f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1161*x1*x1*x4 + 19.84*x1*x1*x3;
  //  	f=0.6224*x1*x3*x4 + 1.7781*x2*x3*x3 + 3.1661*x1*x1*x4 + 19.84*x1*x1*x3;

		ff=constraint(xs,function);
	
			if (ff.f[1]>0) {c= 1+pow(10,10)*ff.f[1];  f=f*c*c; }
			if (ff.f[2]>0) {c=1+ff.f[2]; f=f*c*c;  }
			if (ff.f[3]>0) {c=1+ff.f[3]; f=f*c*c;  }
