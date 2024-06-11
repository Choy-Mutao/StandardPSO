#include "sub_programs.h"
#include "kiss_rng.h"

#include "stdlib.h"
#include "math.h"

// ======================随机生成一个 double======================
//double alea(double a, double b, int randOption)
//{				
//	// random number (uniform distribution) in  [a b]
//	// randOption is a global parameter
//	double r;
//	if (randOption == 0)
//		r = a + (double)rand_kiss() * (b - a) / RAND_MAX_KISS;
//	else
//		r = a + (double)rand() * (b - a) / RAND_MAX;
//	//printf("\nalea 944 r %f",r);
//	return r;
//}

double alea(double a, double b, int randOption)
{
	// random number (uniform distribution) in  [a b]
	// randOption is a global parameter
	double r;
	if (randOption == 0)
		r = a + (double)rand_kiss() * (b - a) / RAND_MAX_KISS;
	else
		r = a + (double)rand() * (b - a) / RAND_MAX;
	//printf("\nalea 944 r %f",r);
	return r;
}

// ======================随机生成一个 int======================
extern int alea_integer(int a, int b, int randOption)
{				// Integer random number in [a b]
	int ir;
	double r;

	r = alea(0, 1, randOption);
	ir = (int)(a + r * (b + 1 - a));

	if (ir > b)	ir = b;

	return ir;
}

// ======================伪高斯分布的随机数========================
extern double alea_normal(double mean, double std_dev, int randOption)
{
	/*
	 Use the polar form of the Box-Muller transformation to obtain a pseudo
	 random number from a Gaussian distribution
	 */
	double x1, x2, w, y1;
	// double y2;

	do
	{
		x1 = 2.0 * alea(0, 1, randOption) - 1.0;
		x2 = 2.0 * alea(0, 1, randOption) - 1.0;
		w = x1 * x1 + x2 * x2;
	} while (w >= 1.0);

	w = sqrt(-2.0 * log(w) / w);
	y1 = x1 * w;
	// y2 = x2 * w;
	if (alea(0, 1, randOption) < 0.5) y1 = -y1;
	y1 = y1 * std_dev + mean;
	return y1;
}


// ==========================生成随机速度========================
extern Velocity aleaVector(int D, double coeff, int randOption)
{
	Velocity V;
	int d;
	int i;
	int K = 2;	// 1 => uniform distribution in a hypercube
				// 2 => "triangle" distribution
	double rnd;

	V.size = D;

	for (d = 0; d < D; d++)
	{
		rnd = 0;
		for (i = 1; i <= K; i++) rnd = rnd + alea(0, 1, randOption);
		V.v[d] = rnd * coeff / K;
	}

	return V;
}

// =======================计算向量的L范数===========================
extern double normL(Velocity v, double L)
{   // L-norm of a vector
	int d;
	double n;

	n = 0;

	for (d = 0; d < v.size; d++)
		n = n + pow(fabs(v.v[d]), L);

	n = pow(n, 1 / L);
	return n;
}

// =========================计算两点的L范数=========================
extern double distanceL(Position x1, Position x2, double L)
{  
	// Distance between two positions
	// L = 2 => Euclidean	
	int d;
	double n;

	n = 0;

	for (d = 0; d < x1.size; d++)
		n = n + pow(fabs(x1.x[d] - x2.x[d]), L);

	n = pow(n, 1 / L);
	return n;
}

//===========================M1 * M2===============================
extern Matrix matrixProduct(Matrix M1, Matrix M2)
{
	// Two square matrices of same size
	Matrix Product;
	int D;
	int i, j, k;
	double sum;
	D = M1.size;
	for (i = 0; i < D; i++)
	{
		for (j = 0; j < D; j++)
		{
			sum = 0;
			for (k = 0; k < D; k++)
			{
				sum = sum + M1.v[i][k] * M2.v[k][j];
			}
			Product.v[i][j] = sum;
		}
	}
	Product.size = D;
	return Product;
}

////======================求转置矩阵===========================
//extern Matrix matrixRotation(Velocity V)
//{
//	/*
//	 Define the matrice of the rotation V' => V
//	 where V'=(1,1,...1)*normV/sqrt(D)  (i.e. norm(V') = norm(V) )
//
//	 */
//	Velocity B = { 0 };
//	int i, j, d, D;
//	double normB, normV, normV2;
//	//Matrix reflex1; // Global variable
//	Matrix reflex2;
//	Matrix rotateV;
//	double temp;
//
//	D = V.size;
//	//normV = normL(V, 2); normV2 = normV * normV;
//	normV = normL(V, 2);
//	reflex2.size = D;
//
//	// Reflection relatively to the vector V'=(1,1, ...1)/sqrt(D)	
//	// norm(V')=1
//	// Has been computed just once  (global Matrix reflex1)	
//
//	//Define the "bisectrix" B of (V',V) as an unit vector
//	B.size = D;
//
//	temp = normV / sqrt(D);
//
//	for (d = 0; d < D; d++)
//	{
//		B.v[d] = V.v[d] + temp;
//	}
//	normB = normL(B, 2);
//
//	if (normB > 0)
//	{
//		for (d = 0; d < D; d++)
//		{
//			B.v[d] = B.v[d] / normB;
//		}
//	}
//
//	// Reflection relatively to B
//	for (i = 0; i < D; i++)
//	{
//		for (j = 0; j < D; j++)
//		{
//			reflex2.v[i][j] = -2 * B.v[i] * B.v[j];
//		}
//	}
//
//	for (d = 0; d < D; d++)
//	{
//		reflex2.v[d][d] = 1 + reflex2.v[d][d];
//	}
//
//	// Multiply the two reflections
//	// => rotation				
//	rotateV = matrixProduct(reflex2, reflex1);
//	return rotateV;
//}

//==========================================================
extern Velocity matrixVectProduct(Matrix M, Velocity V)
{
	Velocity Vp;
	int d, j;
	int Dim;
	double sum;
	Dim = V.size;
	for (d = 0; d < Dim; d++)
	{
		sum = 0;
		for (j = 0; j < Dim; j++)
		{
			sum = sum + M.v[d][j] * V.v[j];
		}
		Vp.v[d] = sum;
	}
	Vp.size = Dim;
	return Vp;
}

// =========================取符号==============================
extern int sign(double x)
{
	if (x == 0)	return 0;
	if (x < 0)	return -1;
	return 1;
}

// ===========================================================
extern Position quantis(Position x, SS SS)
{
	/*
	 Quantisation of a Position
	 Only values like x+k*q (k integer) are admissible
	 */
	int d;
	double qd;
	Position quantx;

	quantx = x;
	for (d = 0; d < x.size; d++)
	{
		qd = SS.q.q[d];

		if (qd > 0)	// Note that qd can't be < 0
		{
			//qd = qd * (SS.max[d] - SS.min[d]) / 2;	      
			quantx.x[d] = qd * floor(0.5 + x.x[d] / qd);
		}
	}
	return quantx;
}