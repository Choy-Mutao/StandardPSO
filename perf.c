#include "math.h"
#include "pso_structures.h"

#define	pi acos((long double)-1)

int sign(double x)
{
	if (x == 0)return 0;
	if (x < 0) return -1;
	return 1;
}

// ===========================================================
double perf(struct Position xs, int function, struct SS SS, double objective) 
// Evaluate the fitness value for the particle of rank s 
{
	double f = 0;

	switch (function)
	{
	case 0:		
		f = Parabola(xs);
		break;
	case 1:
		f = Griewank(xs);
		break;
	case 2:
		f = Rosenbrock(xs);
		break;
	case 3:
		f = Rastrigin(xs);
		break;
	default:
		break;
	}
}

// ============================Function=======================
double Parabola(struct Position xs)			// Parabola(Shpere)
{
	int d = 0;
	double f = 0, xd;

	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f + xd * xd;
	}

	return f;
}

double Griewank(struct Position xs)			// Griewank
{
	int d = 0;
	double f = 0, p = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f + xd * xd;
		p = p * cos(xd / sqrt((double)(d + 1)));
	}
	f = f / 4000 - p + 1;

	return f;
}

double Rosenbrock(struct Position xs)		// Rosenbrock
{
	int d = 0;
	double f = 0;
	double t0, t1, tt;

	t0 = xs.x[0] + 1; // Solution on (0,...0) when offset = 0
	for (d = 1; d < xs.size; d++)
	{
		t1 = xs.x[d] + 1;
		tt = 1 - t0;
		f += tt * tt;
		tt = t1 - t0 * t0;
		f += 100 * tt * tt;
		t0 = t1;
	}

	return f;
}

double Rastrigin(struct Position xs)		// Rastrigin
{
	int d = 0, k = 10;
	double f = 0, xd;

	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f + xd * xd - k * cos(2 * xd * pi);
	}
	f = f + xs.size * k;

	return f;
}

double Tripod_2D(struct Position xs)		// 2D Tripod function
											// Note that there is a big discontinuity right 
											// on the solution point
{
	double x1, x2;
	double s11, s12, s21, s22;
	double f = 0;

	x1 = xs.x[0];
	x2 = xs.x[1];

	s11 = (1.0 - sign(x1)) / 2;
	s12 = (1.0 + sign(x1)) / 2;
	s21 = (1.0 - sign(x2)) / 2;
	s22 = (1.0 + sign(x2)) / 2;

	// f = s21 * (fabs(x1) - x2); // Solution on (0,0)
	f = s21 * (fabs(x1) + fabs(x2 + 50)); // Solution on (0, -50)
	f = f + s22 * (s11 * (1 + fabs(x1 + 50) +
		fabs(x2 - 50)) + s12 * (2 +
			fabs(x1 - 50) +
			fabs(x2 - 50)));

	return f;
}

double Ackley(struct Position xs)				// Ackley
{
	double sum1 = 0, sum2 = 0, xd, f = 0;
	int d = 0, DD = xs.size;
	for (d = 0; d < DD; d++)
	{
		xd = xs.x[d];
		sum1 = sum1 + xd * xd;
		sum2 = sum2 + cos(2 * pi * xd);
	}

	f = -20 * exp(-0.2 * sqrt(sum1 / DD)) - exp(sum2 / DD) + 20 + exp(1);

	return f;
}

double Schwefel(struct Position xs)				// Schwefel
{
	int d = 0;
	double f = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f - xd * sin(sqrt(fabs(xd)));
	}

	return f;
}

double Schwefel_1_2(struct Position xs)			// Schwefel 1.2
{
	int d, k;
	double f = 0, xd, sum1;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		sum1 = 0;
		for (k = 0; k <= d; k++) sum1 = sum1 + xd;
		f = f + sum1 * sum1;
	}

	return f;
}

double Schwefel_2_22(struct Position xs)		// Schwefel 2.22
{
	int d;
	double sum1 = 0, sum2 = 1, xd, f = 0;
	for (d = 0; d < xs.size; d++)
	{
		xd = fabs(xs.x[d]);
		sum1 = sum1 + xd;
		sum2 = sum2 * xd;
	}
	f = sum1 + sum2;

	return f;
}

double Neumaier_3(struct Position xs)			// Neumaier 3
{
	int d;
	double sum1 = 0, sum2 = 1, xd, f;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d] - 1;
		sum1 = sum1 + xd * xd;
	}
	for (d = 1; d < xs.size; d++)
	{
		sum2 = sum2 + xs.x[d] * xs.x[d - 1];
	}

	f = sum1 + sum2;
	return f;
}

double G3(struct Position xs)					// G3 (constrained) min =0 on (1/sqrt(D), ...)
{
	int d;
	double f = 1, sum1 = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f * xd;
		sum1 = sum1 + xd * xd;
	}
	f = fabs(1 - pow(xs.size, xs.size / 2) * f) + xs.size * fabs(sum1 - 1);
	return f;
}

double Schwefel(struct Position xs)					// Schwefel
{
	int d;
	double f = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f - xd * sin(sqrt(fabs(xd)));
	}
	return f;
}

double Goldstein_Price2D(struct Position xs)					// 2D Goldstein-Price function
{
	double x1 = xs.x[0], x2 = xs.x[1];

	double f = (1 + pow(x1 + x2 + 1, 2) * (19 - 14 * x1 + 3 * x1 * x1 - 14 * x2 + 6 * x1 * x2 + 3 * x2 * x2))
		* (30 + pow(2 * x1 - 3 * x2, 2) *
			(18 - 32 * x1 + 12 * x1 * x1 + 48 * x2 - 36 * x1 * x2 + 27 * x2 * x2));
	return f;
}

double Schaffer_F6(struct Position xs)					// Schaffer F6
{
	double x1 = xs.x[0], x2 = xs.x[1];
	double f = 0.5 + (pow(sin(sqrt(x1 * x1 + x2 * x2)), 2) - 0.5) / pow(1.0 + 0.001 * (x1 * x1 + x2 * x2), 2);
	return f;
}

double Schwefel_2_21(struct Position xs)				// Schwefel 2.21
{
	int d;
	double f = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = fabs(xs.x[d]);
		if (xd > f) f = xd;
	}
	return f;
}
double Lennard_Jones(struct Position x)					// Lennard_Jones
{
	/*
		This is for black-box optimisation. Therefore, we are not supposed to know
		that there are some symmetries. That is why the dimension of the Problem is
		3*nb_of_atoms, as it could be 3*nb_of_atoms-6
	*/
	int d;
	int dim = 3;
	double dist;
	double f;
	int i, j;
	int nPoints = x.size / dim;
	struct Position x1, x2;
	double zz;

	x1.size = dim; x2.size = dim;

	f = 0;
	for (i = 0; i < nPoints - 1; i++)
	{
		for (d = 0; d < dim; d++)  x1.x[d] = x.x[3 * i + d];
		for (j = i + 1; j < nPoints; j++)
		{
			for (d = 0; d < dim; d++)  x2.x[d] = x.x[3 * j + d];

			dist = distanceL(x1, x2, 2);
			zz = pow(dist, -6);
			f = f + zz * (zz - 1);
		}
	}
	f = 4 * f;
	return f;
}

double Gear_Train(struct Position xs)					// Gear train
{
	double f = pow(1. / 6.931 - xs.x[0] * xs.x[1] / (xs.x[2] * xs.x[3]), 2);
	return f;
}

double Sine_sine(struct Position xs)					// Sine-sine function
{
	int d;
	double f = 0, xd;
	for (d = 0; d < xs.size; d++)
	{
		xd = xs.x[d];
		f = f - sin(xd) * pow(sin((d + 1) * xd * xd / pi), 20);
	}
	return f;
}

double Perm(struct Position xs)							// Perm function
{
	int k, d;
	double beta = 10, f = 0;
	double sum1, xd;
	for (k = 0; k < xs.size; k++)
	{
		sum1 = 0;
		for (d = 0; d < xs.size; d++)
		{
			xd = xs.x[d];
			sum1 = sum1 + (pow(d + 1, k) + beta) * (pow(xd / (d + 1), k) - 1);
		}
		sum1 = sum1 * sum1;
		f = f + sum1;
	}
	return f;
}