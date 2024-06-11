#ifndef SUB_PROGRAMS_H
#define SUB_PROGRAMS_H

#include "pso_structures.h"
#define	pi acos((long double)-1)

// Sub-programs
extern int sign(double);

double alea(double, double, int);
extern int alea_integer(int, int, int);
extern double alea_normal(double, double, int);
extern Velocity aleaVector(int, double, int);

extern double distanceL(Position, Position, double);
extern Matrix matrixProduct(Matrix, Matrix);
extern Matrix matrixRotation(Velocity);
extern Velocity matrixVectProduct(Matrix, Velocity);
extern double normL(Velocity, double);
extern Position quantis(Position, SS);
#endif // !SUB_PROGRAMS_H
