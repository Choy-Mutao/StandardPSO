Standard PSO 2011

Contact for remarks, suggestions etc.: pso@writeme.com

Updates:
2012-01-06  BW[2] is now for randomness options. In particular, one can use
						quasi-random numbers for initialisation (Sobol or Halton sequences)
						Halton sequences seem to be better.
2011-05-25  Added option 3 for BW[1]: when P=G, then just look around X
2011-04-19  Fixed a small bug in the initialisation of the velocity
2011-05-16  Fixed a small bug in Compression spring (third constraint)
2011-02-01  Added Pressure Vessel
2011-01-21	temporary experimental option (may be removed later): variable ring topology
2011-01-21	aleaIndex() as a function in alea.c (was directly coded in main.c)
2011-01-15	A whole family of stable distributions (CMS method)
						May also replace the Box-Muller method for Gaussian distribution
2011-01-12	Penalized function (code 23)
2011-01-07	BW[1]. Specific methods in case P=G
2011-01-05	BW[0]. Random swarm size 
Initial version: 2011-01-01

Note 1: there are a few optional "bells and whistles" (BW) 
			Not really part of the standard, but may improve the performance
			for some problems (or the contrary ...)
			Example: the random swarm size "around" a mean value
			
Note 2: the code is not optimised, in order to be more readable.

  -------------------------------- Contributors 
 Works and comments of the following persons have been taken
 into account while designing this standard.  Sometimes this is for 
 including a feature, and sometimes for leaving one out. 
 
 Auger, Anne
 Blackwell, Tim
 Bratton, Dan
 Clerc, Maurice
 Croussette, Sylvain 
 Dattasharma, Abhi
 Eberhart, Russel
 Hansen, Nikolaus
 Helwig, Sabine
 Keko, Hrvoje
 Kennedy, James 
 Krohling, Renato
 Langdon, William
 Li, Wentao
 Liu, Hongbo 
 Miranda, Vladimiro
 Poli, Riccardo
 Serra, Pablo
 Spears, William
 Stickel, Manfred
 Yue,Shuai
 
 -------------------------------- Motivation
Quite often, researchers claim to compare their version of PSO 
with the "standard one", but the "standard one" itself seems to vary!
Thus, it is important to define a real standard that would stay 
unchanged for at least a few years.

This PSO version does not intend to be the best one on the market
(in particular, there is no adaptation of the swarm size or of the
coefficients). This is simply similar to the original version (1995),
with few improvements based on some recent works.
------
The above comment was for Standard PSO 2007. It is still valid. However, all
"standard" PSOs (2007 and before) are not satisfying for at least one reason: on a 
given problem, the performance depends on the system of coordinates. This
phenomenon is sometimes called "rotation sensitivity". This happens because for each
particle and at each time step, the set of all "next possible positions" (NPP) is
itself coordinate dependent: typically, a combination of two D-rectangles whose 
sides are parallel to the axes, with uniform probability distribution inside
each rectangle. This combination is itself a D-rectangle, but with a 
non uniform distribution (more dense near the centre).

So the obvious idea is to define an NPP domain that is _not_ coordinate dependent.
There are a lot of reasonable simple shapes. Let X be the current position,
P the best previous one, and G the best previous best in the neighbourhood (informants).
We could for example:
- define  a D-sphere centred somewhere between X, P, G (like in some PSO variants)
- define a D-cylinder along X-P (resp. X-G)
- define a D-Ice_cream_cone (base on X, and the ball around P (resp. G))
- etc.

Experiments suggest that a simple hypersphere H(Gr,||Gr-X||)
is a good compromise, where Gr is the "centre of gravity" of X, P', and G', where 
P' is X+c*(P-X) and G'=X+c*(G-X) with c>1.

******
This is the core of this standard. It combines the same information that
the canonical PSO (i.e. X, P, and G), with the same parameters (w and c),
but in a way that is both simple and not coordinate dependent.
******

However, note that perfect non-sensitivity to rotation is impossible (except for 
very particular problems) whenever the search space itself is not a sphere.


------------------------------------ From SPSO-07 to SPSO-2011
 - normalisation of the search space if possible, which is then [0,normalise]^D
	Reason : for some problems, the range of values is not the same for all
 						variables (the search space is not a hypercube).
 						In such a case, using spheres may lead to bad performance.
 						In the definition of the problem, it is then better to set 
 						pb.SS.normalise to a positive value (typically 1). 
 
- new definition of the next possible position domain : hypersphere H(Gr,||Gr-X||)
Reason: non sensitivity to rotation (as much as possible)

- also, a distribution option has been added (Gaussian distribution instead of the
	uniform one)
	Reason: global convergence property. It does not always mean
					"better performance", but quite often, it indeed works better.

- deterministic formula for swarm size has been removed. Now the user gives  the
value (default 40), or a mean value and for each run the swarm size is randomly 
chosen "around" this value
  Reason: the size given by the formula in SPSO 2007 (depending on the dimension) 
  				often used to come out as far from optimal.

- velocity initialisation is slightly modified
 	Reason: so that no particle leaves the search space at the very first step. 
 	
Note: There do exist far better initialisation schemes than the random one
used for positions in SPSO 2007 (it has not been modified here), but position 
initialisation is not really an inherent part of a search strategy.
Actually, for fair comparison of two iterative algorithms, initial positions 
should be the same. 
         	
- a (non) confinement option has been added (the old historical "let them fly" method)
Reason : All other methods are more or less arbitrary. There is 
				a drawback, though. When the search space is discrete (at least along some dimensions)
				the convergence may be extremely slow when using some "bells and whistles".
				See problem 11 (Network), which is partially binary.

 --------------------------------- Metaphors
swarm: A team of communicating people (particles)
At each time step
    Each particle chooses a few informants , selects the best
	one from this set, and takes into account the information given by
	the chosen particle (the informant).
	If it finds no particle better than itself, it may choose the informant at random.
	
----------------------------------- Parameters/Options
All parameters have been hard coded with default values, but you can of course
modify them.

S := swarm size
K := maximum number of particles _informed_ by a given one
w := first cognitive/confidence coefficient
c := second cognitive/confidence coefficient

 ----------------------------------- Equations
For each particle 
Equation 1:	v(t+1) = w*v(t) + V(t)
Equation 2:	x(t+1) = x(t) + v(t+1)
where
V(t) := Gr(t)- x(t) + r*alea_sphere(Gr,||Gr-X||)
Gr(t)	:= (x(t)+p'(t)+l'(t))/3

v(t) 	:= velocity at time t
x(t) 	:= position at time t

p'(t)	:= x(t)+c*(p(t)-x(t))
l'(t)	:= x(t)+c*(l(t)-x(t))

p(t) 	:= previous best position of the particle. Best position of the particle 
				that has been found so far
l(t) 	:= local best. Best position amongst the previous best positions of the 
				informants of the particle 
w			:= inertia weight
c			:= acceleration coefficient


Note 1: It is tempting to define a weighted sum
Gr(t)	:= w1*x(t)+w2*p(t)+w3*l(t)

In practice, a lot of weighting schemes have been tried. None of them
was really better than the (1/3,1/3,1/3) used here (on average on a large
benchmark. For a specific problem, a particular scheme may indeed be better).

Note 2:
When the particle has no informant better than itself,
it implies p(t) = l(t)
The simplest way is to do nothing special. We then have
Gr(t)	:= (x(t)+2*p'(t))/3 

The standard rule here is
Gr(t)	:= (x(t)+ p'(t))/2

A BW is to simply say that Gr=G.

Another one is to choose a different l(t) at random amongst the other particles
with a probability that depends on S, K, and on the dimension D. 
It usually improves the performance (but not always).
 
 ----------------------------------- Use
 Define the problem (you may add your own one in problemDef() and perf())
 You may modify some parameters.
 Run and enjoy!
================================================================================
A few results over 100 runs (success rate and Mean best).
For reproducible results the seed of the random number generator KISS is set to 
1234567890, and is "warmed up" 10000 times before the first run fo each function.

As 100 runs is not big enough, the success rate is indeed very sensitive to
the RNG and to the seed.

< see the file Results.pdf >

Comments:
1) The performance may be better because the problem is biased (typically the solution point lies on (0,0,...,0))
or because, by chance, the swarm size (which is automatically computed in SPSO 2007 by using a formula depending one the dimension) is a "good" one (for each problem, there exists an optimal swarm size).

3) The performance with the Gaussian distribution is slightly better
than the one with the uniform distribution, except without any confinement 
and some BW.

