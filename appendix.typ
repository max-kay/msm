#set page(numbering: "1")
#set heading(numbering: "A.1.1.")
#set math.equation(numbering: "(1)")
#let comment(c) = {
  text(fill: rgb("#FF3333"))[#c]
}

#text(weight: "bold", size: 24pt)[Appendix]

This document derives and summerizes all math used for our simulation.
Section 1 contains all necessary information to implement the model.
Section 2 contains our optimization strategies.
Section 3 contains the derivations of the math in Section 1.

#outline()

#pagebreak()

= Summary

== Parameters
This section describes all parameters set at the start of the simulation.
#align(
  center,
  table(
    columns: 2,
    align: center + horizon,
    inset: 10pt,
  )[*Description*][ *Value*
  ][particle number][ $N=100$
  ][volume fraction][ $f=1%$
  ][magnetic moment density][ $M = 380 "kA"/ "m"$
  ][relative permitivity in matrix][ $epsilon_m = 2$
  ][relative permitivity in BHF][ $epsilon_p = 10$
  ][external magnetic field][ $norm(bold(H)_0) = 5M$
  ][external electric field][ $norm(bold(E)_0) = 100 "MV"/"m"$
  ][dynamic viscosity][ $eta = 3.5 "Pa" dot "s"$
  ][long semiaxes][ $a=b = 2.5 mu"m"$
  ][short semiaxis][ $c = 2/7 a$
  ][factor for repulsion force][ $beta=40$\ to be tuned after implementation
  ][velocity factor][ $k=1/3$\ to be tuned after implementation
  ][simulation time][ $T~1 "s"$\ to be tuned after implementation
  ],
)

#pagebreak()


== Calculated Constants
This section describes all values used in the simulation which can be calculated in advance.
#align(
  center,
  table(
    columns: 2,
    align: center + horizon,
    inset: 10pt,
  )[*Description*][ *Formula*
  ][Volume of a particle][ $V_p = 4/3 pi a b c$
  ][side length of the representative volume element (RVE)][ $L = root(3, (N V_p)/ f)$
  ][radius used for drag forces and torques, repulsive force and the calculation of $Delta t$ ][ $R = root(3, a^2 c)$
  ][translational drag coefficient][ $D_s = 6 pi eta R$
  ][rotational drag coefficient][ $C_s = 8 pi eta R^3$// Electromechanics of Particles p 88
  ][effective electric suceptibilty along the long semiaxes][
    $chi_x = (epsilon_p - epsilon_m) / (1 + (epsilon_p - epsilon_m)/epsilon_m L_x)$
  ][effective electric suceptibilty along the short semiaxis][
    $chi_z = (epsilon_p - epsilon_m) / (1 + (epsilon_p - epsilon_m)/epsilon_m (1 -2L_x))$
  ][shape factor for the effective electric suceptibility][
    $L_x = (a^2c)/(2 (a^2 - c^2)) [(pi/2)/sqrt(a^2 - c^2)- c/a^2]$
  ][the norm of the dipole moment vector of the particles][
    $m=norm(bold(m)_i)=V_p M$
  ],
)

#pagebreak()

== Positions and Boundary Conditions
Each particle $i$ has a position $bold(x)_i$ and an orienation $bold(d)_i$ with $norm(bold(d)_i)=1$

The vector $bold(r)_(j i)$ is the vector which points from the position of particle $j$ to the position of particle $i$
$
       bold(r)_(j i) & = bold(x)_i - bold(x)_j \
             r_(j i) & = norm(bold(r)_(j i)) \
  hat(bold(r))_(j i) & = bold(r)_(j i)/ r_(j i)
$

We used nearest neighbour boundary condition.
This means that each position $bold(x)$ is equivalent to position $tilde(bold(x)) + L bold(v)$ with $bold(v) in ZZ^3$.

Where ever the translation vector $bold(r_(i j))$ is used in equation we refer to the shortest vector between equivalent positions.
Such that every coefficient in $bold(r_(i j))$ is $in [-L/2, L/2)$.


== The Fields

The magnetic field:
$
  bold(H)_i = bold(H)_0
  + sum_(j!=i) m/(4 pi) 1/r_(j i)^3[
    3 (bold(d)_j dot hat(bold(r))_(j i))hat(bold(r))_(j i)
    - bold(d)_j
  ]
$

The electric field:
$
  bold(E)_i = bold(E)_0
  + sum_(j!=i) 1/(4 pi epsilon_0 epsilon_m) 1/r_(j i)^3[
    3 (bold(p)_j dot hat(bold(r))_(j i))hat(bold(r))_(j i)
    - bold(p)_j
  ]
$

With the polarization $p$
$ bold(p)_i = V_p epsilon_0 [chi_x bold(E)_i + (chi_z - chi_x) (bold(d)_i dot bold(E)_i) bold(d)_i] $

Since the electric field and the dipole moments depend on each other, they need to be approximated by the following iterative process.

Start by calculating the electric field by the assumption
that each particle has the same polarization as in the last time step, denoted $tilde(bold(p))_i$.

$
  bold(E)_i^0 = bold(E)_0
  + sum_(j!=i)1/(4 pi epsilon_0 epsilon_m) 1/r_(j i)^3[
    3 (tilde(bold(p))_j dot hat(bold(r))_(j i))hat(bold(r))_(j i)
    - tilde(bold(p))_j
  ]
$


Then iteratively updating the electric field and electric dipole moments by:
$ bold(p)_i^n = V_p epsilon_0 [chi_x bold(E)_i^n + (chi_z - chi_x) (bold(d)_i dot bold(E)_i^n) bold(d)_i] $
$
  bold(E)_i^(n+1) = bold(E)_0
  + sum_(j!=i)1/(4 pi epsilon_0 epsilon_m) 1/r_(j i)^3[
    3 (bold(p)_j^n dot hat(bold(r))_(j i))hat(bold(r))_(j i)
    - bold(p)_j^n
  ]
$

Note that the first time step has to be treated specially since there are no previous steps.
This can be handeled by calculating the electric dipole moments in the initialisation step of the simulation based on $bold(E)_0$.

#pagebreak()

== Translation

The differential equation for the translational motion of the platelets is:

$ dot(bold(x))_i = 1/D_s (bold(F)^H_i + bold(F)^E_i + bold(F)^R_i) $

With the equations for the forces being:
$
  bold(F)^H_i
  = sum_(j!=i)(3 mu_0 m^2)/(4 pi) 1/r_(j i)^4( & [
                                                   (bold(d)_j dot bold(d)_i)
                                                   - 5(hat(bold(r))_(j i) dot bold(d)_j) (hat(bold(r))_(j i) dot bold(d)_i)
                                                 ]hat(bold(r))_(j i) \
                                               & quad +[
                                                   (hat(bold(r))_(j i) dot bold(d)_i)bold(d)_j
                                                   + (hat(bold(r))_(j i) dot bold(d)_j) bold(d)_i
                                                 ]
                                                 )
$

$
  bold(F)^E_i
  = sum_(j!=i)3/(epsilon_0 epsilon_m 2 pi)1/r_(j i)^4( & [
                                                           (bold(p)_j dot bold(p)_i)
                                                           - 5(hat(bold(r))_(j i) dot bold(p)_j) (hat(bold(r))_(j i) dot bold(p)_i)
                                                         ]hat(bold(r))_(j i) \
                                                       & quad+ [
                                                           (hat(bold(r))_(j i) dot bold(p)_i)bold(p)_j
                                                           + (hat(bold(r))_(j i) dot bold(p)_j) bold(p)_i
                                                         ])
$

$ bold(F)^R_i = sum_(j!=i) (3 mu_0 m^2)/(2 pi (2 R)^4) e^(-beta (r_(j i)/(2 R) - 1)) hat(bold(r))_(j i) $

Note that all summation are over antisymmertic terms. ($bold(v)_(j i)=-bold(v)_(i j)$)

== Rotation

The differential equation for the rotational motion of the platelets is:


$ dot(bold(d))_i = 1/C_s (bold(T)^H_i times bold(d)_i + bold(T)^E_i times bold(d)_i) $

With the equations for the torques being:
$ bold(T)^H_i times bold(d)_i = mu_0 m (bold(H)_i - bold(d)_i (bold(H)_i dot bold(d)_i)) $

$
  bold(T)^E_i times bold(d)_i = V_p epsilon_0 (chi_z - chi_x) (bold(d)_i dot bold(E)_i) [
    bold(E)_i - bold(d)_i (bold(E)_i dot bold(d)_i)
  ]
$

== Evolution of the Positions and Directions

The time step $Delta t$ is determined from the greatest velocity:
$ Delta t = (k R)/ max(norm(dot(bold(x))_i)) $

The time evolution then is:

$ bold(d)_i (t + Delta t) = bold(d)_i (t) + Delta t dot(bold(d))_i (t) $
$ bold(x)_i (t + Delta t) = bold(x)_i (t) + Delta t dot(bold(x))_i (t) $


#pagebreak()

= Optimization Strategies

== Parallelization of the Force Calculation

The following system is designed for a system which can efficiently run $M$ threads in parallel.



The force calculation follows the following formula:
$ bold(F)_i = sum^N_(j=0\ j!=i) bold(f)_(j i) $

with $bold(f)_(i j) = -bold(f)_(j i)$ this can also be represented with the following matrix form:

$ bold(F)_i = bold(A) bold(1) $
Where $bold(A)_(i j) = bold(f)_(j i)$ is an antisymmetric matrix and $bold(1) = vec(1, 1, dots.v, 1)$

We can split the calculation of the matrix $bold(A)$ into block matricies:

$
  bold(A) = mat(
    d_(0,0), dot, dot, dot, dot, dots.h, dot;
    b_0, d_(0,1), dot, dot, dot, dots.h, dot;
    b_1, b_(2m-1), d_(1,0), dot, dot, dots.h, dot;
    b_2, b_(2m), b_(4m-3), d_(1,1), dot, dots.h, dot;
    b_3, b_(2m+1), b_(4m-2), b_(6m-6), d_(2,0), dots.h, dot;
    dots.v, dots.v, dots.v, dots.v, dots.v, dots.down, dots.v;
    b_(2m-2), b_(4m-4), b_(6m-7), b_(8m-11), b_(10m-16), dots.h, d_(m,1);
  )
$
This devides the calculations into square blocks of size $N/(2m) times N/(2m)$.
Combining two of the diagonal blocks $d_(i, 0), d_(i, 1)$ into one to balance the load between blocks.

This results in a total number of blocks $M = 2m^2$, where $M$ is less than the number of available blocks.

E.g. $M = 50$ for Euler for which this system was designed.

def for $k$ is $N = k dot 2m$.

We can use a variant of Newton's short cut within each block to calculate each antisymmertic entry
of the original matrix $bold(A)$ only once. Furthermore this optimizes cache hits since each thread
only needs to access memory of two compact blocks.


The first $2m^2 + m$ blocks are the of diagonal blocks.

To generate the ranges of indecies for which each block is repsonsible for we can use the following algorithm:


```python
row_start = []
col_start = []

current_row = k
current_col = 0

while True:
    row_start.append(current_row)
    col_start.append(current_col)
    current_row += k
    if current_row >= N:
        current_col += k
        current_row = current_col + k
    if current_col >= N - k:
        break
```

And for the diagonals:

```python
first_subblock = []
second_subblock = []
for i in range(m):
    first_subblock.append(2*i*k)
    second_subblock.append((2*i+1)*k)
```

Each "normal" block thread will then produce two arrays of length $k$ as output. One the sum of each row and
the negative of the sum of each column.

For the diagonal blocks the calculation also produces two arrays of lenght $k$ the difference being since the block matrix is
antisymetric we use one of the buffers for the first diagonal block and the other for the second diagonal block.

The summation to achieve the total velocity is then a bit more complex:
```python

lower_blocks = ...
upper_blocks = ...
diag_blocks = ...

velocities = [0] * N

for i in range(N):
    j = 0
    while (j < j//k):
        velocity[i] += lower_blocks[j]


```

== Other Calculation

The other calculations do not profit from a complicated division scheme, since none of them involve antisymmetry.
Thus the division into separate parts for calculation can simply be done by deviding the buffers into equal parts.

Since we have $2m^2$ threads running and $N=2k m$ particles,
we need each thread to be able to run $N/(2m^2) = (k m)/m^2 = k/m$ particles.

To make this simpler we choose $k$ such that $m|k$ i.e. $k = m n$ and thus $N = 2 n m^2$

Since we want each of the complicated "Newton shortcut" blocks and all the equal division blocks


== Avoiding data races

The structure of the simulation is the following:
```c
struct Simulation {
    Params parameters;
    Vec3 *h_field;
    Vec3 *e_field;
    Vec3 *e_dipole;
    Vec3 *velocity;
    Vec3 *position;
    Vec3 *direction;
};

void update_h_field();

void update_e_field();
void update_e_dipole();

void update_velocity();

void find_delta_t();

void update_directions(double delta_t);
void update_positions(double delta_t);
```


#pagebreak()

= Derivations

== Derivation of the electric dipole moment $p$
To calculate $bold(p)$ of each particle consider an axis aligned ellipoid i.e. the region
$(bold(x)^"T" "diag"(1/a^2, 1/b^2, 1/c^2) bold(x)) < 1$ in a homogeneous electric field $bold(E)$

#comment[Electromechanics of Particles pp 111]

$ bold(p) = V_p epsilon_0 (epsilon_p - epsilon_m) bold(E)^"in" $
with:
$ bold(E)^"in"_x = E_x/(1 + (epsilon_p - epsilon_m)/epsilon_m L_x) $
and:
$
  L_x = (a b c)/2 integral_0^oo 1/((s + a^2)R_s) "d"s\
  "where" R_s = sqrt((s + a^2)(s + b^2)(s + c^2))
$

Rewriting using a matrix for the effective suceptibility:
$ bold(p) = V_p epsilon_0 bold(chi) bold(E) $
with:
$
  bold(chi) = "diag"(
    (epsilon_p - epsilon_m)/(1 + (epsilon_p - epsilon_m)/epsilon_m L_x),
    (epsilon_p - epsilon_m)/(1 + (epsilon_p - epsilon_m)/epsilon_m L_y),
    (epsilon_p - epsilon_m)/(1 + (epsilon_p - epsilon_m)/epsilon_m L_z)
  )
$

#comment[Electro mechanics of particles pp 117]
When $a=b>c$:
$ L_x = L_y = (a^2c)/(2 (a^2 - c^2)) [(pi/2)/sqrt(a^2 - c^2)- c/a^2] "and" L_z = 1 - 2 L_x $


To describe a ellipsoid of any orientation we us the unit direction vector $bold(d)$
which for the aligned ellipsoid is $bold(d) = vec(0, 0, 1)$.

Using the orthagonal transformation matix $bold(T)$, which transforms vectors from the particle basis into the reference basis
we can write:

$ bold(p) = V_p epsilon_0 bold(T) bold(chi) bold(T)^T bold(E) $


We know:
$ bold(d) = bold(T) vec(0, 0, 1) "i.e." bold(T) = mat(bold(e)^1, bold(e)^2, bold(d)) $
where $bold(e)^1, bold(e)^2$ have a rotational degree of freedom since the particle is rotationally symetric around $bold(d)$

To finalize the derivation, we utilize the property that for a rotationally symmetric particle ($a=b$),
the susceptibility matrix in the particle basis is:
$ bold(chi) = "diag"(chi_x, chi_x, chi_z) $

Using the basis $bold(T) = mat(bold(e)^1, bold(e)^2, bold(d)_i)$, the transformation $bold(T) bold(chi) bold(T)^T$ can be expanded via outer products:
$ bold(T) bold(chi) bold(T)^T = chi_x (bold(e)^1 bold(e)^1^T + bold(e)^2 bold(e)^2^T) + chi_z (bold(d) bold(d)^T) $

Since the basis is orthonormal, the completeness relation states $bold(e)^1 bold(e)^1^T + bold(e)^2 bold(e)^2^T + bold(d) bold(d)^T = II$. Substituting this to eliminate the arbitrary vectors $bold(e)^1$ and $bold(e)^2$:
$
  bold(T) bold(chi) bold(T)^T & = chi_x (II - bold(d) bold(d)^T) + chi_z (bold(d) bold(d)^T) \
                              & = chi_x II + (chi_z - chi_x) bold(d) bold(d)^T
$

The final expression for the dipole moment $bold(p)$ in the reference basis is:
$ bold(p) = V_p epsilon_0 [chi_x bold(E) + (chi_z - chi_x) (bold(d) dot bold(E)) bold(d)] $

Where:
$ chi_x = (epsilon_p - epsilon_m) / (1 + (epsilon_p - epsilon_m)/epsilon_m L_x) $
$ chi_z = (epsilon_p - epsilon_m) / (1 + (epsilon_p - epsilon_m)/epsilon_m (1 - 2 L_x)) $

#pagebreak()

== Translation

=== Terminal Velocity by Drag Force

Staring with Newtons law plug in the eqaution for the drag force $bold(F)^D_i$:
$ m_p dot.double(bold(x))_i = bold(F)^H_i + bold(F)^E_i + bold(F)^R_i + bold(F)^D_i $
$ bold(F)^D_i = -D_s dot(bold(x))_i $
$ m_p dot.double(bold(x))_i = bold(F)^H_i + bold(F)^E_i + bold(F)^R_i - D_s dot(bold(x))_i $

Set $dot.double(bold(x)) = 0$ for the overdamping:
$ D_s dot(bold(x))_i = bold(F)^H_i + bold(F)^E_i + bold(F)^R_i $

Rearranging we get the equation for the velocity of particle $i$:
$ dot(bold(x))_i = 1/D_s (bold(F)^H_i + bold(F)^E_i + bold(F)^R_i) $

=== Magnetic Force

The derivative of the potential $bold(H)_i dot bold(m)_i$ is the magnetic force acting on particle $i$.
$
  bold(F)^H_i & = mu_0 nabla_i (bold(H)_i dot bold(m)_i) \
              & = mu_0 (nabla_i bold(H)_i^T) bold(m)_i \
$

where the term $nabla_i bold(H)_i$ is a Jacobian matrix.

Using the equation for $bold(H)_i$ we get:

$
  bold(F)^H_i = (nabla_i bold(H)_i^T) bold(m)_i & =
                                                  sum_(j!=i) mu_0/(4 pi) nabla_i (1/r_(j i)^3[
                                                      3 (hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T
                                                      - bold(m)_j^T
                                                    ]) bold(m)_i \
                                                & = sum_(j!=i) bold(f)^H_(j i)
$

where $bold(f)^H_(j i)$ is the force exerted on particle $i$ by particle $j$.
#pagebreak()

Using
$
  nabla_i r_(j i) = hat(bold(r))_(j i) quad "and" quad
  nabla_i hat(bold(r))_(j i)^T = 1/r_(j i) (II - hat(bold(r))_(j i) hat(bold(r))_(j i)^T)
$

$
  (4 pi)/mu_0 bold(f)^H_(j i) & = nabla_i (1/r_(j i)^3[
                                    3 (hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T
                                    - bold(m)_j^T
                                  ]) bold(m)_i \
                              & = (nabla_i 1/r_(i j)^3)[
                                  3 (hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T
                                  - bold(m)_j^T
                                ] bold(m)_i \
                              & quad+ 3/r_(j i)^3 (nabla_i
                                  [(hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T]) bold(m)_i \
                              & = - 3/r_(j i)^4 hat(bold(r))_(j i)[
                                  3 (hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T
                                  - bold(m)_j^T] bold(m)_i \
                              & quad+ 3/r_(j i)^3(
                                  (nabla_i hat(bold(r))_(j i)^T)bold(m)_j hat(bold(r))_(j i)^T
                                  +(hat(bold(r))_(j i)^T bold(m)_j)(nabla_i hat(bold(r))_(j i)^T)
                                ) bold(m)_i \
                              & = - 3/r_(j i)^4 hat(bold(r))_(j i)[
                                  3 (hat(bold(r))_(j i)^T bold(m)_j)hat(bold(r))_(j i)^T
                                  - bold(m)_j^T] bold(m)_i \
                              & quad+ 3/r_(j i)^3(
                                  1/r_(j i) (II - hat(bold(r))_(j i) hat(bold(r))_(j i)^T)bold(m)_j hat(bold(r))_(j i)^T
                                  +(hat(bold(r))_(j i)^T bold(m)_j)1/r_(j i) (II - hat(bold(r))_(j i) hat(bold(r))_(j i)^T)
                                ) bold(m)_i \
                              & = 3/r_(j i)^4[
                                  - 3 hat(bold(r))_(j i) hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i)^T
                                  + hat(bold(r))_(j i) bold(m)_j^T
                                  + II bold(m)_j hat(bold(r))_(j i)^T
                                  - hat(bold(r))_(j i) hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i)^T
                                  + hat(bold(r))_(j i)^T bold(m)_j II
                                  - hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i) hat(bold(r))_(j i)^T
                                ]bold(m)_i
$
distributing $bold(m)_i$
$
  (4 pi)/mu_0 bold(f)^H_(j i) = 3/r_(j i)^4[
    & - 3 hat(bold(r))_(j i) hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i)^T bold(m)_i \
    & quad + hat(bold(r))_(j i) bold(m)_j^T bold(m)_i \
    & quad + II bold(m)_j hat(bold(r))_(j i)^T bold(m)_i \
    & quad - hat(bold(r))_(j i) hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i)^T bold(m)_i \
    & quad + hat(bold(r))_(j i)^T bold(m)_j II bold(m)_i \
    & quad - hat(bold(r))_(j i)^T bold(m)_j hat(bold(r))_(j i) hat(bold(r))_(j i)^T bold(m)_i
  ]
$
and using the dot product:
$
  (4 pi)/mu_0 bold(f)^H_(j i) = 3/r_(j i)^4[
    & - 3 hat(bold(r))_(j i) (hat(bold(r))_(j i) dot bold(m)_j) (hat(bold(r))_(j i) dot bold(m)_i) \
    & quad + hat(bold(r))_(j i) (bold(m)_j dot bold(m)_i) \
    & quad + bold(m)_j (hat(bold(r))_(j i) dot bold(m)_i) \
    & quad - hat(bold(r))_(j i) (hat(bold(r))_(j i) dot bold(m)_j) (hat(bold(r))_(j i) dot bold(m)_i) \
    & quad + (hat(bold(r))_(j i) dot bold(m)_j) bold(m)_i \
    & quad - (hat(bold(r))_(j i) dot bold(m)_j) hat(bold(r))_(j i) (hat(bold(r))_(j i) dot bold(m)_i)
  ]
$

collecting terms:
$
  (4 pi)/mu_0 bold(f)^H_(j i) = & 3/r_(j i)^4[
                                    - 3(hat(bold(r))_(j i) dot bold(m)_j) (hat(bold(r))_(j i) dot bold(m)_i)
                                    - (hat(bold(r))_(j i) dot bold(m)_j) (hat(bold(r))_(j i) dot bold(m)_i)
                                    - (hat(bold(r))_(j i) bold(m)_j) (hat(bold(r))_(j i) bold(m)_i)
                                    + (bold(m)_j dot bold(m)_i)
                                  ]hat(bold(r))_(j i) \
                                & quad + 3/r_(j i)^4[
                                    bold(m)_j (hat(bold(r))_(j i) dot bold(m)_i)
                                    + (hat(bold(r))_(j i) dot bold(m)_j) bold(m)_i
                                  ]
$
finally:

$
  bold(f)^H_(j i) = & (3 mu_0)/(4 pi) 1/r_(j i)^4[
                        (bold(m)_j dot bold(m)_i)
                        - 5(hat(bold(r))_(j i) dot bold(m)_j) (hat(bold(r))_(j i) dot bold(m)_i)
                      ]hat(bold(r))_(j i) \
                    & quad+ (3 mu_0)/(4 pi) 1/r_(j i)^4[
                        (hat(bold(r))_(j i) dot bold(m)_i)bold(m)_j
                        + (hat(bold(r))_(j i) dot bold(m)_j) bold(m)_i
                      ]
$
or:
$
  bold(f)^H_(j i) = & (3 m^2 mu_0)/(4 pi) 1/r_(j i)^4[
                        (bold(d)_j dot bold(d)_i)
                        - 5(hat(bold(r))_(j i) dot bold(d)_j) (hat(bold(r))_(j i) dot bold(d)_i)
                      ]hat(bold(r))_(j i) \
                    & quad+ (3 m^2 mu_0)/(4 pi) 1/r_(j i)^4[
                        (hat(bold(r))_(j i) dot bold(d)_i)bold(d)_j
                        + (hat(bold(r))_(j i) dot bold(d)_j) bold(d)_i
                      ]
$

And so:
$
  bold(F)^H_i & = sum_(j!=i)bold(f)^H_(j i) \
$
Note:
$ bold(f)^H_(j i) = - bold(f)^H_(i j) $

#pagebreak()

=== Electric Force
$
  bold(F)^E_i & = nabla_i (bold(E)_i dot bold(p)_i) \
              & = (nabla_i bold(E)_i^T) bold(p)_i + (nabla_i bold(p)_i^T) bold(E)_i \
$

Lets consider the second term first:

$
  (nabla_i bold(p)_i^T) bold(E)_i
  &= V_p epsilon_0 (nabla_i [chi_x bold(E)_i^T + (chi_z - chi_x) (bold(d)_i dot bold(E)_i) bold(d)_i^T]) bold(E)_i\
  &= V_p epsilon_0 [chi_x (nabla_i bold(E)_i^T) + (chi_z - chi_x) (nabla_i bold(E)_i^T) bold(d)_i bold(d)_i^T] bold(E)_i \
  &= V_p epsilon_0 [chi_x (nabla_i bold(E)_i^T) bold(E)_i + (chi_z - chi_x) (nabla_i bold(E)_i^T) bold(d)_i (bold(d)_i^T bold(E)_i)] \
  &= (nabla_i bold(E)_i^T) ( V_p epsilon_0 [chi_x bold(E)_i + (chi_z - chi_x) (bold(d)_i dot bold(E)_i) bold(d)_i] ) \
  &= (nabla_i bold(E)_i^T) bold(p)_i
$

And so
$
  bold(F)^E_i = 2 (nabla_i bold(E)_i^T) bold(p)_i
$

This term in similar to the dipole interaction term for the magnetic dipoles.
$ bold(F)^E_i = sum_(j!=i) bold(f)^E_(j i) $

$
  bold(f)^E_(j i) = & 2/(epsilon_0 epsilon_m 4 pi) 3/r_(j i)^4[
                        (bold(p)_j dot bold(p)_i)
                        - 5(hat(bold(r))_(j i) dot bold(p)_j) (hat(bold(r))_(j i) dot bold(p)_i)
                      ]hat(bold(r))_(j i) \
                    & quad+ 2/(epsilon_0 epsilon_m 4 pi) 3/r_(j i)^4[
                        (hat(bold(r))_(j i) dot bold(p)_i)bold(p)_j
                        + (hat(bold(r))_(j i) dot bold(p)_j) bold(p)_i
                      ]
$

Note the additional factor of 2 and the change of the constants when compared to the calculations for the magnetic force.

=== Repulsive Force<repulsion>

$ bold(F)^R_i = sum_(j!=i) (3 mu_0 m^2)/(2 pi (2 R)^4) e^(-beta (r_(j i)/(2 R) - 1)) hat(bold(r))_(j i) $

#pagebreak()

== Rotation

=== Terminal Angular Velocity by the Drag Torque

Using the equation for the derivative of angular velocity $dot(bold(omega))_i$ plug in the equation for rotational drag.

$ J dot(bold(omega))_i = bold(T)^H_i + bold(T)^E_i + bold(T)^D_i $

$ bold(T)^D_i = - C_s bold(omega)_i $

$ J dot(bold(omega))_i = bold(T)^H_i + bold(T)^E_i - C_s bold(omega)_i $

Setting $dot(bold(omega))_i = 0$ for overdamping:
$ C_s bold(omega)_i = bold(T)^H_i + bold(T)^E_i $

Using $dot(bold(d)) = bold(omega) times bold(d)$ we get the equation for the derivative:
$ dot(bold(d))_i = bold(omega)_i times bold(d)_i = 1/C_s (bold(T)^H_i + bold(T)^E_i) times bold(d)_i $
or:
$ dot(bold(d))_i = 1/C_s (bold(T)^H_i times bold(d)_i + bold(T)^E_i times bold(d)_i) $

=== Magnetic Torque

Staring from the equation for the torque on a magnetic dipole in $bold(H)$
$
  bold(T)^H_i times bold(d)_i & = mu_0 (bold(m)_i times bold(H)_i) times bold(d)_i \
                              & = mu_0 m (bold(d)_i times bold(H)_i) times bold(d)_i
$

Using the antisymmetry of the cross product twice, the triple product identity and $norm(bold(d)_i) = 1$, we get:
$
  bold(T)^H_i times bold(d)_i & = mu_0 m bold(d)_i times (bold(H)_i times bold(d)_i) \
                              & = mu_0 m (bold(H)_i (bold(d)_i dot bold(d)_i) - bold(d)_i (bold(H)_i dot bold(d)_i)) \
                              & = mu_0 m (bold(H)_i - bold(d)_i (bold(H)_i dot bold(d)_i))
$

=== Electric Torque

Starting with the equation for the torque acting on an electric dipole in the electric field $bold(E)_i$ and using
$bold(v) times bold(v) = 0$ for any vector:
$
  bold(T)^E_i = bold(p)_i times bold(E)_i &= V_p epsilon_0 [chi_x bold(E)_i + (chi_z - chi_x) (bold(d)_i dot bold(E)_i) bold(d)_i] times bold(E)_i\
  &= V_p epsilon_0 (chi_z - chi_x) (bold(d)_i dot bold(E)_i) (bold(d)_i times bold(E)_i)\
$
Multiplying from the right with $bold(d)_i$:
$
  bold(T)^E_i times bold(d)_i &=V_p epsilon_0 (chi_z - chi_x) (bold(d)_i dot bold(E)_i) [(bold(d)_i times bold(E)_i) times bold(d)_i]\
$
similar to the equations for magnetic torque
$
  bold(T)^E_i times bold(d)_i = V_p epsilon_0 (chi_z - chi_x) (bold(d)_i dot bold(E)_i) [
    bold(E)_i - bold(d)_i (bold(E)_i dot bold(d)_i)
  ]
$

