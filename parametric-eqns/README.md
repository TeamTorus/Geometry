This [document](./README.pdf) describes the mathematical approach to our generated
parametric model for a toroidal propeller. A fully parametric propeller
model will allow us to mesh the geometry at any resolution without
having to interpolate, as well as defines a more convenient translation
between input parameter spaces and geometric changes in a propeller.
This document will outline how the cross-sectional 2D shape is defined,
before explaining the remaining math that allows for the airfoil to be
remapped and cast along a 3D curve to create the toroidal propeller
blade.

# Defining the 2D Airfoil Shape

The airfoil shape is defined using the NACA 4-digit series equations,
which provide a standardized method for generating airfoil geometries
based on parameters controlling camber and thickness.

## NACA 4-Digit Airfoil Equations

The mean camber line *y*<sub>*c*</sub>(*t*) is defined piecewise:

$$y_c(t) = \begin{cases}
\frac{m}{p^2}(2pt - t^2), & 0 \leq t \leq p \\
\frac{m}{(1 - p)^2}\[(1 - 2p) + 2pt - t^2\], & p \< t \leq 1
\end{cases}$$

The thickness distribution *y*<sub>*t*</sub>(*t*) is given by:

$$y_t(t) = 5 \times \text{thickness} \times \left\[0.2969\sqrt{t} - 0.1260t - 0.3516t^2 + 0.2843t^3 - 0.1036t^4\right\]$$

where:

-   *t* is the normalized chord length parameter (0 ≤ *t* ≤ 1),

-   *m* is the maximum camber,

-   *p* is the location of maximum camber,

-   thickness is the maximum thickness as a fraction of the chord.

## Constructing Upper and Lower Surfaces

The coordinates of the upper and lower surfaces are calculated as:

$$\begin{aligned}
    x_u(t) &= t - y_t(t) \sin(\theta_c(t)) \\
    y_u(t) &= y_c(t) + y_t(t) \cos(\theta_c(t)) \\
    x_l(t) &= t + y_t(t) \sin(\theta_c(t)) \\
    y_l(t) &= y_c(t) - y_t(t) \cos(\theta_c(t))
\end{aligned}$$

where *θ*<sub>*c*</sub>(*t*) is the angle of the camber line slope:

$$\theta_c(t) = \arctan\left(\frac{dy_c}{dt}\right)$$

The NACA airfoil definition requires the thickness to be applied
perpendicular to the camber line, requiring the above consideration.
However, for the sake of simplicity and computational efficiency
(particularly for the normalization of mesh sizes),
*θ*<sub>*c*</sub>(*t*) can be set to zero, assuming negligible camber
line slope or that the model already represents a sufficient portion of
the airfoil’s geometry space without needing this consideration. When
*θ*<sub>*c*</sub>(*t*) = 0, the thickness is just applied in the
*y*-direction.

## Parametrization Over *t*

The parameter *t* is divided into two intervals:

-   \[0, 1) for the upper surface,

-   (1, 2\] for the lower surface (shifted to maintain continuity).

The complete airfoil coordinates are defined using Heaviside functions:

$$\begin{aligned}
    x\_{\text{airfoil}}(t) &= x_u(t) \cdot H(1 - t) + x_l(2 - t) \cdot H(t - 1) \\
    y\_{\text{airfoil}}(t) &= y_u(t) \cdot H(1 - t) + y_l(2 - t) \cdot H(t - 1)
\end{aligned}$$

where *H* is the Heaviside step function. This ensures that the airfoil
function, **f**(*t*) = (*x*(*t*), *y*(*t*)), is an injective function
with respect to parameter *t*.

# Transformations Along the Blade Centerline

To account for changes in blade orientation and cross-sectional size
along the centerline, rotation and scaling transformations are applied
as functions of the parameter *s*, which parametrizes the centerline
curve.

## Rotation Transformation

The rotation matrix *R*(*s*) for an angle of attack *α*(*s*) is:

$$R(s) = \begin{bmatrix}
\cos\alpha(s) & -\sin\alpha(s) \\
\sin\alpha(s) & \cos\alpha(s)
\end{bmatrix}$$

The angle of attack *α*(*s*) can be defined as a polynomial function of
*s*:

*α*(*s*) = *a*<sub>*α*</sub>*s*<sup>4</sup> + *b*<sub>*α*</sub>*s*<sup>3</sup> + *c*<sub>*α*</sub>*s*<sup>2</sup> + *d*<sub>*α*</sub>*s* + *e*<sub>*α*</sub>

## Scaling Transformation

Scaling factors *S*<sub>*x*</sub>(*s*) and *S*<sub>*y*</sub>(*s*) modify
the airfoil dimensions:

$$\begin{aligned}
    S_x(s) &= a\_{sx}s^4 + b\_{sx}s^3 + c\_{sx}s^2 + d\_{sx}s + e\_{sx} \\
    S_y(s) &= a\_{sy}s^4 + b\_{sy}s^3 + c\_{sy}s^2 + d\_{sy}s + e\_{sy}
\end{aligned}$$

These functions allow for parametric scaling along the blade’s path, as
they are parameterized by *s*, the parameter for the centerline of the
toroidal blade. The two scaling functions scale the airfoil’s *x* and
*y* dimensions. Traditional propellers and other definitions of toroidal
propellers may represent the airfoil transformations along the blade by
performing rotations in all three axis. Instead, for the sake of
simplicity, as well as being able to place every 2D airfoil perfectly
perpendicular to the centerline, the rotation about *z* (the yaw) is
represented with the angle of attack function, and for roll and pitch
rotations, the scaling functions are used as a proxy. This is because
for any roll or pitch, you can represent it by its projection on the
plane perpendicular to the centerline by simply scaling either the *x*
or *y* dimension, respectively, by *c**o**s*(*ϕ*), where *ϕ* is the
angle to have been rotated by in that axis.

## Applying Transformations

The transformed airfoil coordinates are:

$$\begin{bmatrix}
X\_{\text{rotated}}(s, t) \\
Y\_{\text{rotated}}(s, t)
\end{bmatrix}
= R(s) \cdot
\begin{bmatrix}
x\_{\text{airfoil}}(t) \\
y\_{\text{airfoil}}(t)
\end{bmatrix}$$

Scaling is applied as:

$$\begin{aligned}
    X\_{\text{scaled}}(s, t) &= S_x(s) \cdot X\_{\text{rotated}}(s, t) \\
    Y\_{\text{scaled}}(s, t) &= S_y(s) \cdot Y\_{\text{rotated}}(s, t)
\end{aligned}$$

## Function Space Analysis for Transformation Parameters

In designing the transformations for angle of attack, *x*-scale, and
*y*-scale along the blade centerline, it is essential to select function
forms that can model a wide variety of possible transformation behaviors
while minimizing the number of parameters required. This section details
the analysis performed to determine suitable function forms that can
approximate a large function space relevant to airfoil transformations
in toroidal propellers.

### Objective

The goal is to identify a function form *f*(*s*) that:

-   Can approximate a wide range of possible transformation behaviors
    along the blade centerline.

-   Uses a minimal number of parameters to reduce computational
    complexity and overfitting risks.

To achieve this, a Monte Carlo simulation was performed where randomly
generated functions, representing possible transformation behaviors,
were fitted using candidate function forms. The performance of each
candidate function form was evaluated based on its ability to fit the
randomly generated functions.

### Random Function Generators

A set of random function generators that produce diverse function
behaviors was defined to model a sample space that a function may try to
fit. Here, *X* ∼ 𝒰(*a*, *b*) represents a parameter *X* randomly sampled
from a uniform distribution with lower and upper bounds of *a* and *b*,
respectively.

1.  **Random Sinusoidal Function**:
    *f*<sub>sin</sub>(*s*) = *A*sin (2*π**F**s* + *ϕ*) where
    *A* ∼ 𝒰(0.5, 2) is the amplitude, *F* ∼ 𝒰(1, 3) is the frequency,
    and *ϕ* ∼ 𝒰(0, *π*) is the phase shift.

2.  **Random Polynomial Function**:
    $$f\_{\text{poly}}(s) = \sum\_{i=0}^{d} c_i s^i$$
    where *d* ∼ 𝒰{1, 8} is the degree, and *c*<sub>*i*</sub> ∼ 𝒰(−2, 2)
    are the coefficients.

3.  **Random Logistic Function**:
    $$f\_{\text{logistic}}(s) = \frac{L}{1 + e^{-k(s - s_0)}}$$
    where *L* ∼ 𝒰(1, 50) is the carrying capacity, *k* ∼ 𝒰(8, 100) is
    the growth rate, and *s*<sub>0</sub> ∼ 𝒰(0.1, 0.9) is the midpoint.

4.  **Random Piecewise Linear Function**:
    *f*<sub>interp</sub>(*s*) = Interpolate{(*s*<sub>*i*</sub>, *y*<sub>*i*</sub>)}
    where *s*<sub>*i*</sub> are *n*′ points uniformly sampled in
    \[0, 1\] and *y*<sub>*i*</sub> ∼ 𝒰(−1, 1). To avoid excessively
    chaotic behavior, the points taken from *y*<sub>*i*</sub> ∼ 𝒰(−1, 1)
    is downsampled in size so as to have some continuity across *s*.
    What this means is that regardless of how many sample data points
    are being generated, *n*′ ∼ 𝒰(3, 8), with the rest of the data
    points being the result of linear interpolations between those
    sampled data points.

These function generators were chosen to represent behaviors such as
cyclic patterns, sharp transitions, and random fluctuations, which are
characteristic of possible airfoil transformations along a toroidal
blade.

### Candidate Fitting Functions

The following function forms were considered to fit the randomly
generated data:

1.  **4th-Degree Polynomial**:
    *f*<sub>poly4</sub>(*s*) = *a*<sub>4</sub>*s*<sup>4</sup> + *a*<sub>3</sub>*s*<sup>3</sup> + *a*<sub>2</sub>*s*<sup>2</sup> + *a*<sub>1</sub>*s* + *a*<sub>0</sub>

2.  **6th-Degree Polynomial**:
    *f*<sub>poly6</sub>(*s*) = *a*<sub>6</sub>*s*<sup>6</sup> + *a*<sub>5</sub>*s*<sup>5</sup> + *a*<sub>4</sub>*s*<sup>4</sup> + *a*<sub>3</sub>*s*<sup>3</sup> + *a*<sub>2</sub>*s*<sup>2</sup> + *a*<sub>1</sub>*s* + *a*<sub>0</sub>

3.  **Fourier Series (Order 2)**:
    $$f\_{\text{Fourier2}}(s) = A_0 + \sum\_{k=1}^{2} \left\[A_k \cos(2\pi k s) + B_k \sin(2\pi k s)\right\]$$

4.  **Fourier Series (Order 3)**:
    $$f\_{\text{Fourier3}}(s) = A_0 + \sum\_{k=1}^{3} \left\[A_k \cos(2\pi k s) + B_k \sin(2\pi k s)\right\]$$

5.  **Sinusoidal Function**:
    *f*<sub>sin</sub>(*s*) = *A*sin (*B**s* + *C*) + *D*

These functions were chosen to represent different levels of complexity
and flexibility, with varying numbers of degrees of freedom. Thus, only
fitting functions with similar degrees of freedom will be compared, with
dimensionality being a external factor to consider (4th-Degree
Polynomial vs Order 2 Fourier Series and 6th-Degree Polynomial vs Order
3 Fourier Series). The sinusoidal function is fitted as a reference.

### Fitting Procedure

For a randomly generated function *f*<sub>random</sub>(*s*):

1.  Generate *n* data points {(*s*<sub>*i*</sub>, *y*<sub>*i*</sub>)}
    where *s*<sub>*i*</sub> are uniformly sampled in \[0, 1\] and
    *y*<sub>*i*</sub> = *f*<sub>random</sub>(*s*<sub>*i*</sub>).

2.  Fit each candidate function *f*<sub>candidate</sub>(*s*) to the data
    using least squares optimization to minimize the mean squared error
    (MSE):
    $$\text{MSE} = \frac{1}{n} \sum\_{i=1}^{n} \left\[y_i - f\_{\text{candidate}}(s_i)\right\]^2$$

3.  Record the MSE and fitting parameters for each candidate function.

This process was repeated for a large number of iterations (*N* = 10000)
to sample a wide variety of random functions.

## Results and Analysis

The average MSE for each candidate function over all iterations is
summarized below:

|  **Candidate Function**  | **Mean Error (MSE)** | **Number of Parameters** |
|:------------------------:|:--------------------:|:------------------------:|
|  4th-Degree Polynomial   |        2.064         |            5             |
|  6th-Degree Polynomial   |        0.845         |            7             |
| Fourier Series (Order 2) |        2.441         |            5             |
| Fourier Series (Order 3) |        1.193         |            7             |
|   Sinusoidal Function    |        5.748         |            4             |

Mean squared error (MSE) for candidate fitting functions over 10,000
iterations.

### Interpretation

-   The **6th-degree polynomial** achieved the lowest average MSE,
    indicating a superior ability to fit the random functions, at the
    cost of increased complexity (7 parameters).

-   The **4th-degree polynomial** had a moderate MSE while using fewer
    parameters (5 parameters), striking a balance between flexibility
    and simplicity.

-   The **Fourier series** functions, while theoretically capable of
    modeling periodic behaviors, performed worse than the polynomials in
    this context, possibly due to the non-periodic nature of some random
    functions generated.

-   The **sinusoidal function** had the highest MSE, suggesting it is
    insufficiently flexible to model the diverse function space
    considered.

Considering the trade-off between model complexity and fitting accuracy,
the **4th-degree polynomial** was somewhat holistically selected for
modeling the transformation functions. It provides adequate flexibility
to approximate a wide range of behaviors with a manageable number of
parameters. Higher-degree polynomials or functions with more parameters
could lead to overfitting and increased computational burden without
substantial gains in modeling capability for the intended application.

# Generating the Centerline Curve with Splines

A Non-Uniform Rational B-Spline (NURBS) curve is used to define the
blade’s centerline due to its flexibility and precision in modeling
complex natural shapes.

## Coordinate System and Control Point Definition

The generation of the blade centerline begins by dividing the lateral
face of the cylindrical hub into *n* equal sections, where *n*
corresponds to the number of blades. For a section, a local 2D
coordinate system is established:

-   The positive *x*-axis (+*x*) extends horizontally across the
    section.

-   The positive *y*-axis (+*y*) extends vertically downward towards the
    bottom of the cylinder.

This coordinate system facilitates the placement and manipulation of
control points within each blade’s sector. The origin, (0, 0) is placed
1 unit down from the top left of the section. This leaves some padding
so that the blade isn’t intersecting the top of the cylindrical hub.
This 2D coordinate system can be extended to 3D by adding a *z*-axis,
creating a coordinate system denoted *L*.

-   The positive *z*-axis (+*z*) extends radially out from the cylinder.
    For any point (*x*, *y*), the scaled *z*-vector will always normal
    to the lateral cylinder face.

While extending the *z*-axis radially creates a non-Euclidean coordinate
frame, it creates a natural relation between the *z*-coordinate and the
radial distance of control points, which can help with optimization.

## Centerline Inset and Hub Intersection

To facilitate the boolean union operation between the blade and the hub,
the centerline is inset slightly into the hub’s volume. This is achieved
by setting:

$$R\_{\text{inset}} = R\_{\text{hub}} - \frac{R\_{\text{hub}}}{8}$$

This inset ensures that the generated blade intersects with the hub,
allowing for a seamless geometric union that results in a single,
manifold mesh suitable for simulation and analysis.

## Control Points Placement

A cubic NURBS spline requires control points that define the curve’s
shape. For each blade, four control points are defined:

-   *P*<sub>1</sub>: (0, 0, 0)<sub>*L*</sub>

-   *P*<sub>2</sub>:
    (*x*<sub>2</sub>, *y*<sub>2</sub>, *z*<sub>2</sub>)<sub>*L*</sub>

-   *P*<sub>3</sub>:
    (*x*<sub>3</sub>, *y*<sub>3</sub>, *z*<sub>3</sub>)<sub>*L*</sub>

-   *P*<sub>4</sub>: (*x*<sub>4</sub>, *y*<sub>4</sub>, 0)<sub>*L*</sub>

1.  **First Control Point *P*<sub>1</sub>**: The coordinates of
    *P*<sub>1</sub> in the global coordinate system are:
    $$P_1 = \left(R\_{\text{inset}}, \\ 0, \\ \frac{L\_{\text{hub}}}{2} - 1\right)$$
    This is taken such that the local origin for a blade’s lateral face
    section lies on the global *x*-axis. Thus the coordinate for
    *P*<sub>1</sub> in the global cylindrical and Cartesian frame are
    the same numerically.

2.  **Local to Global Coordinate Transformation**: Let *C* be the global
    cylindrical frame, represented by coordinates
    (*r*, *θ*, *z*)<sub>*C*</sub>. Coordinate frame *L* is simply a
    remapping of this cylindrical frame, where change in
    *x*<sub>*L*</sub> is the same as a change in *θ*<sub>*C*</sub>, a
    change in *y*<sub>*L*</sub> is the same as a change in
    −*z*<sub>*C*</sub>, and a change in *z*<sub>*L*</sub> is the same as
    a change in *r*<sub>*C*</sub>. Thus, the change of basis matrix
    becomes:

$$\mathbf{M} = \begin{bmatrix}
\mathbf{e}\_{x\_{\mathit{L}}} & \mathbf{e}\_{y\_{\mathit{L}}} & \mathbf{e}\_{z\_{\mathit{L}}}
\end{bmatrix} =
\begin{bmatrix}
0 & 0 & 1 \\
1 & 0 & 0 \\
0 & -1 & 0
\end{bmatrix}$$

Therefore, for local control point *P*<sub>*i*</sub>, its relative
global cylindrical coordinate is found with:
**P**<sub>**i**</sub><sub>*C**′*</sub> = **M** ⋅ **P**<sub>**i**</sub><sub>*L*</sub>
To get the global coordinates for control point *P*<sub>*i*</sub>,
the offset of origins needs to be considered, so adding the origin
vector, which is just **P**<sub>**1**</sub>, will give the absolute
coordinates.
**P**<sub>**i**</sub><sub>*C*</sub> = **P**<sub>**i**</sub><sub>*C′*</sub> + **P**<sub>**1**</sub><sub>*C*</sub>

## NURBS Curve Construction

With the control points defined, the cubic NURBS spline *C*(*s*) is
constructed as:

$$C(s) = \frac{\sum\_{i=0}^{3} N\_{i,3}(s) w_i P_i}{\sum\_{i=0}^{3} N\_{i,3}(s) w_i}$$

where:

-   *N*<sub>*i*, 3</sub>(*s*) are the cubic B-spline basis functions,

-   *w*<sub>*i*</sub> are the weights associated with each control
    point. For this case, so as to reduce the parameters needed for
    optimization, the weights are fixed with *w*<sub>*i*</sub> = 1

-   *s* is the normalized parameter along the spline (0 ≤ *s* ≤ 1).

The knot vector {*u*<sub>*i*</sub>} for a cubic NURBS with clamped ends
is:

{*u*<sub>*i*</sub>} = {0, 0, 0, 0, 1, 1, 1, 1}

## Spline Fitting and Symmetry Considerations

The cubic NURBS spline provides the necessary flexibility to define
smooth and continuous centerline curves for each blade while maintaining
rotational symmetry across the hub. By defining control points in a
non-Euclidean coordinate space, where the *z*-axis is radial with
respect to the hub, the resulting centerline ensures that each airfoil
section remains perpendicular to the centerline.

The spline fitting process involves adjusting the control points
*P*<sub>1</sub> and *P*<sub>2</sub> to achieve the desired curvature and
alignment of the blade. By manipulating these control points, complex
blade geometries can be accurately modeled with a minimal number of
parameters.

# Extruding the Airfoil Using Frenet-Serret Frames

The Frenet-Serret frame provides an orthonormal basis at each point
along the curve, allowing the airfoil to be correctly oriented
perpendicular to the centerline. This is a relative path coordinate
system that moves along the centerline to orient and define the airfoil
through.

## Calculating the Frenet-Serret Frame

At each point *s* along the curve *C*(*s*):

1.  **Tangent Vector *T*(*s*)**:
    $$T(s) = \frac{C'(s)}{\left\| C'(s) \right\|}$$

2.  **Normal Vector *N*(*s*)**:
    $$N(s) = \frac{T'(s)}{\left\| T'(s)\right\|}$$

3.  **Binormal Vector *B*(*s*)**:
    *B*(*s*) = *T*(*s*) × *N*(*s*)

## Positioning the Airfoil

The final 3D coordinates of the airfoil are:

$$\begin{aligned}
    X\_{\text{final}}(s, t) &= X_c(s) + X\_{\text{scaled}}(s, t) N_x(s) + Y\_{\text{scaled}}(s, t) B_x(s) \\
    Y\_{\text{final}}(s, t) &= Y_c(s) + X\_{\text{scaled}}(s, t) N_y(s) + Y\_{\text{scaled}}(s, t) B_y(s) \\
    Z\_{\text{final}}(s, t) &= Z_c(s) + X\_{\text{scaled}}(s, t) N_z(s) + Y\_{\text{scaled}}(s, t) B_z(s)
\end{aligned}$$

where (*X*<sub>*c*</sub>, *Y*<sub>*c*</sub>, *Z*<sub>*c*</sub>) are the
coordinates of the centerline *C*(*s*), and
(*N*<sub>*x*</sub>, *N*<sub>*y*</sub>, *N*<sub>*z*</sub>),
(*B*<sub>*x*</sub>, *B*<sub>*y*</sub>, *B*<sub>*z*</sub>) are components
of *N*(*s*) and *B*(*s*).

# Generating Evenly Spaced Mesh Points

Using a linear spacing for the airfoil parameter, *t*, yields a mesh
with uneven cell sizes. Since
$$\frac{\mathrm{d}^2 \mathbf{f}}{\mathrm{d} t^2} \neq 0 \text{ where } \mathbf{f}(t) = (x(t), y(t))$$
and **f** is the vector airfoil function, the derivative of **f**
changes along the airfoil’s surface. This means that an even linearly
sampled set for the input will have the output’s spacing be
inconsistent. This would cause the outputted mesh to have differential
cell sizes, with the leading edge, or regions where
$\|\frac{\mathrm{d} \mathbf{f}}{\mathrm{d} t}\|$ is the largest, having
larger cells compared to the trailing edge.

Uniformly distributed mesh points along the airfoil perimeter ensure
consistent cell sizes in the final mesh, improving CFD and meshing
accuracy. The goal is to find some set or vector of parameter values,
**t**, that when plugged into the symbolic equations will create even
cell sizes in the mesh. Even sizes of meshes are determined by covering
the same arc distance along the airfoil.

## Computing Arc Length Parameterization

The arc length *s(t*) along the airfoil is:

$$s(t) = \int\_{t_0}^{t} \sqrt{\left(\frac{dx}{dt}\right)^2 + \left(\frac{dy}{dt}\right)^2} \\ dt$$

Due to the complexity of the NACA equations, numerical integration is
employed.

## Generating *t* Values Corresponding to Equal Arc Length Increments

1.  **Total Arc Length *S*<sub>total</sub>**:
    *S*<sub>total</sub> = *s*(*t*<sub>end</sub>)

2.  **Arc Length Increment *Δs*:
    $$\Delta s = \frac{S\_{\text{total}}}{n}$$
    where *n* is the desired number of intervals (mesh resolution).

3.  **Iterative Calculation of *t* Values**:

    Starting from *t*<sub>0</sub>, solve for *t*<sub>*i*</sub> such
    that:
    *s*(*t*<sub>*i*</sub>) − *s*(*t*<sub>*i* − 1</sub>) = *Δs*
    This is not analytically solvable, so numerical methods are used to
    approximate the solution. This is done using a resolution of
    *n*<sub>*fine*</sub> iterations over the entire airfoil edge,
    however this can be refined further to be more accurate, albeit with
    non-negligible computational costs. This integration at each
    iteration is done using the technique from the QUADPACK Fortran
    library, called by a SciPy proxy function. It uses adaptive Gaussian
    quadrature methods to approximate the integral.

4.  **Generate Mapping Function *t*(*s*)**:

    A linear interpolation between the *n*<sub>*fine*</sub> (by
    default 1000) integrated data points yields a function *t*(*s*),
    which when inputted the current cumulative arc distance, outputs the
    t value that would produce the current coordinate point. Evaluating

$$t(\mathbf{\overline{s}}) = \mathbf{t}$$
where
$$\mathbf{\overline{s}} = {s_{i} = i \times \Delta s | i \in [0, n]\}$$

yields the vector **t**, whose elements are *t*-values that will
produce even mesh spacing.

# Assembling the Complete Propeller

## Duplicating Blades

The blade geometry is duplicated *n* times (as specified by the number
of blades parameter), rotated evenly around the hub’s axis:

$$\theta_i = \theta_0 + \frac{2\pi i}{n}$$

## Generating the Hub

The hub is modeled as a cylinder with radius *R*<sub>hub</sub> and
length *L*<sub>hub</sub>. A parametric representation in cylindrical
coordinates is used:

$$\begin{aligned}
    X\_{\text{hub}}(\theta, z) &= R\_{\text{hub}} \cos\theta \\
    Y\_{\text{hub}}(\theta, z) &= R\_{\text{hub}} \sin\theta \\
    Z\_{\text{hub}}(\theta, z) &= z
\end{aligned}$$

# Mesh Generation and Boolean Operations

## Creating Mesh Faces

Using the computed coordinates, mesh faces are generated by connecting
adjacent points in the parameter grids. This is connected in an order
such that the cross product for the normals of the face will point
outwards.

## Triangulation

At this point, there are *n* + 1 bodies, where *n* corresponds to the
number of blades and the extra body is the hub. In order to merge them
into a contiguous mesh, PyVista is used. This library requires bodies to
first be triangulated, which is done by simply cutting each rectangular
mesh face along the diagonal.

## Geometric Union

The blades and hub are combined into a single mesh using boolean union
operations in PyVista, resulting in a manifold mesh suitable for
simulations.

# Appendix: Mathematical Expressions and Parameters

Default parameters for now

## Airfoil Parameters

-   Maximum camber *m* = 0.02

-   Location of maximum camber *p* = 0.4

-   Maximum thickness thickness = 0.5

## Transformation Coefficients

-   Angle of attack coefficients
    *a*<sub>*α*</sub>, *b*<sub>*α*</sub>, *c*<sub>*α*</sub>, *d*<sub>*α*</sub>, *e*<sub>*α*</sub>

-   Scaling coefficients
    *a*<sub>*s**x*</sub>, *b*<sub>*s**x*</sub>, *c*<sub>*s**x*</sub>, *d*<sub>*s**x*</sub>, *e*<sub>*s**x*</sub>
    and
    *a*<sub>*s**y*</sub>, *b*<sub>*s**y*</sub>, *c*<sub>*s**y*</sub>, *d*<sub>*s**y*</sub>, *e*<sub>*s**y*</sub>

## NURBS Parameters

-   Degree *p* = 3

-   Knot vector {0, 0, 0, 0, 1, 1, 1, 1}

-   Control points *P*<sub>*i*</sub>

-   Weights *w*<sub>*i*</sub> = 1 (uniform)
