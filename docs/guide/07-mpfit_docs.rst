MPFIT Documentation
===================

The text below is adapted from the ``mpfit.py`` documentation and gives
a high-level overview of the theory and math behind how the code works.
More details concerning input parameters and keywords can be found in
the `~eispac.extern.mpfit` doc string.

Description
-----------

MPFIT [#]_ uses the Levenberg-Marquardt technique to solve the least-squares
problem. In its typical use, MPFIT will be used to fit a user-supplied
function (the "model") to user-supplied data points (the "data") by
adjusting a set of parameters. MPFIT is based upon MINPACK-1 (LMDIF.F) [#]_
by More' and collaborators [#]_ [#]_.

For example, a researcher may think that a set of observed data points
is best modelled with a Gaussian curve. A Gaussian curve is parameterized
by its mean, standard deviation and normalization. MPFIT will, within
certain constraints, find the set of parameters which best fits the data.
The fit is "best" in the least-squares sense; that is, the sum of the
weighted squared differences between the model and data is minimized.

The Levenberg-Marquardt technique is a particular strategy for iteratively
searching for the best fit. This particular implementation is drawn from
MINPACK-1 (see NETLIB), and is much faster and more accurate than the
version provided in the Scientific Python package in
Scientific.Functions.LeastSquares. This version allows upper and lower
bounding constraints to be placed on each parameter, or the parameter can
be held fixed.

The user-supplied Python function should return an array of weighted
deviations between model and data. In a typical scientific problem the
residuals should be weighted so that each deviate has a gaussian sigma
of 1.0. If ``X`` represents values of the independent variable, ``Y``
represents a measurement for each value of ``X``, and ``ERR`` represents the
error in the measurements, then the deviates could be calculated as follows:

::

    DEVIATES = (Y - F(X)) / ERR

where ``F`` is the analytical function representing the model. You are
recommended to use the convenience functions MPFITFUN and MPFITEXPR,
which are driver functions that calculate the deviates for you.
If ``ERR`` are the 1-sigma uncertainties in ``Y``, then

::

    TOTAL( DEVIATES^2 )

will be the total chi-squared value. MPFIT will minimize the chi-square
value. The values of ``X``, ``Y`` and ``ERR`` are passed through MPFIT
to the user-supplied function via the ``FUNCTKW`` keyword.

Simple constraints can be placed on parameter values by using the
``PARINFO`` keyword to MPFIT. See below for a description of this keyword.

MPFIT does not perform more general optimization tasks. See TNMIN
instead. MPFIT is customized, based on MINPACK-1, to the least-squares
minimization problem.


User Function
-------------

The user must define a function which returns the appropriate values
as specified above. The function should return the weighted deviations
between the model and the data. It should also return a status flag
and an optional partial derivative array. For applications which use
finite-difference derivatives -- the default -- the user function should
be declared in the following way:

::

    def myfunct(p, fjac=None, x=None, y=None, err=None)
	    # Parameter values are passed in "p"
	    # If fjac==None then partial derivatives should not be
	    # computed.  It will always be None if MPFIT is called with
	    # default flag.
	    model = F(x, p)
	    # Non-negative status value means MPFIT should continue,
        # negative means stop the calculation.
	    status = 0
	    return([status, (y-model)/err]

See below for applications with analytical derivatives.

The keyword parameters ``X``, ``Y``, and ``ERR`` in the example above are
suggestive but not required. Any parameters can be passed to MYFUNCT by using
the ``functkw`` keyword to MPFIT. Use MPFITFUN and MPFITEXPR if you need
ideas on how to do that. The function *must* accept a parameter list, ``P``.

In general there are no restrictions on the number of dimensions in ``X``,
``Y``, or ``ERR``. However the deviates *must* be returned in a
one-dimensional Numeric array of type Float.

User functions may also indicate a fatal error condition using the status
return described above. If status is set to a number between -15 and -1
then MPFIT will stop the calculation and return to the caller.


Analytic Derivatives
--------------------

In the search for the best-fit solution, MPFIT by default calculates
derivatives numerically via a finite difference approximation. The
user-supplied function need not calculate the derivatives explicitly.
However, if you desire to compute them analytically, then the
``AUTODERIVATIVE=0`` keyword must be passed to MPFIT. As a practical matter,
it is often sufficient and even faster to allow MPFIT to calculate the
derivatives numerically, and so ``AUTODERIVATIVE=0`` is not necessary.

If ``AUTODERIVATIVE=0`` is used then the user function must check the parameter
``FJAC``, and if ``FJAC!=None`` then return the partial derivative array in the
return list.

::

    def myfunct(p, fjac=None, x=None, y=None, err=None)
	    # Parameter values are passed in "p"
	    # If FJAC!=None then partial derivatives must be comptuer.
	    # FJAC contains an array of len(p), where each entry
	    # is 1 if that parameter is free and 0 if it is fixed.
	    model = F(x, p)
	    # Non-negative status value means MPFIT should continue,
        # negative means stop the calculation.
	    status = 0
	    if (dojac):
	        pderiv = zeros([len(x), len(p)], Float)
	        for j in range(len(p)):
		        pderiv[:,j] = FGRAD(x, p, j)
	    else:
	        pderiv = None
	    return([status, (y-model)/err, pderiv]

where ``FGRAD(x, p, i)`` is a user function which must compute the derivative
of the model with respect to parameter ``P[i]`` at ``X``. When finite
differencing is used for computing derivatives (i.e., when AUTODERIVATIVE=1),
or when MPFIT needs only the errors but not the derivatives the parameter
FJAC=None.

Derivatives should be returned in the ``PDERIV`` array. ``PDERIV`` should be
an MxN array, where M is the number of data points and N is the number of
parameters. ``dp[i,j]`` is the derivative at the ith point with respect to
the jth parameter.

The derivatives with respect to fixed parameters are ignored; zero is an
appropriate value to insert for those derivatives. Upon input to the user
function, ``FJAC`` is set to a vector with the same length as ``P``, with
a value of 1 for a parameter which is free, and a value of zero for a
parameter which is fixed (and hence no derivative needs to be calculated).

If the data is higher than one dimensional, then the *last* dimension
should be the parameter dimension. Example: fitting a 50x50 image, "dp"
should be 50x50xNPAR.

.. _sec-parinfo:

Constraining Parameter Values
-----------------------------

The behavior of MPFIT can be modified with respect to each parameter to
be fitted. A parameter value can be fixed; simple boundary constraints
can be imposed; limitations on the parameter changes can be imposed;
properties of the automatic derivative can be modified; and parameters
can be tied to one another.

These properties are governed by the ``PARINFO`` structure, which is
passed as a keyword parameter to MPFIT.

``PARINFO`` should be a list of dictionaries, one list entry for each
parameter. Each parameter is associated with one element of the array,
in numerical order. The dictionary can have the following keys (none are
required, keys are case insensitive):

- **value** - the starting parameter value (but see the XALL parameter
  for more information).

- **fixed** - an "boolean" integer of 0 or 1 that determined whether the
  parameter is to be held fixed or not. If set (1), the parameter value
  will be held fixed. Fixed parameters are not varied by MPFIT, but are
  passed on to MYFUNCT for evaluation.

- **limited** - a two-element "boolean" integer array. If the first/second
  element is set, then the parameter is bounded on the lower/upper side.
  A parameter can be bounded on both sides. Both LIMITED and LIMITS must
  be given together.

- **limits** - a two-element float array. Gives the parameter limits on
  the lower and upper sides, respectively. A value will only be used as a
  limit if the corresponding value of LIMITED is set (=1) Both LIMITED and
  LIMITS must be given together.

- **parname** - a string, giving the name of the parameter. The fitting code
  of MPFIT does not use this tag in any way. However, the default iterfunct
  will print the parameter name if available.

- **step** - the step size to be used in calculating the numerical derivatives.
  If set to zero, then the step size is computed automatically. Ignored
  when AUTODERIVATIVE=0.

- **mpside** -  the sidedness of the finite difference when computing numerical
  derivatives. This field can take four values:

  0 -  one-sided derivative computed automatically

  1 -  one-sided derivative (f(x+h) - f(x) )/h

  -1 -  one-sided derivative (f(x) - f(x-h))/h

  2 -  two-sided derivative (f(x+h) - f(x-h))/(2*h)

  Where "h" is the STEP parameter described above. The "automatic"
  one-sided derivative method will chose a direction for the finite
  difference which does not violate any constraints. The other methods
  do not perform this check. The two-sided method is in principle more
  precise, but requires twice as many function evaluations. Default: 0.

- **mpmaxstep** - the maximum change to be made in the parameter value.
  During the fitting process, the parameter will never be changed by more
  than this value in one iteration. A value of 0 indicates no maximum.
  Default: 0.

- **tied** - a string expression which "ties" the parameter to other free or
  fixed parameters. Any expression involving constants and the parameter
  array P are permitted. Example: if parameter 2 is always to be twice
  parameter 1 then use the following: parinfo(2).tied = ’2 \* p(1)’.
  Since they are totally constrained, tied parameters are considered to
  be fixed; no errors are computed for them. [NOTE: the PARNAME can’t
  be used in expressions.]

- **mpprint** - if set to 1, then the default iterfunct will print the parameter
  value. If set to 0, the parameter value will not be printed. This tag
  can be used to selectively print only a few parameter values out of
  many. Default: 1 (all parameters printed)

Future modifications to the ``PARINFO`` structure, if any, will involve
adding dictionary tags beginning with the two letters "MP". Therefore
programmers are urged to avoid using tags starting with the same letters;
otherwise they are free to include their own fields within the ``PARINFO``
structure, and they will be ignored.

PARINFO Example:

::

    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0],
              'limits':[0.,0.]} for i in range(5)]
    parinfo[0]['fixed'] = 1
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limits'][0]  = 50.
    values = [5.7, 2.2, 500., 1.5, 2000.]
    for i in range(5):
      parinfo[i]['value']=values[i]

A total of 5 parameters, with starting values of 5.7, 2.2, 500, 1.5,
and 2000 are given. The first parameter is fixed at a value of 5.7,
and the last parameter is constrained to be above 50.

Example
-------

::

    import mpfit
    import numpy.oldnumeric as Numeric
    x = arange(100, float)
    p0 = [5.7, 2.2, 500., 1.5, 2000.]
    y = (p[0] + p[1]*[x] + p[2]*[x**2] + p[3]*sqrt(x) +
    	 p[4]*log(x))
    fa = {'x':x, 'y':y, 'err':err}
    m = mpfit('myfunct', p0, functkw=fa)
    print 'status = ', m.status
    if (m.status <= 0):
        print 'error message = ', m.errmsg
    print 'parameters = ', m.params

Minimizes sum of squares of MYFUNCT.  MYFUNCT is called with the X, Y,
and ERR keyword parameters that are given by FUNCTKW. The results can
be obtained from the returned object m.

Theory Of Operation
-------------------

There are many specific strategies for function minimization. One very
popular technique is to use function gradient information to realize the
local structure of the function. Near a local minimum the function value
can be taylor expanded about x0 as follows:

::

	f(x) = f(x0) + f'(x0) . (x-x0) + (1/2) (x-x0) . f''(x0) . (x-x0)
           -----   ---------------   ------------------------------- (1)
	Order	0th		     1st                       2nd

Here f'(x) is the gradient vector of f at x, and f''(x) is the Hessian
matrix of second derivatives of f at x. The vector x is the set of
function parameters, not the measured data vector. One can find the
minimum of f, f(xm) using Newton's method, and arrives at the following
linear equation:

::

    f''(x0) . (xm-x0) = - f'(x0)                                     (2)

If an inverse can be found for f''(x0) then one can solve for (xm-x0),
the step vector from the current position x0 to the new projected minimum.
Here the problem has been linearized (ie, the gradient information is known
to first order). f''(x0) is symmetric NxN matrix, and should be positive
definite.

The Levenberg-Marquardt technique is a variation on this theme. It adds an
additional diagonal term to the equation which may aid the convergence
properties:

::

    (f''(x0) + nu I) . (xm-x0) = -f'(x0)                             (2a)

where I is the identity matrix. When nu is large, the overall matrix is
diagonally dominant, and the iterations follow steepest descent. When nu
is small, the iterations are quadratically convergent.

In principle, if f''(x0) and f'(x0) are known then xm-x0 can be determined.
However the Hessian matrix is often difficult or impossible to compute.
The gradient f'(x0) may be easier to compute, if even by finite difference
techniques. So-called quasi-Newton techniques attempt to successively
estimate f''(x0) by building up gradient information as the iterations
proceed.

In the least squares problem there are further simplifications which assist
in solving eqn (2). The function to be minimized is a sum of squares:

::

    f = Sum(hi^2)                                                    (3)

where hi is the ith residual out of m residuals as described above.
This can be substituted back into eqn (2) after computing the derivatives:

::

    f'  = 2 Sum(hi  hi')
    f'' = 2 Sum(hi' hj') + 2 Sum(hi hi'')                            (4)

If one assumes that the parameters are already close enough to a minimum,
then one typically finds that the second term in f'' is negligible (or,
in any case, is too difficult to compute). Thus, equation (2) can be solved,
at least approximately, using only gradient information.

In matrix notation, the combination of eqns (2) and (4) becomes:

::

    hT' . h' . dx = - hT' . h                                        (5)

Where h is the residual vector (length m), hT is its transpose, h' is the
Jacobian matrix (dimensions n x m), and dx is (xm-x0). The user function
supplies the residual vector h, and in some cases h' when it is not found
by finite differences (see MPFIT_FDJAC2, which finds h and hT'). Even if
dx is not the best absolute step to take, it does provide a good estimate
of the best *direction*, so often a line minimization will occur along the
dx vector direction.

The method of solution employed by MINPACK is to form the Q . R
factorization of h', where Q is an orthogonal matrix such that QT . Q = I,
and R is upper right triangular. Using h' = Q . R and the ortogonality of
Q, eqn (5) becomes

::

	(RT . QT) . (Q . R) . dx = - (RT . QT) . h
                 RT . R . dx = - RT . QT . h                         (6)
                      R . dx = - QT . h

where the last statement follows because R is upper triangular. Here, R,
QT, and h are known so this is a matter of solving for dx. The routine
MPFIT_QRFAC provides the QR factorization of h, with pivoting, and
MPFIT_QRSOLV provides the solution for dx.

Authors
-------

The original version of this software, called LMFIT, was written in FORTRAN
as part of the MINPACK-1 package by XXX.

Craig Markwardt converted the FORTRAN code to IDL.
The information for the IDL version is:

    Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
    craigm@lheamail.gsfc.nasa.gov
    UPDATED VERSIONs can be found on my WEB PAGE: http://cow.physics.wisc.edu/~craigm/idl/idl.html

Mark Rivers created this Python version from Craig's IDL version.

    Mark Rivers, University of Chicago
    Building 434A, Argonne National Laboratory
    9700 South Cass Avenue, Argonne, IL 60439
    rivers@cars.uchicago.edu
    Updated versions can be found at: http://cars.uchicago.edu/software

Sergey Koposov converted the Mark's Python version from Numeric to numpy

    Sergey Koposov, University of Cambridge, Institute of Astronomy,
    Madingley road, CB3 0HA, Cambridge, UK
    koposov@ast.cam.ac.uk
    Updated versions can be found at: https://github.com/segasai/astrolibpy

Modification History
--------------------

- Translated from MINPACK-1 in FORTRAN, Apr-Jul 1998, CM

Copyright (C) 1997-2002, Craig Markwardt
This software is provided as is without any warranty whatsoever.
Permission to use, copy, modify, and distribute modified or unmodified
copies is granted, provided this copyright and disclaimer are included
unchanged.

- Translated from MPFIT (Craig Markwardt's IDL package) to Python,
  August, 2002. Mark Rivers

- Converted from Numeric to numpy (Sergey Koposov, July 2008)

- Added full Python 3 compatibility (Sergey Koposov, Feb 2017)

.. rubric:: References

.. [#] Markwardt, C. B. 2009, in Astronomical Society of the Pacific Conference
   Series, Vol. 411, Astronomical Data Analysis Software and
   Systems XVIII, ed. D. A. Bohlender, D. Durand, & P. Dowler, 251

.. [#] MINPACK-1, Jorge More', available from netlib (www.netlib.org).

.. [#] "Optimization Software Guide," Jorge More' and Stephen Wright, SIAM,
   *Frontiers in Applied Mathematics*, Number 14.

.. [#] More', Jorge J., "The Levenberg-Marquardt Algorithm: Implementation
   and Theory," in *Numerical Analysis*, ed. Watson, G. A., Lecture
   Notes in Mathematics 630, Springer-Verlag, 1977.
