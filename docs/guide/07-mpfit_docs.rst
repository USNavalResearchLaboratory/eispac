Key mpfit documentation
=======================

.. _sec-parinfo:

parinfo keywords
----------------

**Constraining Parameter Values with the PARINFO Keyword**


The behavior of MPFIT can be modified with respect to each parameter to
be fitted. A parameter value can be fixed; simple boundary constraints
can be imposed; limitations on the parameter changes can be imposed;
properties of the automatic derivative can be modified; and parameters
can be tied to one another.

These properties are governed by the PARINFO structure, which is passed
as a keyword parameter to MPFIT.

PARINFO should be a list of dictionaries, one list entry for each
parameter. Each parameter is associated with one element of the array,
in numerical order. The dictionary can have the following keys (none are
required, keys are case insensitive):

- **value - the starting parameter value (but see the START_PARAMS parameter
  for more information).

- **fixed** - a boolean value, whether the parameter is to be held fixed or not.
  Fixed parameters are not varied by MPFIT, but are passed on to MYFUNCT for
  evaluation.

- **limited** - a two-element boolean array. If the first/second element is set,
  then the parameter is bounded on the lower/upper side. A parameter can
  be bounded on both sides. Both LIMITED and LIMITS must be given
  together.

- **limits** - a two-element float array. Gives the parameter limits on the lower
  and upper sides, respectively. Zero, one or two of these values can
  be set, depending on the values of LIMITED. Both LIMITED and LIMITS
  must be given together.

- **parname** - a string, giving the name of the parameter. The fitting code of
  MPFIT does not use this tag in any way. However, the default iterfunct
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

- **mpmaxstep** -  the maximum change to be made in the parameter value. During
  the fitting process, the parameter will never be changed by more than
  this value in one iteration. A value of 0 indicates no maximum.
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
