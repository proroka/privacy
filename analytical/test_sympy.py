from sympy.solvers import solve
from sympy import solve_poly_system
from sympy import Symbol

a1 = 1.0
a2 = 1.0
b1 = 1.0
b2 = 1.0

z1 = 10.0
z2 = 15.0
z4 = 20.0

#x3 = Symbol('x3')
#x4 = Symbol('x4')
#x3 = Symbol('x3', real=True, nonnegative=True)
#x4 = Symbol('x4', real=True, nonnegative=True)

x1 = Symbol('x1', real=True, nonnegative=True)
x2 = Symbol('x2', real=True, nonnegative=True)
x3 = Symbol('x3', real=True, nonnegative=True)
x4 = Symbol('x4', real=True, nonnegative=True)
x5 = Symbol('x5', real=True, nonnegative=True)
solutions = solve(
   [
       a1*(z1-x3-z4+x4)*(z2-x3-z4+x4)-a2*x3,
       b1*x3*x4-b2*(z4-x4),
       z1-x1-x3-x5,
       z2-x2-x3-x5,
       z4-x5,
   ],
   dict=True)

#solutions = solve(
#    [
#        a1*(z1-x3-z4+x4)*(z2-x3-z4+x4)-a2*x3,
#        b1*x3*x4-b2*(z4-x4)
#    ],
#    dict=True)
#
print solutions
