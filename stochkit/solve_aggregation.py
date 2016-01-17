from sympy.solvers import solve
from sympy import solve_poly_system
from sympy import Symbol
import math

a1 = 1.0
a2 = 1.0
b1 = 1.0
b2 = 1.0

z1 = 5.0
z2 = 6.0
z4 = 8.0


def get_stationary_c(z1, z2, z4, a1=1., a2=1., b1=1., b2=1.):
    x3 = Symbol('x3', real=True, nonnegative=True)
    x4 = Symbol('x4', real=True, nonnegative=True)

    solutions = solve(
        [
            a1*(z1-x3-z4+x4)*(z2-x3-z4+x4)-a2*x3,
            b1*x3*x4-b2*(z4-x4)
        ],
        dict=True)

    print solutions
    
    l =len(solutions)
    cnt = 0
    for i in range(l): 
        _x3 = float(solutions[i][x3])
        _x4 = float(solutions[i][x4])
        x1 = z1 - _x3 - z4 + _x4
        x2 = z2 - _x3 - z4 + _x4
        x5 = z4 - _x4
        
        if (x1>0. and x2>0. and x5>0.):
            cnt += 1
            sol = i
            _x1 = x1
            _x2 = x2
            _x5 = x5
            sol = [_x1, _x2, _x3, _x4, _x5]

    if cnt == 1:        
        return(sol)
    else:
        print "Error"
        return
        
        
def get_stationary_distr(z1, z2, z4):
    sol = get_stationary_c(z1, z2, z4)    
    print sol
    
    prob = dict()
    
    # compute all possible values for x_i according to rule:
    # x3 == [0, min(z1, z2)]
    # x4 == [0, z4]
    for x3 in range(int(min(z1,z2))):
        for x4 in range(int(z4)):
            x1 = z1 - x3 - z4 + x4
            x2 = z2 - x3 - z4 + x4
            x5 = z4 - x4
            currx = [x1, x2, x3, x4, x5]
            # need to check that all xi are positive, else skip this loop
            pos = (x1>0 and x2>0 and x3>0 and x4>0 and x5>0)
            if not pos: continue
            
            # compute formula for stationary distribution f(c_i, x_i):
            p = 1
            for i in range(len(sol)):
                ci = sol[i]
                xi = currx[i]
                #p *= ci^xi / xi! * exp(-ci):
                p *= (math.pow(ci,xi) / math.factorial(xi) * math.exp(-ci))
                # use dummy function to test:
                #p *= ci * xi
                
            # return dict with key(x_i values) value(eval pi(x_i values))
            prob[(x3, x4)] = p
            
    return prob
    
p = get_stationary_distr(z1, z2, z4)
print p
