Created on Mon Nov 29 2021

# Creator: Swamini Khurana <swamini.khurana@gmail.com>
"""

"""
## Import libraries
#%%
from sympy import symbols, Eq, Function
from sympy.solvers.ode.systems import dsolve_system

## Initialize functions
#%%
f, g = symbols("f g", cls=Function)
x = symbols("x")
eqs = [Eq(f(x).diff(x), g(x)), Eq(g(x).diff(x), f(x))]

## Solve the equations
#%%
dsolve_system(eqs)

## Solve the equations with initial conditions
#%%
dsolve_system(eqs, ics={f(0): 1, g(0): 0})

## Solve the equations specifying the dependent variables and independent variable
funcs = [f(x), g(x)]
dsolve_system(eqs, funcs=funcs, t=x)

## Implicit system of ODEs
#%%
eqs = [Eq(f(x).diff(x)**2, g(x)**2), Eq(g(x).diff(x), g(x))]
dsolve_system(eqs)
