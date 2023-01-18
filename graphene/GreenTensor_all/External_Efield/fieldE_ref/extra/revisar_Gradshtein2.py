
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Revisar la integral
Ref \cite{Gradshteyn} page 718 6.671 8.
con sympy
"""
#from scipy import integrate
from sympy import *

x = Symbol('x')
a = Symbol('a')
b = Symbol('b')

rta = integrate(cos(b*x)*besselj(0,a*x), (x, 0, oo)) 
init_printing(use_unicode=True)

#%%
