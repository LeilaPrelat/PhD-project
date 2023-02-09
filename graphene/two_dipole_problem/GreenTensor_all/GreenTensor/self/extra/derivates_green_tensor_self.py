#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

derivadas que necesito calcular para obtener
el green tensor direct (el correcto)

"""
from sympy import *
init_printing(use_unicode=True)
#%%

x, b = symbols('x b')

#%%


self_f = diff((sqrt(x**2+b**2)-x)**2/sqrt(x**2+b**2),x)

print_latex(self_f)