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

I0_6_aux1 = diff((x**2+b**2)**(-1/2),x)

I0_6_aux2 = diff(I0_6_aux1,x)

###

I2_6_aux1 = diff((sqrt(x**2+b**2)-x)**2*(x**2+b**2)**(-1/2),x)

I2_6_aux2 = diff(I2_6_aux1,x)

print_latex(I2_6_aux2)