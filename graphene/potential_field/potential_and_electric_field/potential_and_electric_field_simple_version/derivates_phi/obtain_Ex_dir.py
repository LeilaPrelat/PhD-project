#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

derivadas que necesito para el campo electrico

derivadas del potential electrico sin usar hankel 
para el dipolo en xD = yD = 0
"""
from sympy import *
import sys
import os 

init_printing(use_unicode=True)

#%%

#name_this_py = os.path.basename(__file__)
#path = os.path.abspath(__file__) #path absoluto del .py actual
#path_basic = path.replace('/' + name_this_py,'')
#path_constants =  path_basic.replace('/potential_field','')
#print('Importar modulos necesarios para este codigo')

#%%


x = symbols('x')
px, pz = symbols('p_x p_z')
zp, b = symbols('zp b')


term1 = Abs(b)**2*(Abs(b)**2 + Abs(x)**2)**(-3/2)
term2 = (Abs(b)**2 + Abs(x)**2)**(-1/2)
term3 = Abs(b)*(Abs(b)**2 + Abs(x)**2)**(-3/2)


term_1_final = -px*(term1 - term2)/Abs(x)



Efield_x = -diff( term_1_final - pz*sign(b)*term3, x)
    
#print_latex(Efield_x)
print_latex(Efield_x)


# Efield_x2 =  Efield_x.subs(sqrt((x+xD)**2 + (y+yD)**2),R2)
# print_latex(Efield_x)

# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

