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

alfa_p = symbols('\alpha_p')
kp = symbols('k_p')
x, y, z = symbols('x y z')
px, py, pz = symbols('p_x p_y p_z')
zp, b = symbols('zp b')
kp_2 = symbols('k^2_p')
Rp = symbols('R_p')


R = sqrt(x**2 + y**2)
phi = atan(y/x)
J1 = besselj(1,kp*R)
J0 = besselj(0,kp*R)
expo_f = exp(-kp*(2*zp + Abs(z) + Abs(b)))

term0_ref = (Abs(z) + Abs(b))**2
term1_ref = term0_ref*((term0_ref + R**2)**(-3/2)) - (term0_ref + R**2)**(-1/2)

term0_dir = Abs(z - zD)**2
term1_dir = term0_dir*((term0_dir + R**2)**(-3/2)) - (term0_dir + R**2)**(-1/2)
term_px_py = px*cos(phi) + py*sin(phi)

term1_final = -(term1_dir + term1_ref)*term_px_py/R 

term2_final = -Rp*kp_2*term_px_py*J1*expo_f

term3_final = sign(z)*pz*(term0_ref*((term0_ref + R**2)**(-3/2)) + term0_dir*((term0_dir + R**2)**(-3/2)))
term4_final = sign(z)*pz*(Rp*kp_2*J0*expo_f)

term_tot_final = term1_final + term2_final + term3_final + term4_final 

Efield_x = -diff( term_tot_final, x)
Efield_y = -diff( term_tot_final, y)
Efield_z = -diff( term_tot_final, z)
    

    
print_latex(Efield_x)


# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

