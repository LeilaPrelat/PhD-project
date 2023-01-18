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
k1 = symbols('k1')
px, py, pz = symbols('p_x p_y p_z')
zp, b = symbols('zp b')
kp_2 = symbols('k^2_p')
Rp = symbols('R_p')

R = sqrt(x**2 + y**2)
phi = atan(y/x)

term_px_py = px*cos(phi) + py*sin(phi)

#term1 = (Abs(z)**2/((Abs(z)**2 + R**2)**(3/2)) - (Abs(z)**2 + R**2)**(-1/2))/R
#term2 = Abs(z)/((Abs(z)**2 + R**2)**(3/2))

expo_f = exp(-kp*(2*zp -z))

factor_comun = I*pi*kp_2*Rp*expo_f

H0 = hankel1(0, kp*R)*factor_comun
H1 = hankel1(1, kp*R)*factor_comun

#term1_final = -term_px_py*(term1 + H1)
#term2_final = -pz*sign(z)*(term2 + H0)

term1_final = -term_px_py*H1
term2_final = -pz*sign(z)*H0

term_tot_final = term1_final + term2_final

Efield_x = -diff( term_tot_final, x)
Efield_y = -diff( term_tot_final, y)
Efield_z = -diff( term_tot_final, z)
    
#print_latex(Efield_x)
print_latex(Efield_y)


# Efield_x2 =  Efield_x.subs(sqrt((x+xD)**2 + (y+yD)**2),R2)
# print_latex(Efield_x)

# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

