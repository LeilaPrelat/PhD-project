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
kx = symbols('kx')
ky = symbols('ky')
px, py, pz = symbols('p_x p_y p_z')
zp, b = symbols('zp b')
Rp = symbols('R_p')
kparallel = symbols('k_\parallel')

aux = sqrt(kx**2 + ky**2)
exp1 = exp(-aux*Abs(z))
exp2 = exp(-2*aux*zp)
rp = Rp*exp2*kp/(aux - kp)

term1 = -1j*kx*px/aux - 1j*ky*py/aux + pz*sign(z)*ky
term2 = 1j*px*kx*rp/aux + 1j*py*ky*rp/aux 
term3 = -pz*sign(z)*ky*rp

exp_x = exp(1j*kx*x)
exp_y = exp(1j*ky*y)


term_tot_final = (term1 + term2 + term3)*exp_x*exp_y*exp1 

Efield_x = -diff( term_tot_final, x)
Efield_y = -diff( term_tot_final, y)
Efield_z = -diff( term_tot_final, z)
    
term_aux = exp(-kparallel*Abs(z))
term_aux_f = -diff(term_aux,z)

print_latex(term_aux_f)


# Efield_x2 =  Efield_x.subs(sqrt((x+xD)**2 + (y+yD)**2),R2)
# print_latex(Efield_x)

# para el potential 
#print_latex(term_tot_diff_x)
#print('')
#print_latex(term_tot_diff_z)

