
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: leila

graficar el campo externo directo con la convencion de z hacia abajo
en z = 0
graficar mapa de color x,y
"""

import numpy as np
import os 
import matplotlib.pyplot as plt



#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'background'
if not os.path.exists(path_save):
    print('Creating folder to save graphs')
    os.mkdir(path_save)


sigma = 0.2
yx_ratio = 0.5

theta1 = 30*np.pi/180
theta2 = 60*np.pi/180
theta3 = 90*np.pi/180
theta4  = 0

def background_theta(x,y):
    
    def gaussian1(theta,sigma):
        exponente = (x*np.cos(theta) - y*yx_ratio*np.sin(theta) )**2/(sigma**2)
        return np.exp(-exponente)


    def gaussian2(theta,sigma):
        exponente = (x*np.cos(theta) + y*yx_ratio*np.sin(theta))**2/(sigma**2)
        return np.exp(-exponente)

    def gaussian3(sigma):
        exponente = (x**2 + y**2)/(sigma**2)
        return np.exp(-exponente)

    def gaussian4(sigma):
        exponente = (x + y)**2/(sigma**2)
        return np.exp(-exponente)
    
   
    rta1 = gaussian1(theta1,sigma) + gaussian1(theta2,sigma) + gaussian1(theta3,sigma)
    rta2 = gaussian2(theta1,sigma) + gaussian2(theta2,sigma)  
    rta3 = gaussian3(sigma)  

    amplitud = np.sqrt(x**2 + (y/yx_ratio)**2)
    
    amplitud2 = np.sqrt(x**2 + y**2)

#    a = np.sqrt(0.7)
#    b = np.sqrt(0.95)
#
#    if (x/a)**2 + (y/b)**2 <= 1 :
#        cte = 0.5
#    else:
#        cte = 0
    
    rta_final = (rta1 + rta2)

        
    return rta_final
    


#%%
    
    
tamfig = [2.75, 4]
tamletra = 7
tamtitle  = 8
tamnum = 6
tamlegend = 6
labelpady = 2
labelpadx = 3
pad = 2.5
mk = 1
ms = 1
hp = 0.5
length_marker = 1.5
dpi = 500



#%% 
N = 1e3
listx = np.linspace(-1.5,1.5,N)
listy = np.linspace(-2.5,2.5,N)

X, Y = np.meshgrid(listx, listy)
f_num = np.vectorize(background_theta)
Z_num = f_num(X, Y)

#%% 
import matplotlib.colors as mcolors


### primera paleta https://coolors.co/gradient-palette/e6e6e6-ffdb3a?number=5
paleta1 = ["#E6E6E6","#E6E6E6","#ECE3BB","#F3E190","#F9DE65","#FFDB3A"]


           ## agregar mas amarillo a la paleta 1 ######## y repetir el gris que es el primero 
           ## https://coolors.co/gradient-palette/e6e6e6-fcce04?number=7
           
paleta1 = [ "#E6E6E6", "#E6E6E6", "#EAE2C0", "#EDDE9B", "#F1DA75", "#F5D64F", "#F8D22A" ]
                  #E6E6E6, #EAE2C0, #EDDE9B, #F1DA75, #F5D64F, #F8D22A, #FCCE04
                  
###  segunda paleta https://coolors.co/gradient-palette/fcce04-f75c14?number=8
paleta2 = ["#FCCE04","#FBC106","#FBB509","#FAA80B","#F99C0D","#F88F0F","#F88312","#F77614"]

colors = paleta1 + paleta2

cmap = mcolors.ListedColormap(colors)

norm = mcolors.Normalize(vmin = np.min(Z_num) ,vmax = np.max(Z_num) )


limits = [np.min(listx) , np.max(listx), np.min(listy) , np.max(listy)]
plt.figure(figsize=tamfig,frameon=False)
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    right=False,      # ticks along the bottom edge are off
    left=False,         # ticks along the top edge are off
    labelleft=False)


#for x in [x1,x2,x3,x4]:
#    plt.plot(ejey_aux1,x*np.ones(10),'--',color = 'grey' )
#norm = mcolors.DivergingNorm(vmin=np.min(Z_num), vmax = np.max(Z_num),vcenter = 0.5)
im = plt.imshow(Z_num, extent = limits, cmap=cmap, aspect='auto', interpolation = 'none',origin = 'lower',norm=norm) 
#plt.plot(maxis,0.18*np.ones(len(maxis)),'o',color = 'green' )
#cbar = plt.colorbar(im, fraction=0.046, pad=0.04, orientation = 'vertical')
#plt.clim(np.min(ticks_z),np.max(ticks_z))
plt.tight_layout()
os.chdir(path_save)
plt.box(False)
plt.savefig( 'background_plane'  + '.png', format='png',bbox_inches='tight',pad_inches = 0,dpi = dpi, transparent=True)   

#%%         
        












    
