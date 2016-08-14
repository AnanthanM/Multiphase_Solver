#!/usr/bin/python

import pylab as pl
import numpy as np
pl.rc('legend',fontsize=12)
#pl.rc('text',usetex=True);
#pl.rc('font',family='serif');

data_ghia_u=np.loadtxt('re1000u.txt')
data_ghia_v=np.loadtxt('re1000v.txt')
data_u=np.loadtxt('line_data_u.dat')
data_v=np.loadtxt('line_data_v.dat')

x_ghia = data_ghia_v[:,0]
v_ghia = data_ghia_v[:,1]
y_ghia = data_ghia_u[:,1]
u_ghia = data_ghia_u[:,0]

x_fvm = data_v[:,0]
v_fvm = data_v[:,1]
y_fvm = data_u[:,0]
u_fvm = data_u[:,1]

pl.figure(1,figsize=(8,4))
pl.subplot(1,2,1)
pl.plot(u_ghia,y_ghia,'k^',linewidth=1.5,label='Ghia et. al.')
pl.plot(u_fvm,y_fvm,'r-',linewidth=1.5,label=r'FVM $(\Delta x = 1/140)$')
lg=pl.legend(loc=4, handlelength=2)
lg.draw_frame(False)
#pl.show
pl.xlabel('U ',fontsize=12)
pl.ylabel('y',fontsize=12)
pl.xticks(fontsize=12)
pl.yticks(fontsize=12)
pl.xlim(-0.5,1.0)
pl.ylim(0.0,1.0)
#pl.xticks(0.1)
pl.subplot(1,2,2)
pl.plot(x_ghia, v_ghia, 'k^',linewidth=1.5,label='Ghia et. al.')
pl.plot(x_fvm, v_fvm,'r-',linewidth=1.5,label=r'FVM $(\Delta x =1/140)$')
lg=pl.legend(loc=3,handlelength=2)
lg.draw_frame(False)
#pl.show
pl.xlabel('x ',fontsize=12)
pl.ylabel('V',fontsize=12)
pl.xticks(fontsize=12)
pl.yticks(fontsize=12)
pl.xlim(0.0,1.0)
pl.ylim(-0.6,0.4)


pl.tight_layout()
pl.savefig("uv_1000.pdf")



pl.show()
