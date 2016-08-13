#!/usr/bin/python

import pylab as pl
import numpy as np
pl.rc('legend',fontsize=12)
#pl.rc('text',usetex=True);
#pl.rc('font',family='serif');
fig, ax1 = pl.subplots(1)

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

#pl.figure(1,figsize=(8,4))
#pl.subplot(1,2,1)
ax1.plot(u_ghia,y_ghia,'k^',linewidth=1.5,label='Ghia et. al.')
ax1.set_xlabel('U ',fontsize=12, color='b')
ax1.set_ylabel('y',fontsize=12, color='b')
for t1 in ax1.get_yticklabels():
    t1.set_color('b')
for t1 in ax1.get_xticklabels():
    t1.set_color('b')

ax2 = ax1.twinx()
#ax2 = ax1.twiny()
ax2 = fig.add_axes(frameon=False)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')

ax2.plot(x_ghia, v_ghia,'k^',linewidth=1.5,label='Ghia et. al.')
ax2.set_ylabel('V',fontsize=12, color='r')
ax2.set_xlabel('x',fontsize=12, color='r')
for t1 in ax2.get_yticklabels():
    t1.set_color('r')
for t1 in ax2.get_xticklabels():
    t1.set_color('r')

pl.show()

