#! /usr/bin/python2 
# modules
import matplotlib.pyplot as plt;
import numpy as np;
# parameters 
# ---- your code here -----
f=open('out/r_av','r')
data=np.array([[jj for jj in ii.split()]for ii in f],dtype='float')
f.close()
f2=open('out/eta','r')
eta=np.array([[jj for jj in ii.split()]for ii in f2],dtype='float')[0,0]
f2.close()
# ---- your code here -----
# plot
plt.grid();
plt.plot(data[:,0]/1000,data[:,1],lw=3);
plt.plot(data[:,0]/1000,data[:,2],lw=3);
plt.plot(data[:,0]/1000,data[:,3],lw=3);

# plot settings
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=16)
plt.title(r'Equilibration: Mean Position at $\eta$ = '+str(eta),fontsize=24);
plt.xlabel(r'$\frac{MC\,Steps}{1000}$',fontsize=24)
plt.ylabel(r"$\bar{r}_i$",fontsize=30)
#  plt.ylim([0,np.max(Diag_u_k)]);
#  plt.xlim([-pi/a,pi/a]);
plt.legend([r"$\bar{x}$",r"$\bar{y}$",r"$\bar{z}$"],loc=4);
plt.tight_layout();
#  plt.savefig("eq.png", transparent = True,dpi=300)
plt.show()
