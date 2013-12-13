#! /usr/bin/python2 
# modules
import matplotlib.pyplot as plt;
import numpy as np;
# parameters 
# ---- your code here -----
f=open('out/g_2.dat','r')
data=np.array([[jj for jj in ii.split()]for ii in f],dtype='float')
f.close()
f2=open('out/eta','r')
eta=np.array([[jj for jj in ii.split()]for ii in f2],dtype='float')[0,0]
f2.close()
fr=open('ref/g2_ref.dat','r')
ref=np.array([[jj for jj in ii.split()]for ii in fr],dtype='float')
fr.close()
# calculate error
print "the relative error is : "+str(np.sum(np.abs(data[:,1]-ref[:,1]))/np.sum(np.abs(ref[:,1])))
# ---- your code here -----
# plot
plt.grid();
plt.plot(data[:,0],data[:,1],lw=3);
plt.plot(ref[:,0],ref[:,1],'x');

# plot settings
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=16)
plt.title(r'RDF of Hard Sphere at $\eta$ = '+str(eta),fontsize=24);
plt.xlabel(r"$\frac{r}{D}$",fontsize=24)
plt.ylabel(r"$g_2(r)$",fontsize=30)
#  plt.ylim([0,np.max(Diag_u_k)]);
#  plt.xlim([-pi/a,pi/a]);
plt.legend(['MC','Exact PY']);
plt.tight_layout();
plt.savefig("g2.png", transparent = True,dpi=300)
#  plt.show()
