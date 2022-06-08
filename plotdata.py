# Template for nice plot 
# 
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt 
import scipy.special
import datetime
#
# Input and organize data
# store the data in data

# files to import. 
files=['r_w.dat', 'fft_exact.dat', 'fft_final.dat', 'simpleRK4.dat' ] 

#example data to plot  
data = [] 

for datfile in files:
    dat = np.genfromtxt(datfile, delimiter=",")
    # print(csv) 
    data.append( dat )   




###########plottting................
mycolors = [ 'k', 'r', 'b', 'm', 'g' ] 
# parameters. 
ncol=1
nrow=2
fignum = 1
tlen = 15.0 # length of tick marks  
twid = 1.5 # width of tick marks
xtick = [] 
ytics = [] 
# limits. 
xlims = [] 
ylims = [] 
# labels 
# xlabs = [] 
# ylabs = [] 

xlabs = [ [], '$\omega$'  ]   # blank [] means it plots without a label.  
ylabs = [ 'Re $y$', 'Im $y$'   ] 
plttitle = datetime.date.today().strftime("%b %d, %Y") \
        + " Fourier Transform of $\exp\left[ -(t-5)^2/2  \\right] $" 
subplttitle=[]
legstr = [ [], [] ] 
fig = plt.figure(num=fignum, figsize=(6,8.5) )
mngr = plt.get_current_fig_manager()  # size on screen?
mngr.window.wm_geometry( "+000+400" ) 

# help string for labels 
strhlp = [ 'Re', 'Im' ] 

# 
for col in range(ncol):
    for row in range(nrow): 
        ax = plt.subplot( nrow, ncol, col+ (ncol)*(row)+1 )
        # modify the axes. 
        ax.xaxis.set_tick_params(length=tlen,width=twid)
        ax.yaxis.set_tick_params(length=tlen,width=twid) 
        p=[]  # this holds plots for legend. 
        # specify each row.  
        #     tmp,=plt.plot( data[3][:,0], data[3][:,1], color='k'  ) 
        #     p.append(tmp) 
        #     legstr[row].append('numeric' )

        tmp,=plt.plot( data[0], data[1][:,row], color='r', linewidth=2.0  ) 
        p.append(tmp) 
        legstr[row].append('exact {}'.format(strhlp[row]) )

        tmp,=plt.plot( data[0], data[2][:,row], color='k', linewidth=1.0  ) 
        p.append(tmp) 
        legstr[row].append('numeric  {}'.format(strhlp[row]))



        # tick marks. 

        # limits 
        if xlims:  # put on the limits if not empty
            if xlims[row]: 
                plt.xlim( xlims[row] )
        if ylims:
            if ylims[row]:
                plt.ylim( ylims[row] )
        # labels 
        if xlabs:
            if xlabs[row]:
                plt.xlabel( xlabs[row] )
        if ylabs:
            if ylabs[row]:
                plt.ylabel( ylabs[row] )
        if legstr:
            if legstr[row]:
                plt.legend(p,legstr[row],frameon=False, fontsize=8)
        # # overall plot title. 
        if row==0 and col==0: 
            plt.title( plttitle )

        # ax.set_xscale('log') 
        # ax.set_yscale('log') 



### plt.ion()
plt.tight_layout()
plt.show()
# time.sleep(2)
# plt.savefig(  )  # for saving the figure  
# plt.close()













