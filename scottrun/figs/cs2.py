import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.13,0.8,0.8])
#y00=-100
#y01=40
#y02=-5
#a00=0.75
#a01=1.5
#a02=4.5

colors=['black','red','green','blue','cyan','violet']

mydata = np.loadtxt('../results_noint/T100.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=0
plt.plot(rho,cs2,linestyle='--',linewidth=2,color=colors[iplot],marker=None,label='no int., T=100')

mydata = np.loadtxt('../results_noint/T120.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=0
plt.plot(rho,cs2,linestyle='dotted',linewidth=2,color=colors[iplot],marker=None,label='no int. T=120')

mydata = np.loadtxt('../results_noint/T140.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=0
plt.plot(rho,cs2,linestyle='-',linewidth=2,color=colors[iplot],marker=None,label='no int. T=140')
   

mydata = np.loadtxt('../results_eos1/T100.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=1
plt.plot(rho,cs2,linestyle='--',linewidth=2,color=colors[iplot],marker=None,label='+int., T=100')

mydata = np.loadtxt('../results_eos1/T120.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=1
plt.plot(rho,cs2,linestyle='dotted',linewidth=2,color=colors[iplot],marker=None,label='+int., T=120')

mydata = np.loadtxt('../results_eos1/T140.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=1
plt.plot(rho,cs2,linestyle='-',linewidth=2,color=colors[iplot],marker=None,label='+int., T=140')

mydata = np.loadtxt('../results_eos2/T100.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=2
plt.plot(rho,cs2,linestyle='--',linewidth=2,color=colors[iplot],marker=None,label='+int., T=100')

mydata = np.loadtxt('../results_eos2/T120.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=2
plt.plot(rho,cs2,linestyle='dotted',linewidth=2,color=colors[iplot],marker=None,label='+int., T=120')

mydata = np.loadtxt('../results_eos2/T140.txt',skiprows=1,unpack=True)
rho=mydata[0]
cs2=mydata[1]
iplot=2
plt.plot(rho,cs2,linestyle='-',linewidth=2,color=colors[iplot],marker=None,label='+int., T=140')

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,1.2,0.1), minor=False)
ax.set_xticklabels(np.arange(0,1.2,0.1), minor=False, family='serif')
ax.set_xticks(np.arange(0,1.2,0.05), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,1.0)

ax.set_yticks(np.arange(0,1.0,0.25), minor=False)
ax.set_yticklabels(np.arange(0,1.0,0.25), minor=False, family='serif')
ax.set_yticks(np.arange(0,1.0,0.05), minor=True)
plt.ylim(0,0.7)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\\rho$ [fm$^{-3}$]', fontsize=18, weight='normal')
plt.ylabel('$c_s^2$ ',fontsize=18,labelpad=0)

legend(loc="upper left")

#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('cs2.pdf',format='pdf')
os.system('open -a Preview cs2.pdf')
#plt.show()
quit()
