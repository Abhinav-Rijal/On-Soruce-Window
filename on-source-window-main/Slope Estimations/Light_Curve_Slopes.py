# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 20:47:27 2022

@author: dymet
"""

#ASASSN-18qk_TESS


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import pandas as pd

def line(x,a,b,c,d,e):
    y = a * x**4 + b * x**3 + c * x**2 + d * x + e
    return y

Slopes = []

#FLUX IS ACTUALLY MAG
#Load in SN Data and Fit Slope
Name = 'ASASSN-18qk'
df1 = np.genfromtxt("ASASSN-18qk_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]

x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2   
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ASASSN-18qk_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()


Slopes = np.append(Slopes,params2[3])

Name = 'ASASSN-18qv'
df1 = np.genfromtxt("ASASSN-18qv_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2 
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ASASSN-18qv_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])

Name = 'ASASSN-19jy'
df1 = np.genfromtxt("ASASSN-19jy_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2   
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ASASSN-19jy_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])

Name = 'ASASSN-19or'
df1 = np.genfromtxt("ASASSN-19or_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2 
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ASASSN-19or_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])

Name = 'ZTF18abzscns'
df1 = np.genfromtxt("ZTF18abzscns_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2  
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ZTF18abzscns_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])

Name = 'ZTF19abqhobb'
df1 = np.genfromtxt("ZTF19abqhobb_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2 
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.xlim(-1,20)
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ZTF19abqhobb_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])

Name = 'ZTF20aahbamv'
df1 = np.genfromtxt("ZTF20aahbamv_TESS.txt")
time = df1[1:,0]
time = time - time.min()
flux = df1[1:,3]
x1 = time[0:]
y1 = flux[0:]

index = np.argmin(y1)
print(index)

x = time[0:index]
y = flux[0:index]


#######PLOTTING
params2, _ = scipy.optimize.curve_fit(line, x, y)
a2, b2, c2, d2, e2 = params2
yfit2 = line(x, a2, b2, c2, d2, e2)
X = np.arange(-1,x1[index],.1)
Y = a2*X**4 + b2*X**3 + c2*X**2 + d2*X + e2 
plt.figure()
plt.plot(x1,y1, '.', label = 'Data Points')
plt.plot(X,Y, label = 'Fitted Quartic Function')
plt.title(Name)
plt.grid()
plt.legend()
plt.xlim(-1,20)
plt.ylabel('Magnitude')
plt.xlabel('Time (Days)')
plt.savefig('ZTF20aahbamv_TESS.png', dpi = 300, bbox_inches = "tight")
plt.gca().invert_yaxis()
plt.show()
plt.close()

Slopes = np.append(Slopes,params2[3])
Slopes = (Slopes) *-1 #Since Magnitude measures inversely 

#Slopes from Plot on Wiki
Slope = np.array([.18,.33,.19,.23,.29,.19])*1.87
Slopes = np.append(Slopes,Slope)

#Adding Kepler Slopes (+/- 10%)
KSlope = np.array([.344])*1.87
Slopes = np.append(Slopes,KSlope)

#Histogram of Slopes:
plt.figure()
plt.hist(Slopes,bins = len(Slopes))
mu = np.mean(Slopes)
sigma = np.std(Slopes)
plt.vlines(mu,0,3, linewidth = 3, color = "Green", label = r"$\mu=$" + str(np.round(mu,4)))
plt.vlines(mu+sigma,0,(2.5), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)))
plt.vlines(mu-sigma,0,(2.5), color = "Red")
plt.vlines(.344*1.87,0,3, color = 'purple', label = 'Kepler Estimated Slope')
plt.vlines(.5676*1.87,0,(2.5), color = 'Black', label = 'Kepler Plus 65% Error')
plt.vlines(.1204*1.87,0,(2.5), color = 'Black', label = 'Kepler Minus 65% Error')
#plt.vlines(.3784*1.87,0,(2.5), color = 'orange', label = 'Kepler Plus 10% Error')
#plt.vlines(.3096*1.87,0,(2.5), color = 'orange', label = 'Kepler Minus 10% Error')
plt.xlabel("Slope of Rise Time")
plt.grid(which = "both")
plt.minorticks_on()
plt.title('Histogram of Rise Time Slopes of Different SNe')
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(fontsize = "small")
plt.savefig('Slopes_HIST_2.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()  

#Scale = 1.87


