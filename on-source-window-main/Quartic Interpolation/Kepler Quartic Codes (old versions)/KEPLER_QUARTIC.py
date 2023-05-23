# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 09:37:42 2022

@author: dymet
"""

### Quartic Fit to SN Light Curve Data with Known Shock Breakout Time ###
### By: Dymetris Ramirez and Colter Richardson ###
# This is used to derive the sigma value which will be used in the QuarticLightCurveFit package below (See README)

### Loading in Python Packages Used Throughout the Code
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from astropy.timeseries import TimeSeries
from scipy import optimize
import time as timer
Start = timer.time()

### Data Import
ts1 = TimeSeries.read('kplr008480662-2011271113734_llc.fits', format='kepler.fits')
time = np.array(ts1.time.jd)
sap_flux = np.array(ts1['sap_flux'])

### Removing Any Data That is Found as 'nan'. Since This Data is Irellevant in the Fit 
index = np.argwhere(np.isnan(sap_flux))
sap_flux = np.delete(sap_flux,index)
time = np.delete(time,index)

### Range of Model_time for How Long the Fitted Curve Extends in the Graph
model_time = np.arange(-20,15,0.001)

### Shift Light Curve to 0 So That the Zero Point (0,0) is What Will be Used to Define Time of Shock Breakout
initial_time = ts1.time.jd[1320]

### Range of Ramp-Up in Data
# Rough estimate of breakout and peak indecies (the interval for which you will use the data to fit the quartic too)
shock_breakout = 1575 #Half Day Expolartion starts at 1350 #7 Day Exploration starts at 1575
peak_lum = 2500

# Indexing of data
fitting_points = sap_flux[shock_breakout : peak_lum]
fitting_time = ts1.time.jd[shock_breakout : peak_lum]

# Estimated by-eye slope of the ramp up of the light curve
d = 50

#Putting an initial condition for the Uncertainty in both time and flux
Ic = np.average(sap_flux[0:1700])
sap_flux_a = sap_flux - Ic
time_1 = time[1700:]
sap_flux_1 = sap_flux[1700:]
    
### Random sample in the ramp-up range
#Number or Random Values you want
N = 10
sample = random.randint(shock_breakout,peak_lum,N)
selected_elements_flux = [sap_flux[index] for index in sample]
selected_elements_time = [ts1.time.jd[index] for index in sample]

### Fitting With Quartic Equation
def model(x,a,b,c,d,e):
    y = a * x**4 + b * x**3 + c * x**2 + d * x + e
    return y

#We use the derivative of the Quartic Polynomial, as it allows us to restrain slopes within a percent error of the Quartic Interpolation
def model_derivative(x,a,b,c,d):
    y = 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d
    return y

# Extracting the parameters from curve fitting the data of selected_element_time and selected_elements_flux to the Quartic Model
par_opt, par_cov = optimize.curve_fit(model,selected_elements_time- initial_time, selected_elements_flux)

# Taking the minimum of the argument for each parameter at the set point of when shock breakout occurs
B = 4915 #<====== estimated sap_flux at the time of shock breakout 
index = np.argmin(np.abs(np.array(model(model_time, par_opt[[0]], par_opt[[1]], par_opt[[2]], par_opt[[3]], par_opt[[4]]))-B))

# Number of iterations you are wanting to run through and the amount of random points you want to choose
random_sample_number = 10000
random_sample = 10

# Creating an empty array to store the calculated shock_breakout_estimate in
shock_breakout_estimate = np.zeros(random_sample_number)

### Plot
fig, ax = plt.subplots(1,1,figsize=([10,10]))
ax.plot(time - initial_time, sap_flux, '.', markersize = 8, color = 'black', label = 'Data Points') 
#ax.plot(time - initial_time, sap_flux)

# For loop to run all the number of iterations
for i in range(random_sample_number):
    sample = random.randint(shock_breakout,peak_lum,random_sample)
    selected_elements_time = [time[index] for index in sample]
    selected_elements_flux = [sap_flux[index] for index in sample]

    # Extracting the parameters from curve fitting the data of selected_element_time and selected_elements_flux to the Quartic Model
    par_opt, par_cov = optimize.curve_fit(model,selected_elements_time - initial_time, selected_elements_flux, bounds = (-np.inf, [0, np.inf, np.inf, np.inf, np.inf]))
    
    # Plotting Data fits that have a derived slope between plus/minus 10 percent of the by-eye slope value of the interval
    if d - d*(10/100)<= model_derivative(model_time[[20000]],par_opt[[0]],par_opt[[1]],par_opt[[2]],par_opt[[3]]) <= d + d*(10/100): #<=== The value inside the model_time is the time point where you are getting your slope value (d) from
        ax.plot(model_time, model(model_time, par_opt[[0]], par_opt[[1]], par_opt[[2]], par_opt[[3]], par_opt[[4]]), alpha = 0.5)
        index = np.argmin(np.abs(np.array(model(model_time, par_opt[[0]], par_opt[[1]], par_opt[[2]], par_opt[[3]], par_opt[[4]]))-B))
        shock_breakout_estimate[i] = model_time[index]
    else:
        pass
        
ax.set_xlabel('Julian Date (Days)')
ax.set_ylabel('SAP Flux (e-/s)')
ax.set_ylim(4900,5120)
ax.set_title('Kepler Light Curve')
ax.grid(which = "both")
ax.minorticks_on()
ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
ax.legend(loc = "lower right")
plt.savefig('Kepler_Quartic.png', dpi = 300, bbox_inches = "tight")

# Delete the value of shock_breakout_estimates that produced zero as there value
shock_breakout_estimate = np.delete(shock_breakout_estimate, np.where(shock_breakout_estimate == 0))

# Plotting a histogram of the derived shock_breakout_estimate values
plt.figure()
y, x, _ = plt.hist(shock_breakout_estimate)
mu = np.average(shock_breakout_estimate) #<==== Average value of shock_breakout_estimate
sigma = np.std(shock_breakout_estimate) #<==== Standard deviation of the shock_breakout_estimate
plt.vlines(mu,0,y.max(), linewidth = 3, color = "Green", label = r"$\mu=$"+ str(np.round(mu,2)) + " days")
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,2)) + " days")
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(1,1,1, color = "White") 
plt.xlabel("Days from Shock Break-Out")
plt.grid(which = "both")
plt.minorticks_on()
plt.title('Kepler Histogram')
plt.xlim(-1.5,1)
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.legend(fontsize = "small")
plt.savefig('Kepler_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

#This sigma value is always changing since the code is choosing random points in between the interval of 25 days
print('The Sigma Value of Kepler is: ', sigma)
End = timer.time()
# Printing the time it took to run this code
print(End - Start, "seconds.")
