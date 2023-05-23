# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 09:58:01 2022

@author: dymet
"""
###############################################################################
###############################################################################
############### Quartic Interpolation of Light Curve Data  ####################
############### By: Dymetris Ramirez and Colter Richardson ####################
###############################################################################
###############################################################################

### Quartic Fit to SN Light Curve Data with Unknown Shock Breakout Time ###
def QuarticLightCurveFit(SN_dat, left_bound, right_bound, mass_known, unc, left_X, right_X):
    
    '''
    *Make sure that the data file and this python file are saved to the same folder
    
    Definition of Inputs:
    SN_data: The name of the text file or dat file inputed into the function. Example - "SN_Name.dat" or "SN_Name.txt"
    left_bound: Farthest left data you are wanting to fit to in the data set
    right_bound: Farthest right data you are wanting to fit to in the data set 
    mass_known: If the mass of the progenitor is known type in 0, if it is not type in 1
    unc: The value taken from the mass graph that gives an uncertainty value pertaining to the mass of the progenitor (See README)
    left_X and right_X: where you want the interpolation line to start and end in your fitted curve
    '''
    
    ### Import of Packages Used
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.optimize
    import pandas as pd
    
    
    ### Loading in Data 
    data = pd.read_csv(str(SN_dat),sep="\s+", names=['time', 'mag', 'dmag','flux','dflux']) #<=== When importing the data, if the txt, csv, or dat file are missing a file
                                                                                           #just remove the column title from the name section, if there is more, simply just add the title to the names
    
    # Shift time data for the first point to start at zero
    time = data['time'].values
    time = time - time.min()
    
    flux = data['flux'].values    
    
    ### Interpolation of SN
    # Setting the interval for the Quartic Interpolation
    x = time[int(left_bound):int(right_bound)]
    y = flux[int(left_bound):int(right_bound)]
        
    ### Functions Used to Do Chi^2 Anaylsis
    #The purpose of the Chi^2 is to minimize the parameters of the Quartic Function to best fit the data set 
    def chisqr(ymeasure, yfit, error):
        return np.sum((ymeasure - yfit)**2/(error)**2)

    def calcunc(ymeasure, yfit, DOF):
        return np.sqrt(((1 / (len(yfit) - DOF)) * np.sum((yfit - ymeasure)**2)))

    # Fitted Equation (Quartic)
    def line(x, a, b, c, d, e):
        return a*x**4 + b*x**3 + c*x**2 + d*x + e 
    
    ### Fitting The Data
    # Extracting the parameters from curve fitting the data of time and flux to the quartic fit
    params, _ = scipy.optimize.curve_fit(line, x, y)

    # Extracting parameters for the Function
    a, b, c, d, e = params
    
    # Generating a function of the fitted data with the derived parameters
    yfit = line(x, a, b, c, d, e)

    # Calculating error from data
    dof = len(y) - 1
    error = calcunc(y, yfit, dof) * np.ones(len(x))

    # Chi squared calculation function
    x_2 = chisqr(y, yfit, error)
    #print('Chi Squared = ',(x_2))

    # Reduced chi squared calculation
    #x2r = 1 / (len(y) - 2) * x_2
    #print('Reduced Chi Squared = ', (x2r))
    print('a = %.5f, b = %.5f, c = %.5f, d = %.5f, e = %.5f' % (a, b, c, d, e))

    # Range of data for fitted function inputted into fitted function
    X1 = np.arange(int(left_X),int(right_X),.0001)  
    Y1 = a*X1**4 + b*X1**3 + c*X1**2 + d*X1 + e
    
    # Plotting Single Fit of Data
    plt.plot(X1,Y1/np.max(flux),'--', color = 'red', label = 'Quartic Interpolation' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    plt.plot(time,flux/np.max(flux), '.', color = 'black', label = 'Data for ' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    
    # Finding where the quartic fit intersects the x-axis which will serve as our estimated point of shock breakout
    index = np.argmin(np.abs(Y1))
    print(index)
    mu = X1[index]
    
        
   # Estimated time of Collapse and OSW
    tcollapse = mu-int(unc)
    print('Time of Shock Breakout = ', np.abs(mu))
    print('time of collapse estimate in days = ', np.abs(tcollapse))

    ### Computing for t1 and t2 values for the OSW
    # Sigma Value taken from KEPLER QUARTIC python code
    #sigma = .25 #(While contiously looking at the vertical error of Kepler for half a day gives the updated sigma of .25)
    #delta_x = 0.2 additional uncertainty looking at Sean Couch graph of shock breakout time versue mass of progenitor
    #if mass_known == 0:
    #    v1 = mu -int(unc)-(sigma*2) - delta_x
    #    v2 = mu-int(unc)+(sigma*2) + delta_x
    
    #if mass_known == 1:
    #    v1 = mu-int(unc)-3 #We subtracted 3 days additional to account for the relationship of the shock breakout delay as a function of progenitor mass. Since our mass can range from 8-20 solar masses the range can extended an additional 3 days then the original 1.5
    #    v2 = mu -int(unc) + (sigma*2) + delta_x
    
    #print('t1 = ', np.abs(v1))
    #print('t2 = ', np.abs(v2))
    #print('OSW value in days = ', np.abs((v1 - v2)))
    
    ### Plotting the Data to Show the OSW on the Graph
    #plt.axvline(tcollapse, ymax = .2, color = 'blue', label = 'Time of Collapse')
    #plt.axvline(v1, ymax = .2, color = 'orange', label = 't1 value')
    #plt.axvline(v2, ymax = .2, color = 'orange', label = 't2 value')
    # # Prints the values of the OSW on the graph, varies depending on where your OSW values are (this is not neccessary just a to make the graph look clean)
    # # plt.text(-9.9,.00285,"t1 =" + str(np.round(v1,5)))
    # # plt.text(-9.9,.00255,"t2 =" + str(np.round(v2,5)))
    plt.ylabel('Normalized Flux')
    plt.xlabel('Days')
    plt.title(str(SN_dat))
    plt.minorticks_on()
    # X and Y limits vary depending on the light curve
    plt.ylim(0, 1.1)
    plt.xlim(-12,30)
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(fontsize = "xx-small", loc = 4)
    plt.savefig(str(SN_dat) + '.png', dpi = 300, bbox_inches = "tight")
    plt.show()
    plt.close()

### Using the Quartic Package on Supernovae Data:
#The interval for all the quartic fits below are for the first 25 days because the sigma value we are getting from Kepler is the fit of up to 25 days. If you wish to change your interval for the supernovae data, you will have to increase your max interval in Kepler to the same amount to get a corresponding sigma value
QuarticLightCurveFit('SN2019ehk_r.dat', 0,18,0, 1.5,-40,30) 
QuarticLightCurveFit('SN2019gaf_r.txt', 0,18,1,1.5,-40,30)
QuarticLightCurveFit('SN2020oi_r.dat', 0,34,1,1.5,-40,30)


    
    
    
    
    
