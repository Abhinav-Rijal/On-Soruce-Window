# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 09:58:01 2022

@author: dymet
"""
###############################################################################
###############################################################################
############### Quartic Interpolation of Light Curve Data  ####################
######################### By: Dymetris Ramirez ################################
###############################################################################
###############################################################################

### Quartic Fit to SN Light Curve Data with Unknown Shock Breakout Time ###
def QuarticLightCurveFit(SN_dat):
    '''
    *Make sure that the data file and this python file are saved to the same folder
    
    Definition of Inputs:
    SN_data: The name of the text file or dat file inputed into the function. Example - "SN_Name.dat" or "SN_Name.txt"
    Back_Data: The name of the text file or dat file inputed into the function. Example - "SN_Name.dat" or "SN_Name.txt"
    '''
    
    ### Import of Packages Used
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.optimize
    import pandas as pd
    
    
    #Background Data
    Average_Mag_r = 20.6375 #SN2019gaf 20.6375 #SN2019ehk 20.481
    ### Loading in Data 
    data = pd.read_csv(str(SN_dat),sep="\s+", names=['time', 'mag', 'dmag', 'band']) #['time', 'date', 'hour time','band','mag', 'unc sign', 'dmag', 'status', 'images1', 'image2', 'images3']) #<=== When importing the data, if the txt, csv, or dat file are missing a file
                                                                                            #just remove the column title from the name section, if there is more, simply just add the title to the names
    
    #print(data)
    # Shift time data for the first point to start at zero
    time = data['time'].values
    mag = data['mag'].values  
    Dmag = data['dmag'].values
    band = data['band'].values
        
    Time_r = time - time[0]
    Mag_r = mag
                
    ### Interpolation of SN
    # Setting the interval for the Quartic Interpolation
    x = Time_r[0:11] #int(Time_r[np.argmin(Mag_r)])]  #SN2019ehk = [0:13] #SN2019gaf = [0:]
    #print(x)
    y = Mag_r[0:11]#int(Time_r[np.argmin(Mag_r)])]
    #print(y)
    Dmag = Dmag[0:11]
    # Fitted Equation (Quartic)
    def line(x, a, b, c, d, e):
        return a*x**4 + b*x**3 + c*x**2 + d*x + e 
    
    ### Fitting The Data
    # Extracting the parameters from curve fitting the data of time and flux to the quartic fit
    params, _ = scipy.optimize.curve_fit(line, x, y)

    # Extracting parameters for the Function
    a, b, c, d, e = params

    # Range of data for fitted function inputted into fitted function
    X1 = np.arange(-80,12,.0001)  
    Y1 = a*X1**4 + b*X1**3 + c*X1**2 + d*X1 + e
    
    # Plotting Single Fit of Data
    plt.plot(X1,Y1,'--', color = 'red', label = 'Quartic Interpolation' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    plt.plot(x,y, '.', color = 'black', label = 'Data for ' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    plt.plot(Time_r,Mag_r, '.', color = 'black', label = 'Data for ' + str(SN_dat))
    
    # Finding where the quartic fit intersects the x-axis which will serve as our estimated point of shock breakout
    index = np.argmin(np.abs(Y1 - Average_Mag_r))
    #print(index)
    mu = X1[index]
    #print(mu)
    
        
   # Estimated Time of ShockBreakout
    print('Time of Shock Breakout = ', np.abs(mu))
    mu = np.round(mu,4)
    
    ### Plotting the Data to Show the OSW on the Graph
    plt.axvline(mu, ymax = .2, color = 'blue', label = 'tsbo value ' + str(mu) + ' Days')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.xlabel('Days')
    plt.title(str(SN_dat))
    plt.minorticks_on()
    # X and Y limits vary depending on the light curve
    plt.ylim(Average_Mag_r,y[-1]-3.25)
    plt.xlim(-10,40)
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(fontsize = "small", loc = 2)
    plt.savefig(str(SN_dat) + '.png', dpi = 300, bbox_inches = "tight")
    plt.show()
    plt.close()
    
    def derivative(x,a,b,c,d):
      return 4*a*x**3 + 3*b*x**2 + 2*c*x + d

    
    Tsbo = mu
    Trise = np.abs(derivative(mu,a,b,c,d)) #Tangential Slope
    Slope_Rise = np.abs((Y1[-1]-Y1[index])/(X1[-1]-X1[index]))
    Rise_Time = X1[-1]-X1[index]
    #print(Y1[index])
    
    selected_elements_time = x
    random_sample_number = 1000
    Range = 10
    shock_breakout_estimate = np.zeros(random_sample_number)
    sig = []
    Ind = np.argmin(y)
    d = (y[Ind]-Y1[index])/(x[Ind]-X1[index])
    
    #Limits
    SIG_2 = 2*0.4603
    Lower = [-0.7175 - SIG_2]
    Higher =[-0.7175 + SIG_2]
    #print(Lower)
    
    for i in range(Range):
        fig, ax = plt.subplots(1,1,figsize=([10,10]))
        plt.gca().invert_yaxis()
        ax.errorbar(x, y, yerr = Dmag, fmt = '.', color = 'black', label = 'Original Data Points with Uncertainty')
        Mag = np.zeros(((random_sample_number),0))
        for j in range(random_sample_number):
            for i in range(0,len(Mag_r)):
                hi = y-Dmag
                low = y+Dmag
            
            magR = np.random.uniform(low,hi)
            Mag = np.append(Mag,magR)
            #print(len(Mag))
            
        
        # Fitting Function 
        ax.legend()
        shock_breakout_estimate = []
        
        for j in range(0,len(Mag), len(y)):
            Mag2 = Mag[j:j+len(y)]
            ax.plot(x[0:len(y)],Mag2[0:len(y)], '.', color = 'red', label = 'Data Points picked from Uncertainty Bound')
            
        
            params2, _ = scipy.optimize.curve_fit(line, selected_elements_time[0:len(y)], Mag2[0:len(y)])
            a2, b2, c2, d2, e2 = params2
            print(params2[[3]])
            #print(d - d*(65/100))
            #print(Higher)
            
            if Lower <= params2[[3]] <= Higher:               
                Y2 = a2*X1**4 + b2*X1**3 + c2*X1**2 + d2*X1 + e2   
                ax.plot(X1,Y2)
            
                #Finding where the quartic fit intersects the x-axis which will serve as our estimated point of shock breakout
                index = np.argmin(np.abs(Y2 - Average_Mag_r))
                shock_breakout_estimate.append(X1[index])
        
        
        
        ax.set_xlabel('Days')
        ax.set_ylabel('Flux')
        ax.set_title(str(SN_dat) + ' Interpolation')
        ax.grid(which = "both")
        #ax.legend()
        plt.ylim(Average_Mag_r,y[-1]-2.25)
        plt.xlim(mu-1,40)
        ax.minorticks_on()
        ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.savefig(str(SN_dat) + '_ERROR.png', dpi = 300, bbox_inches = "tight")
        plt.show()
        plt.close()  
        #ax.legend(fontsize= 'medium', loc = "lower right")
                
            
        #Estimated Time of ShockBreakout
        print('Time of Shock Breakout = ', np.abs(mu))
        
        shock_breakout_estimate = np.delete(shock_breakout_estimate, np.where(shock_breakout_estimate == 0))
    
        
        plt.figure()
        plt.hist(shock_breakout_estimate)
        mu = np.average(shock_breakout_estimate) 
        sigma = np.std(shock_breakout_estimate)
        plt.vlines(mu,0,y.max(), linewidth = 3, color = "Green", label = r"$\mu=$"+ str(np.round(mu,7)) + " days")
        plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,7)) + " days")
        plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
        plt.vlines(1,1,1, color = "White") #label = str(len(shock_breakout_estimate)) + " fits, out of " + str(random_sample_number) + " tests")
        plt.xlabel("Days from Shock Break-Out")
        plt.grid(which = "both")
        plt.minorticks_on()
        plt.title('Histogram of Interpolations')
        plt.xlim(mu-.3,mu+.3)
        plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.legend(fontsize = "medium")
        plt.savefig(str(SN_dat) + '_HIST.png', dpi = 300, bbox_inches = "tight")
        plt.show()
        plt.close()  
            
        print('The estimated time of shock breakout: ', abs(int(x[0])-mu))
        
        sig = np.append(sig, sigma)
       
    plt.figure()
    y, x, _ = plt.hist(sig)
    mu = np.average(sig) 
    sigma = np.std(sig)
    plt.vlines(mu,0,y.max(), linewidth = 3, color = "Green", label = r"$\mu=$"+ str(np.round(mu,7)) + " days")
    plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,7)) + " days")
    plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
    plt.vlines(1,1,1, color = "White") #label = str(len(shock_breakout_estimate)) + " fits, out of " + str(random_sample_number) + " tests")
    plt.xlabel("Days from Shock Break-Out")
    plt.grid(which = "both")
    plt.minorticks_on()
    plt.title('Histogram of Sigmas')
    plt.xlim(mu-1,mu+1)
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(fontsize = "medium")
    plt.savefig(str(SN_dat) + '_Sigma_HIST.png', dpi = 300, bbox_inches = "tight")
    plt.show()
    plt.close()  
    
    return [Tsbo,Trise,Slope_Rise,Rise_Time,Average_Mag_r, Time_r]
    
#Storing Multiple Data Samples    
import numpy as np
import matplotlib.pyplot as plt

Tsbo = [] 
Trise = [] 
Slope_Rise = []
Rise_Time = []
Average_Mag_r = []

### Using the Quartic Package on Supernovae Data:
SN = QuarticLightCurveFit('SN2019gaf.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Time = SN[5]


