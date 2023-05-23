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
def QuarticLightCurveFit(SN_dat, Back_Data):
    '''
    *Make sure that the data file and this python file are saved to the same folder
    
    Definition of Inputs:
    SN_data: The name of the text file or dat file inputed into the function. Example - "SN_Name.dat" or "SN_Name.txt"
    Back_Data: The name of the text file or dat file inputed into the function. Example - "SN_Name.dat" or "SN_Name.txt"
    Range: Number of Ranges
    '''

    ### Import of Packages Used
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.optimize
    import pandas as pd
    
    
    #Background Data
    B_Data = pd.read_csv(str(Back_Data),sep="\s+", names=['time','band','mag', 'status']) #<=== When importing the data, if the txt, csv, or dat file are missing a file

    time = B_Data['time'].values
    #time = time - time.max()
    mag = B_Data['mag'].values  
    band = B_Data['band'].values
    
    #Seperating Background Data:
    Back_timer = []
    Back_timeg = []
    Back_Magr = []
    Back_Magg = []
    
    for i in range(0,len(time)-1,1):
        if band[i] == 'r':
            Back_Magr = np.append(Back_Magr, mag[i])
            Back_timer = np.append(Back_timer, time[i])
        elif band[i] == 'g':
            Back_Magg = np.append(Back_Magg, mag[i])
            Back_timeg = np.append(Back_timeg, time[i])
    
    #Background Magnitude Average Before Shockbreakout
    Average_Mag_r = np.round(np.mean(Back_Magr),2)
    print('Background is: ', Average_Mag_r)
    #Average_Mag_r = 18.33 #20.5
    ### Loading in Data 
    data = pd.read_csv(str(SN_dat),sep="\s+", names=['time', 'date', 'hour time','band','mag', 'unc sign', 'dmag', 'status', 'images1', 'image2', 'images3']) #<=== When importing the data, if the txt, csv, or dat file are missing a file
                                                                                            #just remove the column title from the name section, if there is more, simply just add the title to the names
    
    # Shift time data for the first point to start at zero
    time = data['time'].values
    mag = data['mag'].values  
    band = data['band'].values
    dmag = data['dmag'].values
    
    #Seperating Data bands:
    Time_r = []
    Time_g = []
    Mag_r = []
    Mag_g = []
    Dmag = []
    
    for i in range(0,len(time)-1,1):
        if band[i] == 'r':
            Mag_r = np.append(Mag_r, mag[i])
            Time_r = np.append(Time_r, time[i])
            Dmag = np.append(Dmag,dmag[i])
        elif band[i] == 'g':
            Mag_g = np.append(Mag_g, mag[i])
            Time_g = np.append(Time_g, time[i])
            
    Time_r = Time_r - Time_r[0]
    j = 0
    for i in range(len(Time_r)):
        if Time_r[i] <= 25:
            j = j+1
    print(j)
    Peak_lum = j
         
    ### Interpolation of SN
    # Setting the interval for the Quartic Interpolation
    x = Time_r[0:Peak_lum]
    print(x)
    y = Mag_r[0:Peak_lum]
    print(y)
    Dmag = Dmag[0:Peak_lum]
    
    M2 = y[np.argmin(y)]
    M1 = Average_Mag_r
    Flux_ratio = 1/ (100**((M1-M2)/5))
    Peak_Lum_Time = x[np.argmin(y)]
    
    print('Flux Ratio: ', Flux_ratio)
    print('Peak Luminosity Time: ', Peak_Lum_Time)
    
    
     # Fitted Equation (Quartic)
    def line(x, a, b, c, d, e):
         
         return a*x**4 + b*x**3 + c*x**2 + d*x + e 
    
     ### Fitting The Data
     # Extracting the parameters from curve fitting the data of time and flux to the quartic fit
    params, _ = scipy.optimize.curve_fit(line, x, y)
    
     # Extracting parameters for the Function
    a, b, c, d, e = params
    
     # Range of data for fitted function inputted into fitted function
    X1 = np.arange(-20,x[-1],.0001)  
    Y1 = a*X1**4 + b*X1**3 + c*X1**2 + d*X1 + e
    
     # Plotting Single Fit of Data
    plt.plot(X1,Y1,'--', color = 'red', label = 'Quartic Interpolation' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    plt.plot(Time_r,Mag_r, '.', color = 'black', label = 'Data for ' + str(SN_dat)) #<=== The flux is normalized to allow for us to use the point of x as our estimated time of shock breakout (for y=0)
    
    
     # Finding where the quartic fit intersects the x-axis which will serve as our estimated point of shock breakout
    index = np.argmin(np.abs(Y1 - Average_Mag_r))
     #print(index)
    mu = X1[index]
     #print(mu)
    
    
    # Estimated Time of ShockBreakout
    print('Time of Shock Breakout = ', np.abs(mu))
    
     ### Plotting the Data to Show the OSW on the Graph
    plt.axvline(mu, ymax = .2, color = 'blue', label = 'tsbo value')
    plt.ylabel('Magnitude')
    plt.gca().invert_yaxis()
    plt.xlabel('Days')
    plt.title(str(SN_dat))
    plt.minorticks_on()
    # X and Y limits vary depending on the light curve
    plt.ylim(Average_Mag_r,y[-1]-1.25)
    plt.xlim(mu-1,40)
    plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
    plt.legend(fontsize = "x-small", loc = 4)
    #plt.savefig(str(SN_dat) + '.png', dpi = 300, bbox_inches = "tight")
    plt.show()
    plt.close()
    
    def derivative(x,a,b,c,d):
       return 4*a*x**3 + 3*b*x**2 + 2*c*x + d
    
    
    Tsbo = mu
    Trise = np.abs(derivative(mu,a,b,c,d)) #Tangential Slope
    Slope_Rise = np.abs((Y1[-1]-Y1[index])/(X1[-1]-X1[index]))
    Rise_Time = X1[-1]-X1[index]
    Discovery_Time = x[0] - mu
    
    return [Tsbo,Trise,Slope_Rise,Rise_Time, Average_Mag_r, Flux_ratio, Discovery_Time]
    
#Storing Multiple Data Samples    
import numpy as np
import matplotlib.pyplot as plt

Tsbo = [] 
Trise = [] 
Slope_Rise = []
Rise_Time = []
Average_Mag_r = []
Flux_ratio = []
Discovery_Time = []

### Using the Quartic Package on Supernovae Data:
SN = QuarticLightCurveFit('SN2022fuc.txt','SN2022fuc_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])
    
SN = QuarticLightCurveFit('SN2022hrs.txt','SN2022hrs_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2022mxv.txt','SN2022mxv_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2022ngb.txt','SN2022ngb_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2022jzc.txt','SN2022jzc_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2022ihz.txt','SN2022ihz_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2022pgf.txt','SN2022pgf_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2021rhu.txt','SN2021rhu_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2021acvl.txt','SN2021acvl_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2020hvf.txt','SN2020hvf_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2020hvp.txt','SN2020hvp_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

SN = QuarticLightCurveFit('SN2020jfo.txt','SN2020jfo_Background.txt') 
Tsbo = np.append(Tsbo,SN[0])
Trise = np.append(Trise,SN[1])
Slope_Rise = np.append(Slope_Rise,SN[2])
Rise_Time = np.append(Rise_Time,SN[3])
Average_Mag_r = np.append(Average_Mag_r,SN[4])
Flux_ratio = np.append(Flux_ratio,SN[5])
Discovery_Time = np.append(Discovery_Time, SN[6])

plt.figure()
y, x, _ =plt.hist(Tsbo)
plt.title('Time of Shockbreakout Histogram')
plt.xlabel('Estimated Time of Shockbreakout')
mu=np.mean(Tsbo)
sigma=np.std(Tsbo)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.legend()
plt.grid()
plt.savefig('Tsbo_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

plt.figure() #Since Magnitude is inverted
y, x, _ =plt.hist(Trise)
plt.title('Tangential Slope of Rise Time Histogram')
plt.xlabel('Tangential Slope of Rise Time')
mu=np.mean(Trise)
sigma=np.std(Trise)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.legend()
plt.grid()
plt.savefig('Trise_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

plt.figure() #Since Magnitude is inverted
y, x, _ =plt.hist(Slope_Rise)
plt.title('Slope of Rise Time Histogram')
plt.xlabel('Slope of Rise Time')
mu=np.mean(Slope_Rise)
sigma=np.std(Slope_Rise)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.vlines(0.344,0,(3/5) * y.max(), color = "Green", label = 'Kepler Rise Time Slope')
plt.vlines(0.344*1.65,0,(2/5) * y.max(), color = "purple", label = 'Kepler Rise Time Slope +/- 65%')
plt.vlines(0.1204,0,(2/5) * y.max(), color = "purple")
plt.legend()
plt.grid()
plt.savefig('Slope_Rise_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

plt.figure() #Since Magnitude is inverted
y, x, _ =plt.hist(Rise_Time)
plt.title('Rise Time Histogram')
plt.xlabel('Rise Time')
mu=np.mean(Rise_Time)
sigma=np.std(Rise_Time)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.legend()
plt.grid()
plt.savefig('Rise_Time_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

plt.figure() #Since Magnitude is inverted
y, x, _ =plt.hist(Average_Mag_r)
plt.title('Background r-band Hist')
plt.xlabel('Background r-band')
mu=np.mean(Average_Mag_r)
sigma=np.std(Average_Mag_r)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.legend()
plt.grid()
plt.savefig('Background_r-band_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()

plt.figure() #Since Magnitude is inverted
y, x, _ =plt.hist(Discovery_Time)
plt.title('Discovery Time Histogram')
plt.xlabel('Discovery Time')
mu=np.mean(Discovery_Time)
sigma=np.std(Discovery_Time)
plt.vlines(mu+sigma,0,(2/5) * y.max(), color = "Red", label = r"$\sigma=$" + str(np.round(sigma,4)) )
plt.vlines(mu-sigma,0,(2/5) * y.max(), color = "Red")
plt.vlines(mu,0,(3/5) * y.max(), color = "black", label=r"$\mu=$" + str(np.round(mu,4)) )
plt.legend()
plt.grid()
plt.savefig('Discovery_Time_Hist.png', dpi = 300, bbox_inches = "tight")
plt.show()
plt.close()