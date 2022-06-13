# -*- coding: utf-8 -*-

#!/usr/bin/python
"""
Created on Thu Oct 03 10:22:58 2019

@author: 0124053Z
"""


import time
import sys    
import matplotlib.pyplot as plt # For doing the plots
import os
import pandas as pd
    
#    from Tkinter import *
#    import Tkinter
#    import time
import time
import datetime
import csv
import os
#from tqdm import tqdm
from tqdm import tqdm_notebook as tqdm

import numpy as np # For data manipulation
import scipy # For data manipulation
import random
import matplotlib.pyplot as plt # For doing the plots
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing
from scipy.fftpack import fft
import lmfit
from lmfit.models import GaussianModel
import plotly
import plotly.graph_objs as go
import sys
#import std

#import numpy as np
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve

#import sys
#!{sys.executable} -m pip install rampy
#import rampy as rp #Charles' libraries and function
from openpyxl import Workbook
def lecture_fichier_oct2019(nom_dossier_entree,nom_fichier,nom_dossier_sortie):


#------------------------------------------------------------------------------
#fonction pour convetir les csv 
#------------------------------------------------------------------------------

    ## mise en forme du text sous forme de lignes
    nom_complet=nom_dossier_entree+'/'+nom_fichier
    #print(nom_complet)
    file=open(nom_complet,'r')
    lignes = file.readlines() # on a
    #extraction de l'abcisse et de l'ordonnée 
    x=lignes[21] # abscisse en str
    y=lignes[23] # ordonnée en str
    #conversion du str en txt
    x0=[float(i) for i in x.replace('\n',',').split(',') if i] 
    y0=[float(i) for i in y.replace('\n',',').split(',') if i] 
    x0=np.array(x0)
    y0=np.array(y0)
#    plt.plot(x0,y0)
    '''
    airPLS.py Copyright 2014 Renato Lombardo - renato.lombardo@unipa.it
    Baseline correction using adaptive iteratively reweighted penalized least squares
    
    This program is a translation in python of the R source code of airPLS version 2.0
    by Yizeng Liang and Zhang Zhimin - https://code.google.com/p/airpls
    Reference:
    Z.-M. Zhang, S. Chen, and Y.-Z. Liang, Baseline correction using adaptive iteratively reweighted penalized least squares. Analyst 135 (5), 1138-1146 (2010).
    
    Description from the original documentation:
    
    Baseline drift always blurs or even swamps signals and deteriorates analytical results, particularly in multivariate analysis.  It is necessary to correct baseline drift to perform further data analysis. Simple or modified polynomial fitting has been found to be effective in some extent. However, this method requires user intervention and prone to variability especially in low signal-to-noise ratio environments. The proposed adaptive iteratively reweighted Penalized Least Squares (airPLS) algorithm doesn't require any user intervention and prior information, such as detected peaks. It iteratively changes weights of sum squares errors (SSE) between the fitted baseline and original signals, and the weights of SSE are obtained adaptively using between previously fitted baseline and original signals. This baseline estimator is general, fast and flexible in fitting baseline.
    
    
    LICENCE
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
    '''
    
    def WhittakerSmooth(x,w,lambda_,differences=1):
        '''
        Penalized least squares algorithm for background fitting
        
        input
            x: input data (i.e. chromatogram of spectrum)
            w: binary masks (value of the mask is zero if a point belongs to peaks and one otherwise)
            lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background
            differences: integer indicating the order of the difference of penalties
        
        output
            the fitted background vector
        '''
        X=np.matrix(x)
        m=X.size
        i=np.arange(0,m)
        E=eye(m,format='csc')
        D=E[1:]-E[:-1] # numpy.diff() does not work with sparse matrix. This is a workaround.
        W=diags(w,0,shape=(m,m))
        A=csc_matrix(W+(lambda_*D.T*D))
        B=csc_matrix(W*X.T)
        background=spsolve(A,B)
        return np.array(background)
    
    def airPLS(x, lambda_=100, porder=1, itermax=15):
        '''
        Adaptive iteratively reweighted penalized least squares for baseline fitting
        
        input
            x: input data (i.e. chromatogram of spectrum)
            lambda_: parameter that can be adjusted by user. The larger lambda is,  the smoother the resulting background, z
            porder: adaptive iteratively reweighted penalized least squares for baseline fitting
        
        output
            the fitted background vector
        '''
        m=x.shape[0]
        #print(m)
        w=np.ones(m)
        for i in range(1,itermax+1):
            z=WhittakerSmooth(x,w,lambda_, porder)
            d=x-z
            dssn=np.abs(d[d<0].sum())
            if(dssn<0.001*(abs(x)).sum() or i==itermax):
                if(i==itermax): 
                    print('WARING max iteration reached!')
                    break
            w[d>=0]=0 # d>0 means that this point is part of a peak, so its weight is set to 0 in order to ignore it
            w[d<0]=np.exp(i*np.abs(d[d<0])/dssn)
            w[0]=np.exp(i*(d[d<0]).max()/dssn) 
            w[-1]=w[0]
        return z
    y_base=airPLS(y0,lambda_=75, porder=3, itermax=500)
    #y_base=airPLS(y0,lambda_=75, porder=3, itermax=500)
        #print(len(x0))
        #print(len(y0))
        #y_base=airPLS(y0,lambda_=200, porder=3, itermax=500)
    
    y_corr=y0-y_base
    #plt.plot(x0,y_corr)
    nom_sortie=(nom_dossier_sortie+"/treated_"+nom_fichier+".txt")
    data_treated = np.array([x0, y_corr])
    data_treated=np.transpose(data_treated)
        #print(data_treated)
    
    with open(nom_sortie, 'w') as f:
        csv.writer(f, delimiter=' ').writerows(data_treated)
    return nom_sortie #Variable contenant le nom du dossier contenant le spectre traité





# Ce programme permet de calculer de DE d'un produit à partir de la lecture du fichier csv traité par l'algorithme précédent; ...
#On suppose que le fichier a déjà été traité dans l'algorithme de traitement et sauvegardé dans un dossier "treated_data", sinon le script ne marche pas

import numpy as np 
import sys
import time
    
import matplotlib.pyplot as plt # For doing the plots
import sys
    
    

import time
import time
import datetime
import csv
import os
#from tqdm import tqdm
from tqdm import tqdm_notebook as tqdm

import numpy as np # For data manipulation
import scipy # For data manipulation
import random
import matplotlib.pyplot as plt # For doing the plots
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.spatial import ConvexHull
import scipy.sparse as sparse
from numpy.linalg import norm
from sklearn import preprocessing
from scipy.fftpack import fft
import lmfit
from lmfit.models import GaussianModel
import plotly
import plotly.graph_objs as go
import sys
#import std

#import numpy as np
from scipy.sparse import csc_matrix, eye, diags
from scipy.sparse.linalg import spsolve

#import sys
#!{sys.executable} -m pip install rampy
#import rampy as rp #Charles' libraries and function
from openpyxl import Workbook

def calcul_DE_version_oct2019(nom_fichier,nom_dossier_treated):
    nom_sortie="treated_"+nom_fichier+".txt"
    #print(nom_sortie) #Chemin d'accès du fichier
    nom_treated=nom_dossier_treated+"/"+nom_sortie
    file1=np.loadtxt(nom_treated)
    x=file1 [:,0]
    y=file1 [:,1]
    #plt.plot(x,y,label='spectre_initial')
    def gaussian(x,amp,freq,HWHM): # for spectral fit
            #return amp*np.exp(-np.log(2)*((x-freq)/HWHM)**2)
            return (amp/np.sqrt(2*np.pi*HWHM))*np.exp(-np.log(2)*((x-freq)/HWHM)**2)
            #return amp*np.exp(-((x-freq)/HWHM)**2)
    def residual(pars, x, data=None, eps=None): #Function definition
        # unpack parameters, extract .value attribute for each parameter
        a1 = pars['a1'].value
        a2 = pars['a2'].value
        #a3 = pars['a3'].value
        #a4 = pars['a4'].value
        #a5 = pars['a5'].value
        
        f1 = pars['f1'].value
        f2 = pars['f2'].value
        #f3 = pars['f3'].value
        #f4 = pars['f4'].value
        #f5 = pars['f5'].value 
        
        l1 = pars['l1'].value
        l2 = pars['l1'].value
        #l3 = pars['l3'].value
        #l4 = pars['l4'].value
        #l5 = pars['l5'].value
        
        # Using the Gaussian model function from rampy
        peak1 = gaussian(x,a1,f1,l1)
        peak2 = gaussian(x,a2,f2,l1)
        #peak3 = rp.gaussian(x,a3,f3,l3)
        #peak4 = rp.gaussian(x,a4,f4,l4)
        #peak5 = rp.gaussian(x,a5,f5,l5)
        
        #model = peak1 + peak2 + peak3 + peak4 + peak5 # The global model is the sum of the Gaussian peaks
        #model = peak1 + peak2 + peak3  # The global model is the sum of the Gaussian peaks
        model = peak1 + peak2   # The global model is the sum of the Gaussian peaks
        #model = peak1   # The global model is the sum of the Gaussian peaks
        
        if data is None: # if we don't have data, the function only returns the direct calculation
            #return model, peak1, peak2, peak3, peak4, peak5
            #return model, peak1, peak2, peak3
            return model, peak1, peak2
            #return model, peak1
        if eps is None: # without errors, no ponderation
            return (model - data)
        return (model - data)/eps # with errors, the difference is ponderated
    
    def DE_pics_determination(x,y):
        # signal selection
        #lb = 2850 # The lower boundary of interest
        #hb = 2950 # The upper boundary of interest
        #lb = 2880 # The lower boundary of interest #oct
        #hb = 2955 # The upper boundary of interest #oct
        lb = 2900 # The lower boundary of interest #nov
        #	lb = 2895 # The lower boundary of interest #nov
        hb = 2950 # The upper boundary of interest #nov

        x_fit = x[np.where((x > lb)&(x < hb))]
        #print(len(x_fit))
        y_fit = y[np.where((x > lb)&(x < hb))]
        #ese0 = np.sqrt(abs(y_fit[:,0]))/abs(y_fit[:,0]) # the relative errors after baseline subtraction
        #y_fit[:,0] = y_fit[:,0]/np.amax(y_fit[:,0])*10 # normalise spectra to maximum intensity, easier to handle 
    
        #y_fit[:,0] = y_fit[:,0]/np.amax(y_fit[:,0]) # normalise spectra to maximum intensity, easier to handle 
        #y_fit = y_fit/np.amax(y_fit) # normalise spectra to maximum intensity, easier to handle 
        y_fit = y_fit/(np.amax(y_fit)*np.sqrt(2*np.pi)) # normalise spectra to maximum intensity,
        params = lmfit.Parameters()
        params.add_many(('a1',   0.8,   True,  0.5,      None,  None),
                        #('f1',   2900,   True, 2850,    2915,  None),
                        #('f1',   2900,   True, 2880,    2915,  None),
                        ('f1',   2900,   True, 2892,    2925,  None),
                        #('f1',   900,   True, 850,    915,  None),
                        #('l1',   40,   True,  1,      40,  None),
                        #('l1',   19,   True,  1,      25,  None), #valeur initiale 26/06/2016
                        #('l1',   22,   True,  22,      22.5,  None), #valeur finale 26/06/2016
                        #('l1',   22,   True,  22,      22.0001,  None), #valeur finale 26/06/2016
                        ('l1',   22,   True,  22,      22.0000000001,  None), #valeur finale 26/06/2016
                        #('l1',   40,   True,  1,      60,  None),
                        ('a2',   0.3,   True,  0,      None,  None),
                        ('f2',   2948,  True, 2948,  2965,  None))
                        #('f2',   2948,  True, 2945,  2965,  None))
                        #('f2',   948,  True, 880,  990,  None),
                        #('l2',   30,   True,  1,   120,  None),
                        #('l2',   30,   True,  1,   40,  None))
                        #('l2',   30,   True,  1, None,  None))
                        #('l2',   30,   True,  1,   15,  None))
    
        # we constrain the positions
        params['f1'].vary = False
        params['f2'].vary = False
        algo = 'nelder'  
    
        result = lmfit.minimize(residual, params, method = algo, args=(x_fit, y_fit)) # fit data with  nelder model from scipy
        #And now we release the frequencies:
        # we release the positions but contrain the FWMH and amplitude of all peaks 
        params['f1'].vary = True
        params['f2'].vary = True
        #params['f3'].vary = True
        #params['f4'].vary = True
        #params['f5'].vary = True
        
        
        
        #we fit twice
        result2 = lmfit.minimize(residual, params,method = algo, args=(x_fit, y_fit)) # fit data with leastsq model from scipy
        #result2 = lmfit.minimize(residual, params,method = algo, args=(x_fit, y_fit[:,0])) # fit data with leastsq model from scipy
        model = lmfit.fit_report(result2.params)
        #yout, peak1,peak2,peak3,peak4,peak5 = residual(result2.params,x_fit) # the different peaks
        #yout, peak1,peak2,peak3 = residual(result2.params,x_fit) # the different peaks
        #yout, peak1 = residual(result2.params,x_fit) # the different peaks
        yout, peak1,peak2= residual(result2.params,x_fit) # the different peaks
        #rchi2 = (1/(float(len(y_fit))-15-1))*np.sum((y_fit - yout)**2/sigma**2) # calculation of the reduced chi-square   
    
        #And let's have a look at the fitted spectrum:
    
        #print(result2.params)
    
        position_pic1=result2.params['f1'].value
        #print(position_pic1)
        std_pic1=result2.params['l1'].value
        #print(std_pic1)​x_fit=x0
        #y_fit=y0
        Ampl_pic1=result2.params['a1'].value
        #print(Ampl_pic1)
    
        position_pic2=result2.params['f2'].value
        #print(position_pic2)
        #std_pic2=result2.params['l2'].value
        std_pic2=result2.params['l1'].value
        #print(std_pic2)
        Ampl_pic2=result2.params['a2'].value
        #print(Ampl_pic2)
        
        plt.figure()
        
        #(x,amp,freq,HWHM)
        x_fit1 = x[np.where((x > lb -100)&(x < hb+ 100))]
        #print(len(x_fit))
        y_fit1 = y[np.where((x > lb -100 )&(x < hb+100))]
        y_fit1 = y_fit1/(np.amax(y_fit1)*np.sqrt(2*np.pi)) # normalise spectra to maximum intensity
        plt.plot(x_fit1,y_fit1,label='Experimental'+str(a))
        y_minim1=gaussian(x_fit,Ampl_pic1,position_pic1,std_pic1)
        y_minim2=gaussian(x_fit,Ampl_pic2,position_pic2,std_pic2)
        plt.plot(x_fit,y_minim1,label='Gaussienne 1')
        plt.plot(x_fit,y_minim2,label='Gaussienne 2')
        plt.plot(x_fit,y_minim1+y_minim2,label='Fitting total')
        plt.legend()
        plt.show()
        #Méthode finale
        #15088,92683	-12,66614647	-6,206947943	-155,640838	7,483911183	-7,023615644	135,5197005
    
        #DE_calc=20078.3836+(-5.918841577*position_pic1)+(2.612297294*std_pic1)+(143.5549507*Ampl_pic1)+(-1.035928026*position_pic2)+(0.*std_pic2)+(0.*Ampl_pic2)
        #print("la valeur du DE (méthode 1) de votre échantillon est DE=")
        #print("%6.1f"% (DE_calc))
        #DE_calc=15088.92683+(-12.66614647*position_pic1)+(-6.206947943*std_pic1)+(-155.640838*Ampl_pic1)+(7.483911183*position_pic2)+(-7.023615644*std_pic2)+(135.5197005*Ampl_pic2)
    
        #print("la valeur du DE (méthode 2)  de votre échantillon est DE=")
        #print("%6.1f"% (DE_calc))
    
        #DE_calc=20741.68611+(-9.676020777*position_pic1)+(-2.717640008*std_pic1)+(-173.8527945*Ampl_pic1)+(2.627880916*position_pic2)+(-6.373138495*std_pic2)+(26.41613026*Ampl_pic2)
        #print("la valeur du DE (méthode 2)  de votre échantillon est DE=")
        #print("%6.1f"% (DE_calc))
    
        #						
        #DE_calc=17598.99527+(-6.401846916*position_pic1)+(-6.261864846*std_pic1)+(313.9386413*Ampl_pic1)+(0.296675139*position_pic2)+(0*std_pic2)+(13.07141356*Ampl_pic2)
        #print("la valeur du DE (méthode 3)  de votre échantillon est DE=")
        #print("%6.1f"% (DE_calc))
        #Constante	14751.72877
        #position_pic1 (cm-1)	-11.55777106
        #position_pic2 (cm-1)	6.392540129
        #Ampl_pic1	-8.691608033
        #Ampl_pic2	12.48943612
        #
    #14597,12253	-13,02721498	-28,28297024	7,899992344	31,08603989
    
        #DE_calc=14597.12253+(-13.02721498*position_pic1)+(0*std_pic1)+(0*Ampl_pic1)+(7.899992344*position_pic2)+(0*std_pic2)+(0*Ampl_pic2) #Modèle 29/08/2019
        #DE_calc=14597.12253+(-13.02721498*position_pic1)+(0*std_pic1)+(-28.28297024*Ampl_pic1)+(7.899992344*position_pic2)+(0*std_pic2)+(31.08603989*Ampl_pic2) #Modèle 29/08/2019
        #DE_calc=17237.3201+(-11.43576401*position_pic1)+(0*std_pic1)+(-8.738422635*Ampl_pic1)+(5.425883902*position_pic2)+(0*std_pic2)+(14.53709603*Ampl_pic2) #Modèle 18/09/2019
        #DE_calc=15364.5194+(-12.0427443*position_pic1)+(0*std_pic1)+(0*Ampl_pic1)+(6.66181241*position_pic2)+(0*std_pic2)+(0*Ampl_pic2) #Modèle 18/09/2019
        DE_calc=15364.5194+(-12.0427443*position_pic1)+(0*std_pic1)+(0*Ampl_pic1)+(6.66181241*position_pic2)+(0*std_pic2)+(0*Ampl_pic2) #Modèle 18/09/2019
	#DE_calc=15364.5194+(-12.0427443*position_pic1)+(0*std_pic1)+(0*Ampl_pic1)+(6.66181241*position_pic2)+(0*std_pic2)+(0*Ampl_pic2) #Modèle 18/09/2019 #DEMO OCT
        #DE_calc=14453.80673+(-11.65776746*position_pic1)+(0*std_pic1)+(0*Ampl_pic1)+(6.592640864*position_pic2)+(0*std_pic2)+(0*Ampl_pic2)
        #DE_calc=14751.72877+(-11.55777106*position_pic1)+(0*std_pic1)+(-8.691608033*Ampl_pic1)+(6.392540129*position_pic2)+(0*std_pic2)+(12.48943612*Ampl_pic2)
        #print("les paramètres des gaussiennes permettant de déterminer le DE sont :")
        #print"position_pic1","std_pic1","Ampl_pic1","position_pic2","std_pic2","Ampl_pic2"
        #print position_pic1,std_pic1,Ampl_pic1,position_pic2,std_pic2,Ampl_pic2
        #print("la valeur du DE de votre échantillon est DE=")
        #print("%6.1f"% (DE_calc))
        #print(position_pic1,position_pic2)
        return position_pic1,std_pic1,Ampl_pic1,position_pic2,std_pic2,Ampl_pic2,DE_calc
    
    res=DE_pics_determination(x,y)
#    position_pic1=res[0]
#    std_pic1=res[1]
#    Ampl_pic1=res[2]
#    position_pic2=res[3]
#    std_pic2=res[4]
#    Ampl_pic2=res[5]
    DE_cal=round(res[6],1)
    
    return DE_cal




#fichier='ModbusOutputFile-Ch5_2021-07-13-23.19.46.920_stored_20210713_23_19_54.csv'

#fichier='gluc12prel27septsortiefiltreds28deg60_2019-09-27-15.26.47.797.csv'
a=os.listdir('raw_data/')
# fichier='amidon_mais_E9412_40MS_T35_7s_avecycorr_2022-04-20-04.56.47.183.csv'
# d_entre='raw_data'
# d_sortie='treated_data'
# fichier_treadted=fichier
# print(lecture_fichier_oct2019(d_entre,fichier,d_sortie))
# print(calcul_DE_version_oct2019(fichier_treadted,d_sortie))
# plt.show()
# plt.legend()

Data_tot=np.array([[0 for i in range (0,2)] for j in range(0,np.size(a))] ,dtype=float)
b=pd.DataFrame()
for i in range (0,np.size(a)) : 
    d_entre='raw_data'
    d_sortie='treated_data'
    if 'avecycorr'in a[i] :
        print("Le spectre utilisé est ="+str(lecture_fichier_oct2019(d_entre,a[i],d_sortie)))
        DE=0
        DE=calcul_DE_version_oct2019(a[i],d_sortie)
        print("Le DE calculé est = "+str(DE))
        #df2.at[i,'nom']=a[i][2]
        b.at[i,'nom']=a[i]
        b.at[i,'DE']=DE        
        #Data_tot[i:]=DE,a[i]
        
        
    

