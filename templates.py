import numpy as np
import os as os
import math
from scipy.interpolate import interp1d

redo_path="/media/lorenzo/LORENZO/AjusteLineas/redo_cont/"


############################################
"""
MUY IMPORTANTE

EL FLUJO DEBE DARSE EN F_LAMBDA Y LA LONGITUD DE ONDA
EN AMSTRONGS PARA QUE SHERPA AJUSTE ADECUADAMENTE PUES
LOS ESPECTROS SE ENCUENTRAN EN ESTAS MISMAS UNIDADES.

"""
############################################






def ABS(pars,x) :
    # Modelo de extincion de los datos basado en el de la SMC de Gordon+03 
    # (ver archivo phabs_SMC_A_sigma.dat para mas info).
    # Es un modelo multiplicativo, por lo que si queremos oscurecer la variable A
    # tenemos que ajustar ABS*A
    #
    #  Los parametros a ajustar son dos:
    #
    #     nht= columna de NH, es decir, el propio nivel de oscurecimiento. El valor
    #     de salida hay que multiplicarlo por unos factores para obtener Av o NH
    #      
    #     Av = nht * 2.74 * 0.76
    #     NH = nht * 2.74 * 10**22
    #
    #
    (nht,zz) = pars
    nh=nht * 1e22
    ##Cargamos el archivo como variables globales para optimizar
    if 'AA_absfile' in globals():
        ##Si ya existen estos pasos se saltan
        pass
    else:
        ##Si no existen se carga el archivo, se declaran las variables globales y se cargan longitudes de onda y flujo
        absfile = "/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/templates/phabs_SMC_A_sigma.dat"
        global AA_absfile,sigma_absfile
        AA_absfile, sigma_absfile = np.loadtxt(absfile, usecols=(0, 1), unpack=True)
    
    sigma = np.interp(x, AA_absfile*(1.+zz), sigma_absfile) ##LA interpolación no consume en exceso así que ignoramos su optimización
    return np.exp(-sigma*nh) ##Devolvemos termino de absorción






def cont(pars,x) :
    # Modelo aditivo simple de una ley de potencias y=A*x**B, donde:
    #
    #     A = Factor de normalizacion
    #     B = Indice de la ley de potencias
    #     
    (A,B)=pars
    return A*(x)**(B)

def gauss(pars,x) :
    # Modelo aditivo de una gaussiana adaptado para modelar lineas de un espectro,
    # donde tenemos tres parametros de entrada:
    #
    #     FO = Flujo de la gaussiana
    #     DO = Anchura en velocidad km/s de la linea
    #     LO = Centroide de la gaussiana en Angstrongs
    #     
    (FO,DOf,LO)=pars
    DO=(LO*1e3*DOf)/(3e8*2.35486)
    return FO*np.exp(-0.5*(x-LO)**2/(DO**2))/(np.sqrt(2*3.14159)*DO)
  
def AGNb(pars,x) :
    # Para modelar la emisión del AGN en el continuo se usa
    # la plantilla de Richards 2006 que no tiene en cuenta
    # emisión de lineas o pseudocontinuo.
    
    (ampl,zz)=pars
    ##Cargamos el archivo como variables globales para optimizar
    if 'lambda_A_agnt' in globals():
        pass
    else:
        global lambda_A_agnt,y_agnt
        agnt ="/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/templates/richards_flambda.txt"
        lambda_A_agnt, y_agnt = np.loadtxt(agnt, usecols=(0, 1), unpack=True)
        
    flux = np.interp(x, lambda_A_agnt*(1.+zz), y_agnt)
    return ampl*flux*1e-17/351671.948563006



def Gal(pars,x) :
    (K,zz)=pars
    # Modelo simple de emision estelar de galaxia, con dos parametros:
    #
    #     K = Factor de normalizacion
    #     zz = Redshift del objeto
    #           
    ##Cargamos el archivo como variables globales para optimizar
    if 'xx_mo' in globals():
        pass
    else:
        global xx_mo,yy_mo
        #mo = "/home/ordovas/XSH_def/sherpafits_lines/agncont/ag_j02.dat"
        #mo = "/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/templates/BC03_25.sed"
        mo = "/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/templates/S0_template_norm.sed"
        xx_mo,yy_mo = np.loadtxt(mo, usecols=(0, 1), unpack=True)    
    gg=np.interp(x, xx_mo*(1+zz), yy_mo)
    
    ##Aplicamos extinción
    
    return K*gg

  

def printeaNombres():
    print(temp_name_fe+"\n")
    print(temp_name_BAL+"\n")
    print(temp_name_GAL+"\n")

def borrarGlobales():
    ##Indicamos que se uso
    print('We have been ussing:\n')
    printeaNombres()
    
    ##We clean global variables related to fe and balmer continuums
    if 'flux_fe_template' in globals():
        global flux_fe_template
        global long_fe_template
        global temp_name_fe
        temp_name_fe=""
        flux_fe_template=0
        long_fe_template=0
        del temp_name_fe
        del flux_fe_template
        del long_fe_template
        
    if 'flux_BAL_template' in globals():
        global flux_BAL_template
        global long_BAL_template
        global temp_name_BAL
        temp_name_BAL=""
        flux_BAL_template=0
        long_BAL_template=0
        del temp_name_BAL
        del flux_BAL_template
        del long_BAL_template
        
    if 'flux_GAL_template' in globals():
        global flux_GAL_template
        global long_GAL_template
        global temp_name_GAL
        temp_name_GAL=""
        flux_GAL_template=0
        long_GAL_template=0
        del temp_name_GAL
        del flux_GAL_template
        del long_GAL_template



