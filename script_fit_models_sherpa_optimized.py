##Cargamos sherpa
from sherpa.astro.ui import *
from sherpa.astro.io import read_ascii
#from sherpa.astro import io
#from pychips import *
#from pychips.hlui import *
##Templates
from templates import *
##Paquetes mixtos
import os as os
import numpy as np
import sys

##Para test de tiempo
from timeit import default_timer as timer

##Configuración Matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


def my_stat_func(data, model, staterror, syserror=None, weight=None): 
    # A simple function to replicate chi2
    fvec = ((data - model) / staterror)**2
    stat = fvec.sum()/100.
    return (stat, fvec)

def my_staterr_func(data):
    # A simple staterror function
    return get_staterror()


def fitModel():  
  
  
  ##Ruta archivos en general
  redo_path="/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/"
  ##Ruta espectros
  data="/home/ordovas/privatecode_astro_and_ml/sherpa_cont/scripts_isma/data"
  
  ##Abro un archivo de output de resumen
  abs_out=open(redo_path+"res/abstract.dat","w")
  ##Definimos encabezado
  abs_out.write("#IAU_NAME redshift nh A_AGN A_Fer A_BAL A_Gal deltaup_nh deltaup_A_AGN deltaup_A_Fe deltaup_A_BAL deltaup_A_Gal deltadown_nh deltadown_A_AGN deltadown_A_Fe deltadown_A_BAL deltadown_A_Gal\n")
    
  #La lista de ficheros de ajuste
  with open (redo_path+"files.txt", 'r+', encoding="utf-8") as myfile:
    datos=myfile.readlines()
  
 

  
  #load info of the sources: host galaxy trigger, z and av from sed
  zfile=redo_path+"model.txt"
  iii=-1
  no_host_galaxy,zv,av = np.loadtxt(zfile, usecols=(0,1,2), unpack=True)
  #zfile=redo_path+"fits/flux_est_list.txt"
  #flx,lx = np.loadtxt(zfile, usecols=(0,1), unpack=True)

  ##Nos limitamos a los ajustes de los espectros de ISIS pendiente
  for iii in (0,1):
      
    ##Cargamos el nombre del fichero y el redshift  
    a=datos[iii].split('.dat', 1 )
    st=a[0]
    zz=zv[iii]
    
    ############################################################################################
    ##Se eligen metodos de ajuste y chi2
    load_user_stat("mystat", my_stat_func)
    set_stat(mystat)
    #set_stat("leastsq")
    show_stat()
    set_method("levmar")
    #set_method("moncar")
    show_method() 
    ############################################################################################

    ##Cargamos el espectro
    load_data(1,redo_path+"data/"+st+".dat",3)
    #read_ascii(redo_path+"data/"+st+".dat",3)
    set_conf_opt("soft_limits",True)
    
    ##Indicamos que estamos ajustando
    print("FIT: ",st,zz, iii, no_host_galaxy[iii])
    
    #Cargamos nuestros modelos  
    load_user_model(AGNb, "AGNb")
    load_user_model(ABS, "nh")
    load_user_model(Gal, "Gal")

    
    
    ##Valores por defecto de la normalización de la galxia
    if int(no_host_galaxy[iii]) == 1:
       ##Si no hay galaxia
       gi=0.0
    else:
       ## si hay galaxia
       gi=0.0005

    ##Añadimos componentes
    add_user_pars("AGNb", ["amp_AGN","z1"], parvals=[20,zz], parmins=[0,0.],parmaxs=[1000.,5.], parfrozen=[0,1]) #flx[iii]*1
    add_user_pars("nh", ["av","z1"], parvals=[0.2,zz], parmins=[0,0.],parmaxs=[10,5.], parfrozen=[0,1])
    add_user_pars("Gal", ["amp_F","z1"], parvals=[gi,zz], parmins=[0,zz-0.01],parmaxs=[100,zz+0.01], parfrozen=[0,1])

    

    #ignore abs features
    ignore(6840,6950)
    ignore(7580,7700)
    ignore(8100,8350)
    ignore(8900,9010)
    ignore(9050,9175)
    ignore(9250,9800)
    
    #ignore Broad lines
    ig_b=150.
    ignore((6564.6-ig_b)*(1.+zz),(6564.6+ig_b)*(1.+zz)) #Ha
    ignore((4862-ig_b)*(1.+zz),(4862+ig_b)*(1.+zz))     #Hb
    ignore((2798-ig_b)*(1.+zz),(3100.)*(1.+zz))     #MgII    
    ignore((1909-ig_b)*(1.+zz),(1909+ig_b)*(1.+zz))     #CIII]
    
    #ignore narrow lines   
    ig_n=20.
    ignore((6300-ig_n)*(1.+zz),(6300+ig_n)*(1.+zz)) #OI
    ignore((5876-ig_n)*(1.+zz),(5876+ig_n)*(1.+zz)) #HeI
    ignore((4960-ig_n)*(1.+zz),(5008+ig_n)*(1.+zz)) #OIII
    ignore((3727-ig_n)*(1.+zz),(3727+ig_n)*(1.+zz)) #OII
    ignore((3426-ig_n)*(1.+zz),(3426+ig_n)*(1.+zz)) #NeV
    ignore((6716-ig_n)*(1.+zz),(6732+ig_n)*(1.+zz)) #[SII]
    
          
    ##Máximo Chi2
    set_conf_opt("max_rstat",300)
    
    ##Primero solo galaxia y AGN    
    set_source(1,"nh*(AGNb)+Gal")
    show_model()
    fit()
    result = get_fit_results()

    print('Resultandos obtenidos')


    
    print('Se pintan resultados')
    
    #####################################
    ##     Ploteamos los resultados    ##
    #####################################
      
    ##Volvemos a la carpeta donde se encuentran los espectros
    os.chdir(data+"/")
    
    ##Cargamos el espectro ajustado
    try:
        longi,flux,error=np.loadtxt(data+"/"+st+".dat",unpack=True)
    except Exception:
        longi,flux=np.loadtxt(data+"/"+st+".dat",unpack=True)
        error=np.zeros(len(longi))
    
    ##Rango de longitudes de onda del modelo ajustado
    x=get_resid_plot().x
      
    ##Obtenemos las contribuciones de las diferentes componentes
    ##Absorción
    absorption=ABS((result.parvals[0],zz),x)
    ##Componente continuo AGN
    y_AGN=AGNb((result.parvals[1],zz),x)
    ##Continuo Galaxia
    ##Si no hay host galaxy
    if int(no_host_galaxy[iii]) == 1:
        y_gal=Gal((0,zz),x)
    ##Si hay host galaxy
    else:
        y_gal=Gal((result.parvals[2],zz),x)
    ##Total absorbido
    y_abs=absorption*(y_AGN)
      
    ##Obtenemos los residuos
    residuos=get_resid_plot().y
    """  
    ##Definimos el plot
    fig, axs = plt.subplots(2,sharex=True, gridspec_kw={'hspace': 0,'height_ratios': [4, 1]}, figsize=(8, 6))
    
    ###Regiones no ajustadas o ajustadas excepcionalmente
    ##Sombreamos regiones de absorción
    shadow_zones=[(6840,6950),(7580,7700),(8100,8350),(8900,9010),(9050,9175),(9250,9800)]
    for i,pair in enumerate(shadow_zones):
        axs[0].axvspan(pair[0], pair[1], alpha=0.2, color='blue',label=  "_"*i + "Removed Areas")
    
    ##Sombreamos regiones de lineas estrechas
    shadow_zones=[((6564.6-ig_b)*(1.+zz),(6564.6+ig_b)*(1.+zz)),((4862-ig_b)*(1.+zz),(4862+ig_b)*(1.+zz)),((2798-ig_b)*(1.+zz),(3100.)*(1.+zz)),((1909-ig_b)*(1.+zz),(1909+ig_b)*(1.+zz))]
    for i,pair in enumerate(shadow_zones):
        axs[0].axvspan(pair[0], pair[1], alpha=0.2, color='blue')

    ##Sombreamos regiones de lineas anchas
    shadow_zones=[((6300-ig_n)*(1.+zz),(6300+ig_n)*(1.+zz)),((5876-ig_n)*(1.+zz),(5876+ig_n)*(1.+zz)),((4960-ig_n)*(1.+zz),(5008+ig_n)*(1.+zz))]
    shadow_zones.append(((3727-ig_n)*(1.+zz),(3727+ig_n)*(1.+zz)))
    shadow_zones.append(((3426-ig_n)*(1.+zz),(3426+ig_n)*(1.+zz)))
    shadow_zones.append(((6716-ig_n)*(1.+zz),(6732+ig_n)*(1.+zz)))
    for i,pair in enumerate(shadow_zones):
        axs[0].axvspan(pair[0], pair[1], alpha=0.2, color='blue')
    

      
    #####Espectro + ajuste
    ##Espectro
    axs[0].errorbar(longi,flux,yerr=error,ls='',marker='o', mfc='none', mec='black',ecolor='black',elinewidth=0.11,ms=2,zorder=0,mew=0.1)
    ##Comonentes
    axs[0].plot(x,y_gal,ls='--',color='maroon',label='Galaxy',zorder=30)
    axs[0].plot(x,y_AGN*absorption,ls='-.',color='blue',label='AGN Component',zorder=120)
    ##Total
    axs[0].plot(x,y_abs+y_gal,ls='-',color='red',label='Total',zorder=220)
    ##Legenda
    lgd = axs[0].legend(bbox_to_anchor=(1.04,1), loc="upper left")
    ##Tick y por denstro
    axs[0].tick_params(direction="in",top=True, right=True,labeltop=False,labelright=False, labelsize=7)
    axs[0].minorticks_on()
    axs[0].tick_params(which='minor', direction='in',top=True, right=True,labeltop=False,labelright=False, labelsize=7)
    ##Formato labels y ticks
    axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    axs[0].set_ylabel(r"$f_\lambda$ / ($\frac{\mathrm{erg}}{\mathrm{s\times cm^2\times \AA}}$)", fontsize=10)    
      
    ###Residuos
    axs[1].scatter(x,residuos,marker='o',facecolor='none',edgecolor='black',s=2,zorder=10,linewidths=0.3)
    axs[1].axhline(0,ls='--',color='black',zorder=1)
    ##Ejes
    axs[1].set_ylim([0.98*min(residuos),1.02*max(residuos)])
    ##Formato labels
    axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    axs[1].set_xlabel(r"$\lambda$ / $\AA$", fontsize=10)
    axs[1].set_ylabel(r"Residuos", fontsize=10)
    ##Tick y por denstro
    axs[1].tick_params(direction="in",top=True, right=True,labeltop=False,labelright=False, labelsize=7)
    axs[1].minorticks_on()
    axs[1].tick_params(which='minor', direction='in',top=True, right=True,labeltop=False,labelright=False, labelsize=7)
      
    ##Guardamos y cerramos
    plt.xlim(min(long),max(long))
    plt.savefig(redo_path+"res_redo/plot_"+str(int(iii))+"_"+st[0:-4]+".pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.savefig(redo_path+"res_redo/plot_"+str(int(iii))+"_"+st[0:-4]+".png",dpi=150, bbox_extra_artists=(lgd,), bbox_inches='tight')
    ##Guardamos una versión con zoom
    axs[0].set_ylim([0.05*np.mean(flux),3*np.mean(flux)])
    plt.savefig(redo_path+"res_redo/plot_"+str(int(iii))+"_"+st[0:-4]+"_zoomed.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.savefig(redo_path+"res_redo/plot_"+str(int(iii))+"_"+st[0:-4]+"_zoomed.png",dpi=150, bbox_extra_artists=(lgd,), bbox_inches='tight')
    
    ##Pintamos y cerramos
    plt.draw()
    plt.close()
    """
    ##Nos desplzamos a la carpeta de output
    os.chdir(redo_path+"res/")
    
    ##Obtenemos intervalos de confianza
    print("Calculo intervalos de confianza")
    conf()
    cov1 = get_conf_results()
    sep = "\n \n"
    ##Guardamos el modelo ajustado en un .dat
    print("Guardo ajuste")
    save_model(1, redo_path+"res/fittedx"+str(iii)+"_"+st+".dat", ascii=True)
    """
    ##Obtenemos valores de AV
    print("Calculo valores AV")
    int_proj("nh.av",min=0,max=5,nloop=150)
    intav=get_int_proj("nh.av",min=0,max=1.5,nloop=120)
    
    
    ##Guardamos archivo de AV
    print("Guardo valores AV")
    fileout= open("avint_"+str(iii)+"_"+st+".dat","w+")
    print(sep,file=fileout)
    print(intav,file=fileout)
    print(sep,file=fileout)
    fileout.close()
    """
    ##Guardamos los residuos
    print("Guardo residuos")
    save_resid(1,"resid"+str(iii)+"_"+st+".dat", ascii=True)


    ##Guardamos el archivo de resultados
    print("Guardo archivo de resultados")
    fileout= open("results"+str(iii)+"_"+st+".dat","w+")
    print(sep,file=fileout)
    print(result,file=fileout)
    print(sep,file=fileout)
    print(cov1,file=fileout)
    fileout.close()
    
    print('Modelo salvado')
    """
    ##Rellenamos correspondiente fila en el abstract
    ###ATENCIÓN: Código muy largo para abarcar todos los posibles casos
    ##Nombre AGN y redshift
    abs_out.write(st+" "+str(zz)+" ")
    ##Si no hay host galaxy
    if int(no_host_galaxy[iii]) == 1:
        ##Valores mejor ajuste
        abs_out.write(str(result.parvals[0]*2.74*1e22)+" "+str(result.parvals[1])+" "+str(result.parvals[2])+" "+str(result.parvals[3])+" "+str(0.0)+" ")
        ##Valores error superior (por la derecha)
        for index,up_error_value in enumerate(cov1.parmaxes):
            if index==0:
                try:
                    abs_out.write(" "+str(up_error_value*2.74*1e22))
                except TypeError:
                    abs_out.write(" "+str(-1))
            elif index==4:
                abs_out.write(" "+str(0.0))
            else:
                if str(up_error_value)=="None":
                    abs_out.write(" "+str(-1))
                else:
                    abs_out.write(" "+str(up_error_value))
        ##Valores error inferior (por la izquierda)
        for index,down_error_value in enumerate(cov1.parmins):
            if index==0:
                try:
                    abs_out.write(" "+str(down_error_value*2.74*1e22))
                except TypeError:
                    abs_out.write(" "+str(-1))
            elif index==4:
                abs_out.write(" "+str(0.0))
            else:
                if str(up_error_value)=="None":
                    abs_out.write(" "+str(-1))
                else:
                    abs_out.write(" "+str(up_error_value))
    ##Si hay host galaxy
    else:
        ##Valores mejor ajuste
        abs_out.write(str(result.parvals[0]*2.74*1e22)+" "+str(result.parvals[1])+" "+str(result.parvals[2])+" "+str(result.parvals[3])+" "+str(result.parvals[4]))
        ##Valores error superior (por la derecha)
        for index,up_error_value in enumerate(cov1.parmaxes):
            if index==0:
                try:
                    abs_out.write(" "+str(up_error_value*2.74*1e22))
                except TypeError:
                    abs_out.write(" "+str(-1))
            else:
                if str(up_error_value)=="None":
                    abs_out.write(" "+str(-1))
                else:
                    abs_out.write(" "+str(up_error_value))
        ##Valores error inferior (por la izquierda)
        for index,down_error_value in enumerate(cov1.parmins):
            if index==0:
                try:
                    abs_out.write(" "+str(down_error_value*2.74*1e22))
                except TypeError:
                    abs_out.write(" "+str(-1))
            else:
                if str(up_error_value)=="None":
                    abs_out.write(" "+str(-1))
                else:
                    abs_out.write(" "+str(up_error_value))
                    
    ##Salto de linea
    abs_out.write("\n")
        
    ##Borramos las variables globales de template correspondiente a
    ##las plantillas del hierro y de balmer.
    borrarGlobales()
    
  ##Cerramos el archivo de output
  abs_out.close()
"""


"""
Método para realizar un analisis del rendimiento temporal
del método principal. Pide el tamaño de la muestra usada
para el estudio del rendimiento y el nombre de la CPU para
dar una idea de que se uso.
"""
def simpleTimeTest(size_sample=10,cpu_name="Intel I7 7700K 4.20 GHz"):
    ##Creamos archivo para guardar la info
    time_test=open("time_test.dat","w")
    ##Contamos tiempo de inicio
    start = timer()
    ##Realizamos el ajuste
    fitModel()
    ##Contamos el tiempo de final
    end = timer()
    ##Calculamos tiempo tardado
    totaltime=end-start
    ##Escribimos el reporte
    time_test.write("Usando un procesador"+cpu_name+"\n")
    time_test.write("se ha tardado para un total de "+str(size_sample)+"\n")
    time_test.write("espectros un total de "+str(totaltime)+"segundos. \n\n")
    
    time_test.write("Eso supone:\n")
    time_test.write(str(totaltime/size_sample)+" segundos por objeto. \n")
    time_test.write(str(totaltime/size_sample/60)+" minutos por objeto. \n")
    
    
    


















