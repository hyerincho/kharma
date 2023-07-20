import numpy as np
import pyharm
import h5py
import glob
import matplotlib
import matplotlib.pyplot as plt
import pdb
from astropy.constants import G,c
from astropy import units as u
import os
import matplotlib.image as mpimg

import pyharm.plots.plot_dumps as pd  
import pyharm.plots.plot_results as pr
from pyharm.plots.figures import floors, plot_hst, e_ratio, prims
from pyharm.plots.frame import frame
from pyharm.ana_results import AnaResults
from pyharm import io

import bondi_analytic as bondi
from bondi_analytic import define_globals
###from turbulence import fft_comp

def matplotlib_settings():
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif'] 

def calc_Phi_b(dump,EH=None): # iEH at 2 or 5? 5 according to Chatterjee & Narayan
    if EH is not None: iEH = np.argmin(abs(dump["r1d"]-EH)) # TODO: expand to Kerr
    else: iEH=None
    Phi_b = 0.5 * pyharm.shell_sum(dump, 'abs_B1', at_i=iEH) # from ana/analyses.py
    return Phi_b*np.sqrt(4.*np.pi)#/np.sqrt(Mdot)

def calc_Mdot(dump,EH=None):
    if EH is not None: iEH = np.argmin(abs(dump["r1d"]-EH)) # TODO: expand to Kerr
    else: iEH=None
    Mdot = -pyharm.shell_sum(dump, 'FM', at_i=iEH)
    #Mdot = -pyharm.shell_sum(dump, 'FM', at_i=iEH, j_slice=pyharm.ana.reductions.get_j_slice(dump)) # for disk
    return Mdot

def plot_mdot_r(oz,nz,sz,dirtag="",skip=1,avg=False):
    matplotlib_settings()                 
    fig,ax=plt.subplots(1,1,figsize=(8,6))
      
    GHOST=False #True
    ONEZONE=oz #False #
    KERR=False # True #
    BFIELD=False #True # 
    num_zones=nz #430
    start_zone=sz #0

    colors=plt.cm.gnuplot(np.linspace(0.9,0.3,num_zones))

    # configure
    dirname=glob.glob("../data/"+dirtag+"*0*/")[0]
    fname1=sorted(glob.glob(dirname+"*.rhdf"))[0]
    dump = pyharm.load_dump(fname1,ghost_zones=GHOST)
    mdot_sum=[0]*dump["nzone"]
    num_sum=[0]*dump["nzone"]

    for i,zone in enumerate(range(start_zone,start_zone+num_zones,skip)):
        dirname="../data/"+dirtag+"*{:05d}/".format(zone)

        files=sorted(glob.glob(dirname+"*.phdf"))
        fname1=files[0]
        fname2=files[-1]

        time = dump['t']
        if num_zones<20 or i==0 or i==num_zones-1:
            label="t= {:.2g}".format(time) #"zone {}".format(i)
        else:
            label="__nolegend__"
        if avg:
          dump = pyharm.load_dump(fname1,ghost_zones=GHOST)
          zone_num = int(np.log10(dump["r_in"])/np.log10(int(dump["base"])))
          for file_ in files[len(files)//2:]:
            dump=pyharm.load_dump(file_,ghost_zones=False)
            mdot_all = calc_Mdot(dump)
            mdot_sum[zone_num]+=mdot_all
            num_sum[zone_num]+=1
          r_arr = dump['r1d'] 
        else:
          dump1 = pyharm.load_dump(fname1,ghost_zones=GHOST)
          dump = pyharm.load_dump(fname2,ghost_zones=GHOST)
          r_arr = dump['r1d'] #dump['r'][:,0,0] # only r values      
          mdot_1 = calc_Mdot(dump1)#-pyharm.shell_sum(dump1, 'FM')
          mdot_all = calc_Mdot(dump) #-pyharm.shell_sum(dump, 'FM')
          ax.plot(r_arr,mdot_all,label=label,color=colors[i],lw=4,ls='--',alpha=0.5)
          #ax.plot(r_arr,mdot_1,label=label,color=colors[i],lw=2,ls='-')

        # TODO do a more general transformation
        #rho_ucon_fake = np.zeros(np.shape(dump["ucon"])) # making fake rho_ucon because it nees to have same dimension as ucon in order to do coordinate transformation
        _,C1,_,_=define_globals(dump["rs"])
        #rho_ucon_fake[1] = C1/dump["r"]**2
        #rho_ucon = np.einsum("i...,ij...->j...", rho_ucon_fake, dump['dXdx'])[1]# rho * ucon
        # C1/r**2 is rho*u in base coordinates. should change to native coordinates
        #mdot_analytic = pyharm.shell_sum(dump, rho_ucon)
        if i==0: #ONEZONE and zone==0
            r_bondi = np.logspace(0,np.log10(dump["nzone"]+2),100)
            #ax.plot(r_arr,mdot_1,label="t= 0",color=colors[i],lw=4)
            ax.plot(r_bondi,C1*4.*np.pi*np.ones(np.shape(r_bondi)),'g:')  ## TEST
            #ax.plot(r_arr,mdot_analytic,'k:') # testing for onezone case for now

        
    if avg:
      print(num_sum)
      for i in range(dump["nzone"]):
        dx = np.log(int(dump["base"])**2)/len(r_arr)
        r_in = np.power(8,i)
        r_out = np.power(8,i+2)
        x1=np.arange(np.log(r_in),np.log(r_out),dx)+dx/2.
        r_sum = np.exp(x1)
        if len(r_sum)>len(r_arr):
          r_sum=r_sum[:-1]
        ax.plot(r_sum,mdot_sum[i]/num_sum[i],label=label,color=colors[0],lw=4,ls='--',alpha=0.5)
    ax.set_xlabel(pyharm.pretty('r'))   
    ax.set_ylabel(pyharm.pretty('Mdot'))
    ax.set_title("{} runs {} - {}".format(dirtag.replace("bondi_multizone_",""),start_zone,zone))
    ax.axvline(2)
    ax.legend()

    ax.set_xscale('log')
    #if "JE3" in quantity: #ax.get_ylim()[0] <0:
    #    ax.set_yscale('symlog')
    #else:
    ax.set_yscale('log')

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    save_path="./plots/"+dirtag
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    savefig_name='/mdot_r_final.png'
    #if ONEZONE:
    #    savefig_name=savefig_name.replace("final", "onezone_final")
    if KERR:
        savefig_name=savefig_name.replace("final", "kerr_final")
    plt.savefig(save_path+savefig_name,bbox_inches='tight')
    plt.savefig("./plots/"+savefig_name,bbox_inches='tight') # make a copy here as well

def plot_Ldot_r(sz,dirtag=""):
    matplotlib_settings()                 
    fig,ax=plt.subplots(1,1,figsize=(8,6))
      
    GHOST=True
    start_zone=sz #0

    start=0 #59
    until=-1 # 61#
    uphi=False #True # u^phi instead of Ldot

    dirname="../data/"+dirtag+"/bondi_multizone_{:05d}/".format(start_zone)
    files=glob.glob(dirname+"/*.phdf")
    files = sorted(files)[start:until]
    
    num=len(files)
    colors=plt.cm.gnuplot(np.linspace(0.9,0.3,num))

    for i,fname in enumerate(files):
        if 1: #if i < until or until <0:
            dump = pyharm.load_dump(fname,ghost_zones=GHOST)
            r_arr = dump['r1d']
            if uphi:
                Ldot = pyharm.shell_avg(dump, 'u^phi') # actuall ang vel
            else:
                Ldot = pyharm.shell_sum(dump, 'FL_Fl') # TODO: change it to ask if there's a bfield
            time = dump['t']

            if num<20 or i==0 or i==num-1:
                label="t= {:.5g}".format(time) #"zone {}".format(i)
            else:
                label="__nolegend__"
            ax.plot(r_arr,Ldot,label=label,color=colors[i],lw=2,ls='-',alpha=0.5)
            ax.plot(r_arr,-(Ldot),color=colors[i],lw=2,ls=':',alpha=0.5)

    ax.set_xlabel(pyharm.pretty('r'))
    if uphi:
        ax.set_ylabel(pyharm.pretty('u^phi'))
    else:
        ax.set_ylabel(pyharm.pretty('Ldot'))
    ax.set_title(dirtag+" run = {}".format(start_zone))
    #ax.axvline(2)
    ax.legend()

    ax.set_xscale('log')
    #if "JE" in quantity: #ax.get_ylim()[0] <0:
    #    ax.set_yscale('symlog')
    #else:
    ax.set_yscale('log')

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    save_path="./plots/"+dirtag
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if uphi:
        savefig_name=save_path+'/u^phi_r.png'
    else:
        savefig_name=save_path+'/Ldot_r_final.png'
    plt.savefig(savefig_name,bbox_inches='tight')

def plot_one_run(run,quantity,dirtag=""):
    # plot spherically averaged profiles
    matplotlib_settings()                 
    fig,ax=plt.subplots(1,1,figsize=(8,6))
      
    GHOST=True
    start_zone=run #0

    start=1 #0 #59
    until=-1 # 61#
    avg=True
    to_sum=quantity
    if quantity == "Ldot":
       avg=False
       to_sum="FL_Fl"
    if quantity == "Mdot":
        avg=False
        to_sum="FM"


    dirname="../data/"+dirtag+"/bondi_multizone_{:05d}/".format(start_zone)
    files=glob.glob(dirname+"/*.rhdf")
    files = sorted(files)[start:until]
    
    num=len(files)
    colors=plt.cm.gnuplot(np.linspace(0.9,0.3,num))

    for i,fname in enumerate(files):
        if 1: #if i < until or until <0:
            dump = pyharm.load_dump(fname,ghost_zones=GHOST)
            r_arr = dump['r1d']
            if avg:
                avged = pyharm.shell_avg(dump, to_sum)
            else:
                avged = pyharm.shell_sum(dump, to_sum)
            time = dump['t']

            if num<20 or i==0 or i==num-1:
                label="t= {:.5g}".format(time) #"zone {}".format(i)
            else:
                label="__nolegend__"
            ax.plot(r_arr,avged,label=label,color=colors[i],lw=2,ls='-',alpha=0.5)
            ax.plot(r_arr,-(avged),color=colors[i],lw=2,ls=':',alpha=0.5)

    ax.set_xlabel(pyharm.pretty('r'))
    ax.set_ylabel(pyharm.pretty(quantity))
    ax.set_title(dirtag+" run = {}".format(start_zone))
    #ax.axvline(2)
    ax.legend()

    ax.set_xscale('log')
    #if "JE" in quantity: #ax.get_ylim()[0] <0:
    #    ax.set_yscale('symlog')
    #else:
    ax.set_yscale('log')

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    save_path="./plots/"+dirtag
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    savefig_name=save_path+'/'+quantity+'_at_run{}'.format(start_zone)+'.png'
    print("saving figure at "+savefig_name)
    plt.savefig(savefig_name,bbox_inches='tight')

def r_average(dump,quantity,imin=0,sum_instead=False,mass_weight=True):
  if isinstance(quantity,str):
    to_average=np.copy(dump[quantity])
  else:
    to_average=np.copy(quantity)
  
  if mass_weight:
    to_average *= dump['rho']

  if sum_instead:
      return pyharm.shell_sum(dump,to_average)
  else:
      if mass_weight:
        return pyharm.shell_sum(dump,to_average)/pyharm.shell_sum(dump,dump["rho"])
      else:
          if dump['n3']>1: #3d
            return pyharm.shell_avg(dump,to_average) #np.mean(to_average[imin:,:,:],axis=(1,2)) #
          else:
            return np.mean(to_average[imin:,:],axis=1)

def th_average(dump,quantity,i):
  # average in theta direction given one radius
  volumetric_weight=dump["gdet"]
  if isinstance(quantity,str):
    to_average=dump[quantity]
  else:
    to_average=quantity
  return np.sum(to_average[i,:,:]*volumetric_weight[i,:,:],axis=1)/np.sum(volumetric_weight[i,:,:],axis=1)

def phi_average(dump,quantity,i):
  # average in theta direction given one radius
  volumetric_weight=dump["gdet"]
  if isinstance(quantity,str):
    to_average=dump[quantity]
  else:
    to_average=quantity
  return np.sum(to_average[i,:,:]*volumetric_weight[i,:,:],axis=0)/np.sum(volumetric_weight[i,:,:],axis=0)

def plot_single_dump(ax,dump,quantity,color,norm,lw=2,ls='-',alpha=1.,show_label=True):
    params = dump.params
    convergence = 0 # if converged, return the mean frac difference
    convergence = -1 # don't care about convergence for now. (not supported)

    r = dump['r'][:,dump['n2']//2,0] #0]
    #if dump['n3']>1:
      #r = dump['r'][:,dump['n2']//2,dump['n3']//2]
    #else:
    #  r = dump['r'][:,dump['n2']//2]

    imin = 0
    #while r[imin] < params['r_eh']:
    #    imin += 1

    r = r[imin:]
    if quantity != 'T' and quantity!='g_00' and quantity!="Ldot" and quantity!="beta": # TODO: must look for variable list and then make exceptions
      to_plot=r_average(dump,quantity,imin)
    
    if quantity=='u1': # and i==0 and not 'no_outflow' in run:
      to_plot[abs(to_plot)<1e-5]=np.nan
    if quantity=='u^r':
        convergence = -1 # don't care about convergence
        to_plot = -to_plot # only show infalling velocity
        #to_plot = abs(to_plot) # show magnitude of the velocity
    if quantity=='u^phi':
        to_plot=abs(to_plot)
        #to_plot = to_plot*r # transform from mks to ks vector
    time = dump['t']
    if norm and quantity=='u^r':
      #u=np.mean(dump['u'][imin:,:],axis=1)
      #rho=np.mean(dump['rho'][imin:,:],axis=1)
      gam=dump['gam']

      #n=1./(gam-1.)
      #cs=np.sqrt(u/(rho*n))
      # from Gammie 03
      w=dump['rho']+dump['u']+dump['p']
      cs=r_average(dump,np.sqrt(gam*dump['p']/w),imin)
      to_plot = abs(to_plot/cs)
      #ax.plot(r, abs(to_plot/cs), label='t= {:.2g}'.format(time),color=color)
    elif quantity == 'T':
      gam=dump['gam']
      to_plot=r_average(dump,(dump['u']/dump['rho']*(gam-1.)),imin)
    elif quantity=="g_00":
        to_plot = abs(r_average(dump,-1./2*(1.+(dump["gcov"][0,0])),imin))
    elif quantity=="Ldot":
        to_plot=r_average(dump,"FL_Fl",imin,sum_instead=True)
    elif quantity=="beta":
        to_plot=1./r_average(dump,1./dump["beta"],imin,mass_weight=True)
    #elif quantity == 'vr':
      #mdot = 1. # TODO: how to incorporate with Bondi?
      #to_plot=r_average(dump,(-dump['u^r']/dump['u^t']/(1.-2.*mdot/r[:,0])),imin)
      #ax.plot(r, T, label='t= {:.2g}'.format(time),color=color)
 
    nvals = 10 #5 # only pick out this number of r values
    threshold = 1e-2
    args = np.arange(0,len(r),int(len(r)/nvals)) 
    r_short = r[args]#,0] 
    analytic_sol = bondi.get_quantity_for_rarr(r_short,quantity,rs=dump["rs"])
    
    if convergence != -1 and analytic_sol is not None:
        frac_diff = (analytic_sol-to_plot[args])/analytic_sol
        #print(frac_diff)

        if np.all(abs(frac_diff) < threshold):
            convergence=abs(np.mean(frac_diff))
    else:
        convergece=-1 # negative number zero number
    
    if show_label or convergence>0:
        label="t= {:.2g}".format(time)
    else:
        label="__nolegend__"
    ax.plot(r, (to_plot), label=label,color=color,alpha=alpha,lw=lw,ls=ls)

    return convergence

def get_base_num(dirtag):
    fname0="../data/"+dirtag+"bondi_multizone_00000/*.phdf" # TODO: find

def plot_grid_1d(nz,dirtag,GHOST=True):
    matplotlib_settings()                 
    fig,ax=plt.subplots(1,1,figsize=(8,6))

    for zone in range(nz):
        dirname="../data/"+dirtag+"bondi_multizone_{:05d}/".format(zone)
        fname=glob.glob(dirname+"*.00000.phdf")
        if len(fname)>1:
            print("ERROR: more than one initial file found in "+dirname)
        else:
            fname=fname[0]
        dump = pyharm.load_dump(fname)#,ghost_zones=GHOST)
        to_plot=dump['th1d']
        ax.plot(np.arange(len(to_plot)),to_plot,'--',alpha=0.5)

    ax.set_ylabel(r'$\theta$')
    plt.savefig("./plots/grid_1d.png",bbox_inches='tight')

def plot_grid_3d(dirtag,GHOST=True):
    matplotlib_settings()                 
    fig=plt.figure()
    ax = plt.axes(projection='3d')
    
    zone=0
    dirname="../data/"+dirtag+"bondi_multizone_{:05d}/".format(zone)
    fname=glob.glob(dirname+"*.00000.phdf")
    if len(fname)>1:
        print("ERROR: more than one initial file found in "+dirname)
    else:
        fname=fname[0]
    dump = pyharm.load_dump(fname,ghost_zones=GHOST)

    # grid points
    idx=-1
    r=dump['r1d'][-1]
    theta=dump['th'][idx,:,0]
    phi=dump['phi'][idx,0,:]
    TH,PHI = np.meshgrid(theta,phi)
    X=r*np.sin(TH)*np.cos(PHI)
    Y=r*np.sin(TH)*np.sin(PHI)
    Z=r*np.cos(TH)
    ax.plot_wireframe(X,Y,Z)

    plt.savefig("./plots/grid_3d.png",bbox_inches='tight')

def plots_together(dirtag):
    fig,ax = plt.subplots(1,4,figsize=(32,6))
    plt_path="./plots/"+dirtag
    plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
    for i,quantity in enumerate(["rho", "abs_u^r", "abs_u^th", "abs_u^phi"]):
      plt_name=plt_path+"/multizone_check_"+quantity+".png"
      image = mpimg.imread(plt_name)
      ax[i].imshow(image)
      ax[i].axis('off')
      ax[i].axis("tight")  # gets rid of white border
      ax[i].axis("image")  # square up the image instead of filling the "figure" space
    #plt.tight_layout()
    plt.savefig(plt_path+"/multizone_plots.png", bbox_inches='tight', pad_inches = 0)

def plot_diffs(fname1,fname2,var):
    dump1=pyharm.load_dump(fname1)
    dump2=pyharm.load_dump(fname2)
    tdump1 = int(io.get_dump_time(fname1))
    tdump2 = int(io.get_dump_time(fname2))
    fig,ax=plt.subplots(1,2,figsize=(16,8))
    sz=dump1["r_out"]
    log_r=True
    log=True #False #
    rel=True #False #
    if log_r:
      sz=np.log10(sz)
    if log:
      vmax=1e-3#-5 0
      vmin=1e-10 #-12 -7
      #vmax=None
      #vmin=None
    else:
      vmax=1e-3 #1
      vmin=-1e-3 #-1
    kwargs={**{'window':[-sz,sz,-sz,sz]}, **{'vmin':vmin}, **{'vmax':vmax}, **{'log_r':log_r}, **{'log':log}} # vmin -7 vmax 0
    pd.plot_diff_xz(ax[0], dump1, dump2, var,rel=rel,**kwargs)
    pd.plot_diff_xy(ax[1], dump1, dump2, var,rel=rel,**kwargs)
    title="diff: \n  "+ fname1+" \n"+fname2
    savefig_name="./plots/diff.png"
    if tdump1 == tdump2:
      title="t = {} \n".format(tdump1)+title
      savefig_name=savefig_name.replace(".png","_t{:013d}.png".format(tdump1))
    else:
      title="t1 = {}, t2 ={} \n".format(tdump1,tdump2)+title
    plt.suptitle(title)
    plt.savefig(savefig_name,bbox_inches="tight")
    print("figure saved at "+savefig_name)


def plot_mdot_phib_t(num_zones,start_zone,dirtag,r_list=[2]):
    matplotlib_settings()
    GHOST=False #True
    first_dir=sorted(glob.glob("../data/"+dirtag+"/*/"))[0]
    fname=glob.glob(first_dir+"/*.00000.rhdf")[0]
    dump = pyharm.load_dump(fname,ghost_zones=GHOST)
    #if num_zones>1: #2: #
    #  start_zone=dump["nzone"]-1 # innermost annulus
    #  skip=2*(dump["nzone"]-1) #dump["nzone"] #
    #else:
    #  start_zone=0
    skip=1

    length = int(np.ceil(num_zones/(2.*(dump["nzone"]-1))))#int(num_zones/skip)+(skip!=1)
      
    fig,ax = plt.subplots(2,length,figsize=(24,12),sharey='row',sharex='col')
    if length==1:
      ax=np.array([[ax[0]],[ax[1]]])
    var="phi_b_per"
    colors=plt.cm.rainbow(np.linspace(1,0,len(r_list)))
    #logfmt = matplotlib.ticker.LogFormatterExponent(base=10.0, labelOnlyBase=True)
    #tmins=np.zeros(length)
    #tmaxs=np.zeros(length)
    tmin = None; tmax = 0; rndtrp=0; last_phi_b=0 # initialization
    
    for r_num, r_chosen in enumerate(r_list):
        t=[]; phi_b=[]; Mdots=[]
        for i,zone in enumerate(range(start_zone,start_zone+num_zones,skip)):
            # configure
            dirname="../data/"+dirtag+"*{:05d}/".format(zone)
            fnames=sorted(glob.glob(dirname+"*.phdf")) #*.rhdf")) #
            dump = pyharm.load_dump(fnames[0],ghost_zones=GHOST)
            if i == 0: rndtrp_norm = dump["iteration"]//2
            if (rndtrp != dump["iteration"]//2 - rndtrp_norm or i == num_zones-1) and len(t)>0 and rndtrp<length:
                # if roundtrip # changes, plot and restart
                xlist = np.linspace(0,1,len(phi_b))
                if rndtrp == length-2: label = "r={}".format(r_chosen)
                else: label=None
                ax[0,rndtrp].semilogy(xlist,phi_b, color=colors[r_num],marker='.',label=label) #,ls='None')
                ax[1,rndtrp].semilogy(xlist,Mdots, color=colors[r_num],marker='.') #,ls='None')
                if r_chosen == 2: last_phi_b = phi_b[-1]
                t=[]; phi_b=[]; Mdots=[]; rndtrp+=1


        
        #if USE_HST_FILE:
          # NOT SUPPORTED NOW read from hst file. this is not accurate since it assumes EH is at i=5
          #diag=glob.glob(dirname+"*.hst")[0]
          #ana = AnaResults(diag)
          #t,phi_b=(ana.get_result('t', var))
          #_,Mdots=(ana.get_result('t', "mdot"))
        #else:
          # calculate directly from dump file
            if r_chosen < dump["r_out"] and r_chosen > dump["r_in"]:
                print(zone)
                for fname in fnames:
                    dump=pyharm.load_dump(fname,ghost_zones=GHOST)
                    Mdot=calc_Mdot(dump,EH=r_chosen)
                    Mdots+=[Mdot]
                    phi_b+=[calc_Phi_b(dump,EH=r_chosen)/np.sqrt(Mdot)]
                    t+=[io.get_dump_time(fname)]
                    rndtrp = int(dump["iteration"])//2 - rndtrp_norm # roundtrip
                    del dump
                
                if tmin is None: tmin = t[0] # initialize tmin only if its not initialized
                if t[-1] > tmax: tmax=t[-1] # keep updating

        
        #if length>1:
            #for j in range(2):
                #ax[j,i].axvline(t[-1],color='gray')
                #ax[j,i].axvline(t[0],color='gray')
            #tmins[i]=t[0]
            #tmaxs[i]=t[-1]
        #ax[i].xaxis.set_minor_formatter(logfmt)

    if length>1:
        x_ticks=[tmin, tmax]#[tmins[0],tmaxs[-1]] #(tmins+tmaxs)/2
        for i in range(length):
            if 1:
                #ax[1,i].axhline(0.5,alpha=0.5)
                ax[0,i].axhline(50,lw=5,alpha=0.5)
                ax[1,i].set_xticks([])
            for j in range(2):
                #ax[j,i].set_xticks(x_ticks)
                #ax[j,i].set_xticklabels(['{:.3g}'.format(x) for x in x_ticks])
                ax[j,i].set_xlim(ax[j,i].get_xlim())#[tmins[i],tmaxs[i]])
                ax[j,i].axis('tight')

                if i==0:
                  ax[j,i].spines['right'].set_visible(False)
                elif i==length-1:
                  ax[j,i].spines['left'].set_visible(False)
                else:
                  ax[j,i].spines['left'].set_visible(False)
                  ax[j,i].spines['right'].set_visible(False)
          
    #ax.set_ylim([1e-1,20])
    ax[0,0].set_ylabel(pyharm.pretty(var))
    ax[1,0].set_ylabel(pyharm.pretty('Mdot'))

    plt.subplots_adjust(wspace=0.005)
    ax[1,length//2].set_xlabel(r"$t$")
    ax[0,0].set_ylim([10,1e4])
    ax[1,0].set_ylim([1e-5,10])
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    #ax.grid(True)
    ax[0,length-2].legend(bbox_to_anchor=(1.1, 0.5), loc="center left", borderaxespad=0)

    # save
    save_path="./plots/"+dirtag
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    savefig_name="/plot_phib_mdot.png"
    #if USE_HST_FILE:
      #savefig_name=savefig_name.replace(".png","_hst.png")
    fig.suptitle(dirtag.replace("bondi_multizone_","")+r', t={:.3g}-{:.3g} last $\phi_b$={:.5g}'.format(tmin, tmax,last_phi_b))
    fig.savefig(save_path+savefig_name,bbox_inches='tight')
    fig.savefig("./plots/"+savefig_name,bbox_inches='tight')

def plot_figures_one_run(dirtag,run,fxn=floors):
  # plot plots from pyharm.plots.figures for one annulus run
  matplotlib_settings()
  fnames=sorted(glob.glob("../data/"+dirtag+"/bondi_multizone_{:05d}".format(run)+"/*.phdf"))
  diag=None
  save_path="../data/"+dirtag+"/frames_"+fxn.__name__
  os.makedirs(save_path,exist_ok=True)
  print_len=len(fnames)//10
  for i,fname in enumerate(fnames):
    if i%print_len==0: print("plotting "+fname)
    plotrc={'log_r':True, 'window':None} # 'log': True
    dump=pyharm.load_dump(fname,ghost_zones=False) #True)
    fig=plt.figure(figsize=(36,12))
    tdump = io.get_dump_time(fname)
    fxn(fig, dump, diag, plotrc)
    plt.savefig(save_path+"/frames_%013d.png"%int(tdump),bbox_inches='tight')
    plt.close()

def calc_tot_time(dirtag):
  dirs=sorted(glob.glob("../data/"+dirtag+"/bondi_multizone_*"))
  tot_time=0
  for i,dr in enumerate(dirs):
    run = dr.split('/')[-1][-5:]
    fname=glob.glob("../logs/"+dirtag+"/log_multizone"+run+"_out")[0]

    f = open(fname, 'r')
    lines = f.readlines()
    for line in lines[::-1]:
      if "walltime used" in line:
        time=float(line.split(' ')[-1])
        break
    tot_time+=time

  tot_time=(tot_time*u.s).to('h')

  print(dirtag.replace("bondi_multizone_","")+": ~#"+ run+", the total time is {:g}".format(tot_time))
  return tot_time
  

def plot_quantities(quantities,oz=False,nzone=1,sz=0,dirtag='',skip=1,show_until=None):
    matplotlib_settings()
    GHOST=True
    NORM=False #True # normalize by sound speed
    ONEZONE= oz# True #False # 
    KERR=False #True #
    start_zone=sz
    
    colors=plt.cm.gnuplot(np.linspace(0.9,0.3,int(nzone/skip)+(skip!=1)))
    
    # configure
    dirname=glob.glob("../data/"+dirtag+"*0*/")[0]
    fname1=sorted(glob.glob(dirname+"*.rhdf"))[0]
    dump = pyharm.load_dump(fname1,ghost_zones=GHOST)
    if dump["solver"]=="none" or dump["bz"]==0:
        # if HD, use phdf file
        BFIELD=False
        filetype="phdf"
    else:
        # if MHD, use rhdf file
        BFIELD=True
        filetype="rhdf"

    if not BFIELD:
        if "beta" in quantities: quantities.remove("beta")
        if "bsq" in quantities: quantities.remove("bsq")
        if "b" in quantities: quantities.remove("b")

    axes=[]
    figs=[]
    
    for quantity in quantities:
      fig,ax = plt.subplots(1,1,figsize=(16,12))
      figs+=[fig]
      axes+=[ax]

    for i,zone in enumerate(range(start_zone,start_zone+nzone,skip)):
        dirname="../data/"+dirtag+"*{:05d}/".format(zone)
        
        # load files
        files=sorted(glob.glob(dirname+"*."+filetype))
        fname1=files[0]
        fname2=files[-1]
        dump1 = pyharm.load_dump(fname1,ghost_zones=GHOST)
        dump2 = pyharm.load_dump(fname2,ghost_zones=GHOST)
        
        if nzone<20:
            show_label=True
        else:
            if i==0 or i==nzone-1:
                show_label=True
            else:
                show_label=False
        
        for j,quantity in enumerate(quantities):
            # plot
            ax=axes[j]
            plot_single_dump(ax,dump1,quantity,colors[i],NORM,show_label=show_label)
            plot_single_dump(ax,dump2,quantity,colors[i],NORM,lw=4,ls='--',alpha=0.5,show_label=show_label)
            
    
    save_path="./plots/"+dirtag
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    for i, quantity in enumerate(quantities):
        ax=axes[i]

        # show bondi analytical solutions
        if 1: #"gizmo" not in dirtag:
            rarr= np.logspace(np.log10(2),np.log10(ax.get_xlim()[1]),20) #np.array([2,3,5,8,10,16])#
            analytic_sol=bondi.get_quantity_for_rarr(rarr,quantity,dump["rs"])
            if analytic_sol is not None and NORM==False:
                ax.plot(rarr,analytic_sol,'k:',label='bondi analytic')
        if dump1["use_gizmo"]=="true" and "txt" in dump1["datfn"]:
            try:
                dat_gz=np.loadtxt(dump1["datfn"])#"../data/gizmo/first_test/dat.txt")
            except:
                dat_gz=np.loadtxt("../data/gizmo/031623_100Myr/dat.txt")
            r_gizmo=dat_gz[:,0]
            rho_gizmo = dat_gz[:,1]
            T_gizmo=dat_gz[:,2]
            if quantity=='rho':
                ax.plot(r_gizmo,rho_gizmo,'b:',label='GIZMO')
            elif quantity=="u":
                ax.plot(r_gizmo,T_gizmo*rho_gizmo*3/2,'b:',label="GIZMO")
            elif quantity=="T":
                ax.plot(r_gizmo,T_gizmo,'b:',label="GIZMO")

        # show Kepler velocity
        if "u^phi" in quantity or "u^th" in quantity:
            rarr= np.logspace(np.log10(2),np.log10(ax.get_xlim()[1]),20)
            ax.plot(rarr,np.power(rarr,-3./2),'r:',label=r'$\Omega_K$')

        # show density scalings
        if "rho" in quantity: # and "rot" in dirtag:
            rarr= np.logspace(np.log10(2),np.log10(dump["rs"]**2.),21)
            factor=1e-5#7e-7
            ax.plot(rarr,np.power(rarr/1e3,-1)*factor,'g-',alpha=0.3,lw=10,label=r'$r^{-1}$')
            ax.plot(rarr,np.power(rarr/1e3,-1/2)*factor,'g:',alpha=0.3,lw=10,label=r'$r^{-1/2}$')

        # labels, scales, legend, title
        ax.axvline(2)
        if "superexp" in dirtag and ("u^r" in quantity or "U1" in quantity):
            ax.axvline(4e6)
        ax.set_xlabel(pyharm.pretty('r'))
        if quantity!="g_00" and quantity!="u^r_cs":
            ax.set_ylabel(pyharm.pretty(quantity))
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_title("{} runs {} - {}".format(dirtag.replace("bondi_multizone_",""),start_zone,zone))

        # save figures
        savefig_name="/multizone_check_"+quantity+".png"
        print("saved at "+save_path+savefig_name)
        figs[i].savefig(save_path+savefig_name,bbox_inches='tight')
        figs[i].savefig("./plots/"+savefig_name,bbox_inches='tight') # make a copy here as well

def same_time_comparison(oz_dir,mz_dir):
  fig,ax=plt.subplots(1,1,figsize=(8,6))
  quantity='rho' #'T' #'beta' #
  NORM=False
  show_label=True

  oz_last = sorted(glob.glob("../data/"+oz_dir+"bondi_multizone_*/"))[-1]
  oz_f = sorted(glob.glob(oz_last+"*.phdf"))[-1]
  oz_d = pyharm.load_dump(oz_f)
  oz_time = pyharm.io.get_dump_time(oz_f)

  mz_dirs = sorted(glob.glob("../data/"+mz_dir+"bondi_multizone_*/*final.phdf"))[::-1] # search backwards
  for i, mz_dir in enumerate(mz_dirs):
    mz_f = mz_dir #glob.glob(mz_dir+"*final.phdf")[0]
    mz_time = pyharm.io.get_dump_time(mz_f)
    if mz_time < oz_time:
      print("stop search at ",i, mz_f, mz_time, oz_time)
      break

  i_stop = i
  mz_d = pyharm.load_dump(mz_dirs[i_stop])
  zone_num = int(mz_d["nzone"]-np.log(mz_d["r_in"])/np.log(int(mz_d["base"]))-1)
  out_to_in = np.power(-1,mz_d["iteration"]+1)
  for i in range(mz_d["nzone"]):
    i_other = i_stop+(zone_num-i)*out_to_in # other zones to cover all radii
    print(i, mz_dirs[i_other])
    mz_d = pyharm.load_dump(mz_dirs[i_other])
    plot_single_dump(ax,mz_d,quantity,'b',NORM,ls='--',alpha=0.5,show_label=show_label)

  plot_single_dump(ax,oz_d,quantity,'k',NORM,show_label=show_label)

  ax.set_xlabel(pyharm.pretty('r'))
  ax.set_ylabel(pyharm.pretty(quantity))
  ax.set_xscale('log'); ax.set_yscale('log')
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

  savefig_name="same_time_comparison.png"
  plt.savefig("./plots/"+savefig_name,bbox_inches='tight') # make a copy here as well

def _main():
  ONEZONE=False # True #
  NZONE = 8
  STARTZONE=46711 #
  SKIP=1
  #dirtag="production_runs/bondi_bz2e-8_1e8_96/"
  #dirtag="061023_bflux0_n8_gamma3/"
  #dirtag="061323_diffinit_better/"
  #dirtag="061423_in2out_rst/"
  #dirtag="061623_ozrst/"
  #dirtag="062023_0.02tff/"#_ur0/"
  #dirtag="062623_constinit_difftchar/"
  #dirtag="bondi_multizone_050523_bflux0_1e-4_128^3_n3_noshort/"
  #dirtag="062623_n3_1/"
  #dirtag="062623_0.04tff/"
  #dirtag="062723_b3_0.2tff/"
  #dirtag="062723_0.01tff_8_6/"
  dirtag="071023_beta01/"
  #dirtag="071123_b3n16_beta01/"
  #dirtag="071423_beta03/"
  #dirtag="071723_diffbinit/"


  #ozdir="bondi_multizone_050423_onezone_bflux0_1e-8_2d_n4/"
  #mzdir="bondi_multizone_050423_bflux0_1e-8_2d_n4/"
  ozdir="bondi_multizone_050123_onezone_bflux0_1e-4_64^3/"
  mzdir="bondi_multizone_042723_bflux0_1e-4_64^3/"
  #same_time_comparison(ozdir,mzdir)

  quantities=["rho","u","T","b","beta","abs_u^r","abs_u^phi","abs_u^th","K","u^r","u^phi"] #,"abs_U1"] #
  plot_quantities(quantities,ONEZONE,NZONE,STARTZONE,dirtag=dirtag,show_until=None)
  #plot_mdot_r(ONEZONE,NZONE,STARTZONE,dirtag=dirtag,avg=False)
  #plot_mdot_phib_t(NZONE,STARTZONE,dirtag,r_list=[2,25,200,1500,12000])
  #calc_tot_time(dirtag)
  #plot_figures_one_run(dirtag,35,floors) # e_ratio, prims

  #plots_together(dirtag)
  #check_bondi_multizone('JE3Fl',ONEZONE,NZONE,STARTZONE,dirtag=dirtag,skip=SKIP) # angular momentum density
  #plot_mdot_t()
  #plot_one_run(STARTZONE,'Mdot',dirtag=dirtag)
  #plot_grid_3d(dirtag)

  '''
  #fname1="../data/"+dirtag+"/bondi_multizone_00001/resize_restart_kharma.out1.00000.rhdf"
  fname1="../data/"+dirtag+"/bondi_multizone_00000/bondi.out1.final.rhdf"
  #fname2="../data/"+dirtag+"/bondi_multizone_00002/resize_restart_kharma.out1.now.rhdf"
  fname2="../data/bondi_multizone_033023_bflux_1e-3_lin/bondi_multizone_00000/bondi.out1.final.rhdf"
  #for i in range(6):
  #  fname1="../data/"+dirtag+"/bondi_multizone_00000/bondi.out1.{:05d}.rhdf".format(i)
  #  fname2="../data/bondi_multizone_040123_onezone_1e-3_lin/bondi_multizone_00000/bondi.out1.{:05d}.rhdf".format(i)
  plot_diffs(fname1,fname2,"b")

  #dirtag="bondi_multizone_022723_jit0.3_new_coord/"
  dump=pyharm.load_dump("../data/"+dirtag+"/bondi_multizone_00001/resize_restart_kharma.out1.00000.rhdf")
  fig=plt.figure(figsize=(8,6))
  #print(fft_comp(dump,'u^r').shape)
  #for i,r in enumerate(dump["r1d"]):
  i=30 #11 # somewhere middle
  r=dump["r1d"][i]
  #plt.plot(dump["th1d"],th_average(dump,"rho",i),label=r'$r={:.3g}$'.format(r))
  quantity="bsq" #"rho" #
  plt.plot(dump["phi1d"],phi_average(dump,quantity,i),label=r'$r={:.3g}$'.format(r))
  #plt.xlabel(r'$\theta$'); plt.ylabel(r'$\rho$');
  plt.xlabel(r'$\phi$'); plt.ylabel(pyharm.pretty(quantity));
  plt.legend()
  plt.savefig("./temp.png",bbox_inches='tight')

  '''

if __name__=='__main__':
  _main()

