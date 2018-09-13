#code to filter the lc, eliminating bad epochs.
#import scipy.stats as stat
import numpy as np
from scipy import interpolate
import scipy.signal as sgl
import pyfits as pf
import os
from decimal import Decimal
from astropy import stats
from scipy.integrate import quad
import glob



filt='Z'

field = 'ECDFS' # ECDFS, ELAIS_S1, or XMM_LSS

cat_path = '../cat_filt/'+field+'/'+filt+'/'

lc_path = '../final_catalog_lc/'+field+'/'+filt+'/'




cmd='rm '+lc_path+'agn_*fits'
os.system(cmd)

agn=sorted(glob.glob(lc_path+'*fits'))
agn=np.core.defchararray.replace(agn,lc_path,'')



###################################################################
#cargamos la lista de AGNs

#agn=np.loadtxt('agn_Y.txt',dtype='str').transpose()

color=['b*','g*','r*']
colorf=['ro','bo','go']
aper=['2','3','7']

print color
###################################################################
#hacemos el filtro

for i, item in enumerate(agn):
    #for i in range(1,30):

    print agn[i]
    arch=pf.open(lc_path+agn[i])
    datos=arch[1].data
    head0=arch[0].header
    mag2=datos['MAG_2']
    errmag2=datos['MAGERR_2']
    flux2=datos['FLUX_2']
    errflux2=datos['FLUXERR_2']
    jd0=datos['JD']
    mag=mag2
    m=np.where(mag2<30)
    jd=jd0[m]
    mag2=mag2[m]
    err2=errmag2[m]
    flux2=flux2[m]
    errflux2=errflux2[m]
    errmean2=np.median(err2)
    l=np.where((err2<2*errmean2) & (errflux2>0))
    jd=jd[l]
    mag2=mag2[l]
    err2=err2[l]
    flux2=flux2[l]
    errflux2=errflux2[l]

    #se hace un mean_filter de tres puntos
    #mag_filt=sgl.medfilt(mag2)
    #se ajusta un polinomio de orden 5.
    if len(jd)>0:
        coefficients = np.polyfit(jd, mag2, 5)
        polynomial = np.poly1d(coefficients)
        pol=polynomial(jd)
        #se calcula el rms
        rm=(mag2-pol)**2
        rms=np.sum(rm)
        rms=rms/(len(jd)-1)
        #se hace el filtro por 3 sigmas
        dist=np.abs(mag2-pol)
        sigma=np.sqrt(rms+err2**2)
        n=np.where((dist/sigma)<=2)
        dia=jd[n]
        magf2=mag2[n]
        errf2=err2[n]
        fluxf2=flux2[n]
        errfluxf2=errflux2[n]
        dif=len(jd0)-len(dia)
        print dif


    #se guardan los fits
    c1=pf.Column(name='JD',format='D',array=dia)
    c2=pf.Column(name='MAG_2',format='D',array=magf2)
    c3=pf.Column(name='MAGERR_2',format='D',array=errf2)
    c4=pf.Column(name='FLUX_2',format='D',array=fluxf2)
    c5=pf.Column(name='FLUXERR_2',format='D',array=errfluxf2)
    coldef=pf.ColDefs([c1,c2,c3,c4,c5])
    thdu=pf.new_table(coldef)
    thdu.writeto(lc_path+'agn_'+agn[i])
    arch1=pf.open(lc_path+'agn_'+agn[i],mode='update')
    head=arch1[0].header
    head=arch1[0].header

    head['RA']=    head0['RA']
    head['DEC']=    head0['DEC']
    head['REDSHIFT']=head0['REDSHIFT']
    head['ORG_CAT']=head0['ORG_CAT']
    head['DELPOINT']=dif
    head['FILTER']=filt
    head['T_RANGE']=dia[-1]-dia[0]

    print (dia[-1]-dia[0])
    arch1.flush()

    cmd='rm '+lc_path+'agn_'+agn[i]
    print "dia= %d" % (len(dia))
    if (len(dia)<3): os.system(cmd)
    arch.close()
