import numpy as np
import pyfits as pf
import math
import os
import time
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS
from astropy.coordinates.matching import match_coordinates_sky as match_coord
import glob

###########################################################################
#define the filter band and the location of the catalogs
filt='Ks'

field = 'ECDFS' # ECDFS, ELAIS_S1, or XMM_LSS

fold_catalogs = '../cat_filt/'+field+'/'+filt+'/'

defa='../default_files'

list_cat = sorted(glob.glob(fold_catalogs+'*cat.fits'))

list_cat_name = np.core.defchararray.replace(list_cat,fold_catalogs,'')

namelist=np.core.defchararray.replace(list_cat_name,'.cat.fits','')

###########################################################################
#read the UVISTA catalog


catq=pf.open('../catalogs/VIDEO-'+field+'_2016-04-14_fullcat_errfix_below25.fits')
dat=catq[1].data
raq=dat['ALPHA_J2000']
decq=dat['DELTA_J2000']
if filt!='Ks':
    qmag=dat[filt+'_MAG_APER_2']
    qerr=dat[filt+'_MAGERR_APER_2']
else:
    qmag=dat['K_MAG_APER_2']
    qerr=dat['K_MAGERR_APER_2']

#convert units
#cat_coord=ICRS(raq,decq,unit=(u.degree,u.degree))
cat_coord = SkyCoord(ra=raq*u.degree, dec=decq*u.degree)
#print cat_coord
#SkyCoord(ICRS, ra=ra, dec=dec


############################################################################
#read the list of the catalogs of the images

#os.system('ls v2*cat.fits > cat_list.txt')

#os.system('rm v2*calib.fits >& /dev/null')

#rlist=np.loadtxt('cat_list.txt',dtype='str').transpose()

#namelist=np.core.defchararray.replace(rlist,'.cat.fits','')

for i, item in enumerate(list_cat_name):
#for i in range(61,62):
    print list_cat_name[i]

    if os.path.isfile(fold_catalogs+namelist[i]+'.calib.fits'):
        cmd='rm '+fold_catalogs+namelist[i]+'.calib.fits'
        os.system(cmd)
        print "file %s deleted" % (fold_catalogs+namelist[i]+'.calib.fits')

    arch=pf.open(list_cat[i])
    dat=arch[1].data
    alpha=dat['ALPHA_J2000']
    delta=dat['DELTA_J2000']
    flags=dat['FLAGS']
    flagsf=flags
    clas=dat['CLASS_STAR']
    fwhm=dat['FWHM_IMAGE']
    maginst=dat['MAG_APER'][:,1]
    magQi=maginst
    #print np.amin(maginst),np.amax(maginst)
    #errinst=dat['MAGERR_APER'][:,1]
    fluxoo=dat['FLUX_APER'][:,1]
    errinst=(2.5/np.log(10))*(dat['FLUXERR_2']/fluxoo)
    errmagQi=errinst
    #img_coord=ICRS(alpha,delta,unit=(u.degree,u.degree))
    img_coord = SkyCoord(ra=alpha*u.degree, dec=delta*u.degree)
    ind,ang,dis=match_coord(cat_coord,img_coord,nthneighbor=1)
    #ang=ang*u.degree
    ang0=np.array(ang)
    #n=np.where(ang0<0.00008)
    n=np.where(ang0<0.000277778)
    n=n[0]
    #print ind
    print np.amin(ang0)
    if len(n)>0:
        pos=ind[n]
        flags=flags[pos]
        class_star=clas[pos]
        maginst=maginst[pos]
        errinst=errinst[pos]
        q_mag=qmag[n]
        q_err=qerr[n]
        fn=np.where((class_star>=0.9) & (flags<=4) & (q_mag>15) & (q_mag<20) & (q_err<0.3) & (maginst<23) & (maginst>15) & (errinst<0.3))
        maginst=maginst[fn]
        magtest=maginst
        print "before clean", len(maginst)
        if len(maginst)>20:
            errinst=errinst[fn]
            #print "err_inst",errinst
            q_mag=(qmag[n])[fn]
            q_err=(qerr[n])[fn]
            #print "q_err", q_err
            sigma=np.sqrt(q_err**2+errinst**2)
            errf=sigma**2
            #print "sigma",sigma
            diff=maginst-q_mag

            #print "diff",diff
            coefficients = np.polyfit(maginst, diff, 1)
            polynomial = np.poly1d(coefficients)
            fit=polynomial(maginst)
            #print "fit",fit
            res1=np.abs(diff-fit)
            rms=(np.sum((diff-fit)**2))/(len(maginst)-1)
            sigma=np.sqrt(sigma**2+rms)
            l1=np.where((res1/sigma)<=3)
            r1=np.where((res1/sigma)>3)
            print "first clean",len(l1[0])

            #print "res1",res1
            diff2=diff
            diff=diff[l1]
            diff_rej=diff2[r1]
            sigma=sigma[l1]
            errf=errf[l1]
            maginst=maginst[l1]
            q_mag=q_mag[l1]
            fita=fit[l1]

            '''
            plt.plot(magtest,diff2,'ro')
            plt.plot(maginst,diff,'k*')
            plt.plot(magtest,fit,'b-')
            plt.show()
            '''

            if len(diff)>10:
                while (len(diff)>10):
                    if len(diff_rej)>0:

                        coefficients = np.polyfit(maginst, diff, 1)
                        polynomial = np.poly1d(coefficients)
                        fit=polynomial(maginst)
                        res2=np.abs(diff-fit)
                        rms=(np.sum((diff-fit)**2))/(len(maginst)-1)
                        sigma=np.sqrt(sigma**2+rms)
                        l2=np.where((res2/sigma)<=3)
                        r2=np.where((res2/sigma)>3)
                        diff2=diff
                        diff=diff[l2]
                        diff_rej=diff2[r2]
                        sigma=sigma[l2]
                        errf=errf[l2]
                        maginst=maginst[l2]
                        q_mag=q_mag[l2]
                    else: break



                coefficients = np.polyfit(maginst, diff, 1)
                fit=polynomial(maginst)

                rms=(np.sum((diff-fit)**2))/(len(maginst)-1)

                chi_red=(np.sum(((diff-fit)**2)/errf))/(len(maginst)-2)
                print np.sqrt(rms), chi_red, len(maginst)



                '''
                plt.plot(maginst,diff,'ro')
                plt.plot(maginst,fit,'b-')
                plt.show()
                '''

                final_pol=polynomial
                final_fit=final_pol(magQi)
                Q_cal=magQi-final_fit
                Q_err=np.sqrt(errmagQi**2+rms)

                #writing the new catalogs:
                c1=pf.Column(name='ALPHA_J2000',format='D',array=alpha)
                c2=pf.Column(name='DELTA_J2000',format='D',array=delta)
                c3=pf.Column(name='MAG_2',format='D',array=Q_cal)
                c4=pf.Column(name='MAGERR_2',format='D',array=Q_err)
                c5=pf.Column(name='FLAGS',format='D',array=flagsf)
                c6=pf.Column(name='CLASS_STAR',format='D',array=clas)
                c7=pf.Column(name='FWHM_IMAGE',format='D',array=fwhm)
                c8=pf.Column(name='MAG_inst',format='D',array=magQi)
                coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8])
                thdu=pf.new_table(coldef)
                thdu.writeto(fold_catalogs+namelist[i]+'.calib.fits')

                arch1=pf.open(fold_catalogs+namelist[i]+'.calib.fits',mode='update')
                head=arch1[0].header
                head['Field']=    field
                head['fit_rms']= np.sqrt(rms)
                head['chi_red']= chi_red
                head['num_stars']= len(maginst)
                arch1.flush()
                del arch1

                print 'catalog', namelist[i] ,'generated ###########'
