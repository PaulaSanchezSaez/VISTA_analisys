#codigo que genera las curvas de luz, en archivos fits para cada objeto del catalogo de los colaboradores, a partir de flujos, y guarda tanto los flujos como las magnitudes AB
import numpy as np
#import pyfits as pf
import astropy.io.fits as pf
import math
import os
import time
from scipy.spatial import cKDTree
from astropy import units as u
from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
from astropy.coordinates.matching import match_coordinates_sky as match_coord
from multiprocessing import Pool

start = time.time()

filt='Z'

field = 'ECDFS' # ECDFS, ELAIS_S1, or XMM_LSS

cat_path = '../cat_filt/'+field+'/'+filt+'/'

lc_path = '../final_catalog_lc/'+field+'/'+filt+'/'

defa='../default_files'




dra=0.015 #size of the seccion eliminated from the images in ra
ddec=0.015 #size of the seccion eliminated from the images in dec


###################################################################
#cargamos el documento que tiene los nonmbres de los archivos, el JD, y errzp

lista_cat=np.loadtxt('../stat/datos_cat_'+field+'_'+filt+'.txt',dtype='str').transpose()
cat=lista_cat[0]
jd=lista_cat[2].astype(np.float)
zp=lista_cat[3].astype(np.float)
errzp=lista_cat[4].astype(np.float)
tim=lista_cat[5].astype(np.float)
extinct=lista_cat[6].astype(np.float)
air=lista_cat[7].astype(np.float)




###################################################################
#cargamos el catalogo de AGNs detectados a partir de XR

agn=pf.open('../catalogs/final_catalog_'+field+'.fits')
dat=agn[1].data
ra1=dat['RA']
dec1=dat['DEC']
redshift=dat['REDSHIFT']
org_cat=dat['catalog']


n=int((float(len(ra1))/3.0))

'''

ra1=ra1[0:n]
dec1=dec1[0:n]
redshift=redshift[0:n]
org_cat=org_cat[0:n]



ra1=ra1[n:2*n]
dec1=dec1[n:2*n]
redshift=redshift[n:2*n]
org_cat=org_cat[n:2*n]



ra1=ra1[2*n:]
dec1=dec1[2*n:]
redshift=redshift[2*n:]
org_cat=org_cat[2*n:]

'''

#cat_coord=ICRS(ra1,dec1,unit=(u.degree,u.degree))

cat_coord = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)


#se carga la lista de los objetos que vamos a guardar:
#df=np.empty((len(ra),0)).tolist()


###################################################################
#iteramos por imagen:

def cat_match(cat,jd,zp,errzp,tim,extinct,air,ncores,nthread):

    ncat=int(len(jd)/ncores)


    if nthread<(ncores-1):
        cat=cat[nthread*ncat:(nthread+1)*ncat]
        jd=jd[nthread*ncat:(nthread+1)*ncat]
        zp=zp[nthread*ncat:(nthread+1)*ncat]
        errzp=errzp[nthread*ncat:(nthread+1)*ncat]
        tim=tim[nthread*ncat:(nthread+1)*ncat]
        extinct=extinct[nthread*ncat:(nthread+1)*ncat]
        air=air[nthread*ncat:(nthread+1)*ncat]


    else:
        cat=cat[nthread*ncat:-1]
        jd=jd[nthread*ncat:-1]
        zp=zp[nthread*ncat:-1]
        errzp=errzp[nthread*ncat:-1]
        tim=tim[nthread*ncat:-1]
        extinct=extinct[nthread*ncat:-1]
        air=air[nthread*ncat:-1]




    d=np.empty((len(ra1),0)).tolist()
    dchips=np.empty((len(ra1),0)).tolist()
    largo=len(jd)


    for j, item in enumerate(jd):

        if ("_16" not in cat[j]) and os.path.isfile(cat_path+cat[j]+"_"+filt+".calib.fits"):
            arch=pf.open(cat_path+cat[j]+"_"+filt+".calib.fits")
            datos=arch[1].data
            head=arch[0].header
            fit_rms=head['fit_rms']
            chi_red=head['chi_red']
            num_stars=head['num_stars']
            alpha=datos['ALPHA_J2000']
            delta=datos['DELTA_J2000']
            flags=datos['FLAGS']
            lar_ra=len(alpha)
            #im_coord=np.dstack((alpha,delta))
            #im_coord=im_coord[0]
            #kdt = cKDTree(im_coord)
            #dist,ind=kdt.query(cat_coord, k=1, distance_upper_bound= 0.00007)


            min_ra=np.amin(alpha)
            max_ra=np.amax(alpha)
            min_dec=np.amin(delta)
            max_dec=np.amax(delta)

            #img_coord=ICRS(alpha,delta,unit=(u.degree,u.degree))
            img_coord = SkyCoord(ra=alpha*u.degree, dec=delta*u.degree)
            ind,ang,dis=match_coord(cat_coord,img_coord,nthneighbor=1)
            ang0=np.array(ang)
            #n=np.where(ang0<0.00008)
            n=np.where(ang0<0.000277778)
            #import pdb;pdb.set_trace()
            n=n[0] #posicion en el catalogo de los objetos
            if len(n)>0:
                pos=ind[n] #posicion en la imagen
                alpha0=alpha[pos]
                delta0=delta[pos]
                flags0=flags[pos]
                magap=datos['MAG_2'][pos]
                #magaper=magap-extinct[j]*(air[j]-1.0)
                magaper=magap
                erraper=datos['MAGERR_2'][pos]

                f_ap=10.0**(-0.4*(magaper+48.6))
                errf_ap=0.4*np.log(10)*(10.0**(-0.4*(magaper+48.6)))*erraper


                for i, item in enumerate(pos):
                    if (flags0[i]==0) and (alpha0[i]>(min_ra+dra)) and (alpha0[i]<(max_ra-dra)) and (delta0[i]>(min_dec+ddec)) and (delta0[i]<(max_dec-ddec)):
                        e=[jd[j],magaper[i],erraper[i],f_ap[i],errf_ap[i]]
                        d[n[i]].append(e)
                        dchips[n[i]].append([cat[j],fit_rms,chi_red,num_stars])
                        #print d[ind[n][i]]

            #print "done image %s" % (cat[j])
            print "done image %s    %f     %d" % (cat[j],jd[j], (largo-j))
            arch.close()
            del arch
    print len(d)
    return (d,dchips)

pool=Pool(processes=8)
p0=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,0))
p1=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,1))
p2=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,2))
p3=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,3))
p4=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,4))
p5=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,5))
p6=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,6))
p7=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,8,7))

'''
p8=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,8))
p9=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,9))
p10=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,10))
p11=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,11))
p12=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,12))
p13=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,13))
p14=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,14))
p15=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,15))
p16=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,16))
p17=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,17))
p18=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,18))
p19=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,19))
p20=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,20))
p21=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,21))
p22=pool.apply_async(cat_match,args=(cat,jd,zp,errzp,tim,extinct,air,23,22))
'''










d0=p0.get()
d1=p1.get()
d2=p2.get()
d3=p3.get()
d4=p4.get()
d5=p5.get()
d6=p6.get()
d7=p7.get()
'''
d8=p8.get()
d9=p9.get()
d10=p10.get()
d11=p11.get()
d12=p12.get()
d13=p13.get()
d14=p14.get()
d15=p15.get()
d16=p16.get()
d17=p17.get()
d18=p18.get()
d19=p19.get()
d20=p20.get()
d21=p21.get()
d22=p22.get()
'''

pool.close()
pool.join()

d=np.empty((len(ra1),0)).tolist()
dchips=np.empty((len(ra1),0)).tolist()

for i, item in enumerate(ra1):
    d[i]=d0[0][i]+d1[0][i]+d2[0][i]+d3[0][i]+d4[0][i]+d5[0][i]+d6[0][i]+d7[0][i]#+d8[0][i]+d9[0][i]+d10[0][i]+d11[0][i]+d12[0][i]+d13[0][i]+d14[0][i]+d15[0][i]+d16[0][i]+d17[0][i]+d18[0][i]+d19[0][i]
    dchips[i]=d0[1][i]+d1[1][i]+d2[1][i]+d3[1][i]+d4[1][i]+d5[1][i]+d6[1][i]+d7[1][i]#+d8[1][i]+d9[1][i]+d10[1][i]+d11[1][i]+d12[1][i]+d13[1][i]+d14[1][i]+d15[1][i]+d16[1][i]+d17[1][i]+d18[1][i]+d19[1][i]


print len(d)



#del d0[:],d1[:],d2[:],d3[:],d4[:],d5[:],d6[:],d7[:]#,d8[:],d9[:],d10[:],d11[:],d12[:],d13[:],d14[:],d15[:],d16[:],d17[:],d18[:],d19[:],d20[:],d21[:],d22[:]
###################################################################



print "match done"



print "start light curve generation"


def gen_lc(d,dchips,ra,dec,redshift,org_cat):
    #function to create every light curve
    f=np.array(d)
    fchips=np.array(dchips)
    #print "arreglo",f

    if f.any():
        f=f[f[:,0].argsort()]
        fchips=fchips[f[:,0].argsort()]
        jdo=f[:,0]
        #print jdo
        if (len(jdo)>2):
            c1=pf.Column(name='JD',format='D',array=jdo)
            c2=pf.Column(name='MAG_2',format='D',array=f[:,1])
            c3=pf.Column(name='MAGERR_2',format='D',array=f[:,2])
            c4=pf.Column(name='FLUX_2',format='D',array=f[:,3])
            c5=pf.Column(name='FLUXERR_2',format='D',array=f[:,4])

            c6=pf.Column(name='catalog',format='32A',array=fchips[:,0])
            c7=pf.Column(name='fit_rms',format='D',array=fchips[:,1])
            c8=pf.Column(name='chi_red',format='D',array=fchips[:,2])
            c9=pf.Column(name='num_stars',format='D',array=fchips[:,3])



            #coldef=pf.ColDefs([c1,c2,c3,c4,c5,c6,c7,c8,c9])
            #thdu=pf.new_table(coldef)
            thdu=pf.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6,c7,c8,c9])

            ras="%.6f" % ra
            decs="%.6f" % dec

            thdu.writeto(ras+"_"+decs+"_"+filt+".fits")
            arch1=pf.open(ras+"_"+decs+"_"+filt+".fits",mode='update')
            head=arch1[0].header
            head['RA']=    ra
            head['DEC']=    dec
            head['FILTER']= filt
            head['REDSHIFT']=redshift
            head['ORG_CAT']=org_cat



            #print "hr", hr[j]
            arch1.flush()
            arch1.close()
            del arch1
            cmd="mv "+ras+"_"+decs+"_"+filt+".fits "+lc_path
            os.system(cmd)
            print "done obj = %s \t" % (ras+"_"+decs+"_"+filt+".fits")
    return (ra)


def run_gen_lc(args):
#function necessary to use of pool.map when running gen_lc with different arguments
   return gen_lc(*args)




args_run_gen_lc=[]
#the array with the arguments is generatedm this is necessary to use pool.map
for i, item in enumerate(ra1):
    args_run_gen_lc.append((d[i],dchips[i],ra1[i],dec1[i],redshift[i],org_cat[i]))


print len(args_run_gen_lc)
#gen_lc is run
pool = Pool(processes=7)

results = pool.map(run_gen_lc,args_run_gen_lc)

pool.close()
pool.join()




###################################################################

elapsed = (time.time() - start)
print elapsed
