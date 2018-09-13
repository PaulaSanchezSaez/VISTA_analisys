#codigo que corre sextractor para una lista de imagenes filtrandolas
import glob
import numpy as np
#import pyfits as pf
import math
import os
#import time
import matplotlib.pyplot as plt
from scipy import misc
from scipy import ndimage
from astropy.io import fits

#start = time.clock()

filt = 'Ks'#filter band

field = 'ELAIS_S1' # ECDFS, ELAIS_S1, or XMM_LSS

fold_stacks = '../VIDEO_data/'+field+'_stacks/'+filt+'/'

fold_catalogs = '../cat_filt_Aug2018/'+field+'/'+filt+'/'

defa='../default_files'

# list with the files's paths

list_stacks = sorted(glob.glob(fold_stacks+'*st.fit'))

list_stacks_conf = np.core.defchararray.replace(list_stacks,'st.fit','st_conf.fit')

list_img_name = [x[-22:-4] for x in list_stacks]

print "total number of images: ", len(list_stacks)

###################################################################

#function to calculate error correction factor
def empty_apertures(image, seg, lado):

    hdufits_ima = fits.open(image)
    imag_data = hdufits_ima[0].data

    ver_max,hor_max = imag_data.shape

    hdufits_seg = fits.open(seg)
    segm_data = hdufits_seg[0].data

    filtered_segm_data=ndimage.gaussian_filter(segm_data, 2)

    segm_mask = (filtered_segm_data > 0)

    #segm_mask = (segm_data > 0)

    n = 0
    empty_ap = []

    matrix_used = np.zeros(np.shape(segm_data))

    while n < 2001:
        random_center = np.random.random_integers(350 ,np.amin(np.array([ver_max,hor_max]))-350, 2)
        #segm_mask = (segm_data > 0)
        if (all((~segm_mask[random_center[0]:random_center[0]+lado, random_center[1]:random_center[1] + lado]).flat)) and (np.sum(matrix_used[random_center[0]:random_center[0]+lado, random_center[1]:random_center[1] + lado])==0):
            counts,npix=get_flux(random_center[0], random_center[1], imag_data, int(lado))
            empty_ap.append(counts)
            n += 1
            matrix_used[random_center[0]:random_center[0]+lado, random_center[1]:random_center[1] + lado ]=22

    empty_ap=np.array(empty_ap)

    #plt.hist(empty_ap,50)
    #plt.savefig(image+'_hist_'+str(lado)+'.png')
    #plt.close('all')

    empty_ap2=empty_ap[np.where((empty_ap<(np.mean(empty_ap)+3*np.std(empty_ap))) & (empty_ap>(np.mean(empty_ap)-3*np.std(empty_ap))))]


    return (np.median(empty_ap2),np.mean(empty_ap2),np.std(empty_ap2),npix)


def get_flux(pos_ver, pos_hor, img, lado):
    #if verbose:  print 'Radius %i' % radius

    counts=np.sum(img[pos_ver : pos_ver + lado, pos_hor : pos_hor + lado])
    numpix=lado**2

    return counts,numpix


def get_error_model(img,seg,apmin,apmax,numap):

    hdufits_ima = fits.open(img)
    imag_data = hdufits_ima[0].data

    ver_max,hor_max = imag_data.shape


    hdufits_seg = fits.open(seg)
    segm_data = hdufits_seg[0].data

    filtered_segm_data=ndimage.gaussian_filter(segm_data, 2)

    segm_mask = (filtered_segm_data > 0)

    matrix_used = np.zeros(np.shape(segm_data))

    pix_count=[]
    n=0

    while n < 3001:
        #segm_mask = (segm_data > 0)
        random_center = np.random.random_integers(300 ,np.amin(np.array([ver_max,hor_max]))-300, 2)
        if (all((~segm_mask[random_center[0], random_center[1]]).flat)) and (np.sum(matrix_used[random_center[0], random_center[1]])==0):

            pix_count.append(imag_data[random_center[0]][random_center[1]])
            n += 1
            matrix_used[random_center[0], random_center[1]]=22

    pix_count=np.array(pix_count)

    #pix_count=pix_count[np.where((pix_count<(np.mean(pix_count)+3.5*np.std(pix_count))) & (pix_count>(np.mean(pix_count)-3.5*np.std(pix_count))))]
    sigma1=np.std(pix_count)

    #plt.hist(pix_count,50)
    #plt.savefig(img+'_sigma1_hist.png')
    #plt.close('all')

    apertures=np.linspace(apmin,apmax,num=numap)
    median=[]
    mean=[]
    std=[]
    npix=[]
    for i in xrange(len(apertures)):
        median0,mean0,std0,npix0=empty_apertures(img, seg, apertures[i])
        median.append(median0)
        mean.append(mean0)
        std.append(std0)
        npix.append(npix0)

    npix=np.array(npix)
    npix=npix[np.where(npix>0)]
    median=np.array(median)[np.where(npix>0)]
    mean=np.array(mean)[np.where(npix>0)]
    std=np.array(std)[np.where(npix>0)]

    N=np.sqrt(npix)

    x=np.log10(N)
    y=np.log10(std)-np.log10(sigma1)
    coefficients=np.polyfit(x,y,1)

    alpha=10**(coefficients[1])
    beta=coefficients[0]

    #plt.plot(N,std,'r*')
    #plt.plot(N,sigma1*alpha*(N**beta),'b-')
    #plt.yscale('log')
    #plt.savefig(img+'_model.png')
    #plt.close('all')

    return (sigma1,alpha,beta)


###################################################################
#cargamos la lista que contiene los nombres de los archivos

#datos=np.loadtxt('lista_datos.txt',dtype='str').transpose()

#definimos el seeing al cual se cambiaran todas las imagenes:

sigma0=2.95 #valor del seeing al cual se convertiran las imagenes

###################################################################
#se genera la lista que contiene las ganancias y saturaciones para cada chip.
GAIN=[3.66,4.25,3.95,4.16,4.24,4.11,3.84,4.22,4.53,3.97,4.62,3.95,5.67,4.78,3.99,4.96]
SAT=[33000.0,32000.0,33000.0,32000.0,24000.0,36000.0,35000.0,33000.0,35000.0,35000.0,37000.0,34000.0,33000.0,35000.0,34000.0,34000.0]
SAT=np.array(SAT)
SAT=SAT*0.9

###################################################################
#se crea un archivo que contiene informacion de cada imagen, como el dia juliano, el zp y el error en el zp.
arch2=open("../stat/datos_cat_"+field+"_"+filt+".txt","w")
arch3=open("../stat/failed_cat_"+field+"_"+filt+".txt","w")
#arch2.writelines("imagen \t chip \t JD \t zp \t errzp \t time \t extinct \t air \n")

###################################################################
#separamos los chips en distintos fits, extraemos la info de los headers y corremos sextractor

for i in xrange(len(list_stacks)):
#for i in range(389,len(datos)):
    print list_img_name[i]
    arch = fits.open(list_stacks[i])
    #arch=pf.open(datos[i]+'_'+filt+'.fits')
    arch_conf=fits.open(list_stacks_conf[i])
    #arch_conf=pf.open(datos[i]+'_conf_'+filt+'.fits')
    head0=arch[0].header
    tim=head0['EXPTIME']
    air_e=head0['HIERARCH ESO TEL AIRM END']
    air_s=head0['HIERARCH ESO TEL AIRM START']
    air=(air_e+air_s)*0.5
    print list_img_name[i], tim
    for j in range(1,17):
        try:
            chip=arch[j].data
            head=arch[j].header

            seeing=head['SEEING']
            print seeing
            if (seeing>0) and (seeing<sigma0) :
                try:
                    sigma=(1/(np.sqrt(8*np.log(2))))*np.sqrt((sigma0**2)-(seeing**2))
                    b_chip=ndimage.gaussian_filter(chip, sigma)

                    #aux=head['HIERARCH ESO DET CHIP PXSPACE']
                    #aux=float(aux)
                    #head['HIERARCH ESO DET CHIP PXSPACE']=0
                    fits.writeto("temp/"+list_img_name[i]+'_'+str(j)+'_'+filt+'.fits',data=b_chip,header=head,output_verify='fix')


                    chip_conf=arch_conf[j].data
                    head_conf=arch_conf[j].header
                    b_chip_conf=ndimage.gaussian_filter(chip_conf, sigma)

                    #auxc=head_conf['HIERARCH ESO DET CHIP PXSPACE']
                    #auxc=float(auxc)
                    #head_conf['HIERARCH ESO DET CHIP PXSPACE']=auxc

                    fits.writeto("temp/"+list_img_name[i]+'_conf_'+str(j)+'_'+filt+'.fits',data=b_chip_conf,header=head_conf,output_verify='fix')

                    filter_image='gauss_2.5_5x5.conv'

                    seeing=seeing*0.339
                    seeing0=sigma0*0.339
                    zp=head['MAGZPT']+2.5*np.log10(tim)

                    errzp=head['MAGZRR']
                    jd=head['MJD-OBS']
                    ext=head['EXTINCT']
                    z=j-1
                    #gain=GAIN[z]#*time
                    gain=4.19
                    sat=SAT[z]
                    #print gain



                    #cmd1="sex "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".fits -c "+defa+"/vistapp.sex -CATALOG_TYPE FITS_1.0 -CATALOG_NAME "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits"+" -PARAMETERS_NAME "+defa+"/paula4.param -FILTER_NAME "+defa+"/"+filter_image+" -MAG_ZEROPOINT "+str(zp)+" -WEIGHT_IMAGE "+"temp/"+list_img_name[i]+"_conf_"+str(j)+"_"+filt+".fits -WEIGHT_TYPE MAP_WEIGHT  -SATUR_LEVEL "+str(sat)+" -GAIN "+str(gain)+" -SEEING_FWHM "+str(seeing0)+" -STARNNW_NAME "+defa+"/default.nnw -CHECKIMAGE_NAME "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg  -CHECKIMAGE_TYPE SEGMENTATION -PHOT_APERTURES 4,6,9,20,30"
                    cmd1="sex "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".fits -c "+defa+"/vistapp.sex -CATALOG_TYPE FITS_1.0 -CATALOG_NAME "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits"+" -PARAMETERS_NAME "+defa+"/paula4.param -FILTER_NAME "+defa+"/"+filter_image+" -MAG_ZEROPOINT "+str(zp)+" -WEIGHT_IMAGE "+"temp/"+list_img_name[i]+"_conf_"+str(j)+"_"+filt+".fits -WEIGHT_TYPE MAP_WEIGHT  -SATUR_LEVEL "+str(sat)+" -GAIN "+str(gain)+" -SEEING_FWHM "+str(seeing0)+" -STARNNW_NAME "+defa+"/default.nnw -CHECKIMAGE_NAME "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg,temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".subs  -CHECKIMAGE_TYPE SEGMENTATION,-BACKGROUND -PHOT_APERTURES 4,6,9,20,30"
                    os.system(cmd1)


                    hdufits_cat = fits.open("temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits")
                    cata = hdufits_cat[1].data
                    err_1= cata['FLUXERR_APER'][:,0]
                    err_2= cata['FLUXERR_APER'][:,1]
                    err_3= cata['FLUXERR_APER'][:,2]
                    err_4= cata['FLUXERR_APER'][:,3]
                    err_5= cata['FLUXERR_APER'][:,4]
                    flux_1= cata['FLUX_APER'][:,0]
                    flux_2= cata['FLUX_APER'][:,1]
                    flux_3= cata['FLUX_APER'][:,2]
                    flux_4= cata['FLUX_APER'][:,3]
                    flux_5= cata['FLUX_APER'][:,4]

                    flag= cata['FLAGS']

                    #sigma1,alpha,beta=get_error_model("temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".fits", "temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg",1,7,6)
                    sigma1,alpha,beta=get_error_model("temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".subs", "temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg",2,10,9)

                    print sigma1,alpha,beta


                    npix1=((4*0.5)**2)*np.pi

                    sigma_phot_1=np.sqrt(((sigma1**2)*(alpha**2)*(npix1**beta))+(flux_1/gain))

                    npix2=((6*0.5)**2)*np.pi

                    sigma_phot_2=np.sqrt(((sigma1**2)*(alpha**2)*(npix2**beta))+(flux_2/gain))

                    #print np.amin(np.sqrt((((alpha**2)*(npix2**beta))+(flux_2/((sigma1**2)*gain)))/(npix2+(flux_2/((sigma1**2)*gain)))))

                    #npix3=get_npix(int(9*0.5))
                    #print npix3
                    npix3=((9*0.5)**2)*np.pi
                    #print npix3
                    #sigma_phot_3=np.sqrt((err_3**2)*((((alpha**2)*(npix3**beta))+(flux_3/((sigma1**2)*gain)))/(npix3+(flux_3/((sigma1**2)*gain)))))
                    sigma_phot_3=np.sqrt(((sigma1**2)*(alpha**2)*(npix3**beta))+(flux_3/gain))


                    #npix4=get_npix(int(20*0.5))
                    #print npix4
                    npix4=((20*0.5)**2)*np.pi
                    #print npix4
                    #sigma_phot_4=np.sqrt((err_4**2)*((((alpha**2)*(npix4**beta))+(flux_4/((sigma1**2)*gain)))/(npix4+(flux_4/((sigma1**2)*gain)))))
                    sigma_phot_4=np.sqrt(((sigma1**2)*(alpha**2)*(npix4**beta))+(flux_4/gain))





                    orig_cols=cata.columns

                    c1=fits.Column(name='FLUXERR_OP',format='E',array=sigma_phot_1)

                    c2=fits.Column(name='FLUXERR_2',format='E',array=sigma_phot_2)


                    c3=fits.Column(name='FLUXERR_3',format='E',array=sigma_phot_3)

                    c4=fits.Column(name='FLUXERR_7',format='E',array=sigma_phot_4)

                    c5=fits.Column(name='FLUX_2',format='E',array=flux_2)



                    new_cols=fits.ColDefs([c1,c2,c3,c4,c5])


                    #thdu=fits.new_table(orig_cols+new_cols)
                    thdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
                    thdu.writeto("temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat2.fits")
                    #print np.amin(np.sqrt((((alpha**2)*(npix3**beta))+(flux_3/((sigma1**2)*gain)))/(npix3+(flux_3/((sigma1**2)*gain)))))

                    #plt.plot(err_2,sigma_phot_2,'ro')
                    #plt.show()

                    #f2=np.sqrt(std2**2/(np.mean(err_2**2)-(np.mean(flux_2)/gain)))
                    #f3=np.sqrt(std3**2/(np.mean(err_3**2)-(np.mean(flux_3)/gain)))
                    #f1=std1/np.median(err_1)
                    #f2=std2/np.median(err_2)
                    #f3=std3/np.median(err_3)


                    #print std2, np.mean(err_2)
                    #print std3, np.mean(err_3)
                    #print 1.0/f1,1.0/f2,1.0/f3

                    arch2.writelines("%s \t %d \t %f \t %f \t %f \t %f \t %f \t %f\n" % (list_img_name[i]+"_"+str(j),j,jd,zp,errzp,tim,ext,air))

                    cmd2="mv "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat2.fits "+fold_catalogs+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits"
                    os.system(cmd2)
                    #cmd3="rm "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".fits "+"temp/"+list_img_name[i]+"_conf_"+str(j)+"_"+filt+".fits "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits"
                    cmd3="rm "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".fits "+"temp/"+list_img_name[i]+"_conf_"+str(j)+"_"+filt+".fits "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".seg "+"temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".cat.fits"+" temp/"+list_img_name[i]+"_"+str(j)+"_"+filt+".subs "
                    os.system(cmd3)

                    print "done", list_img_name[i]+"_"+str(j)+"_"+filt+".fits"
                except:
                    arch3.writelines("%s \n " % (list_img_name[i]+"_"+str(j)+"_"+filt+".fits"))
        except:
            arch3.writelines("%s \n " % (list_img_name[i]+"_"+str(j)+"_"+filt+".fits"))
    arch.close()

arch2.close()


#elapsed = (time.clock() - start)
#print elapsed
