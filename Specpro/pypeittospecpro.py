# Made by Alden Beck
# Last edited Summer 2023
# This program takes DEIMOS Coadded 1d and 2d files and cuts them up per slit.

import numpy as np
from astropy.io import fits
import pandas as pd
import sys


#These packages are used to move files around nicely
import pathlib
import shutil
import os

# Do you want to make the image files as well? Default is false
write_images = False


# This is the number of mosaics your files have, set to 4 when running pypeit_run over all of the DIEMOS detector
# If the file only has one detector put 1
num_msc = 4

# Set up to take file names from command line, if none are given it will default to the ones given under else:
inputlist = sys.argv


if len(inputlist) > 1:
    name1d = inputlist[1]
    name2d = inputlist[2]
    textfile = inputlist[3]
    target = inputlist[4]
else:
    name2d = 'spec2d_DE.20140701.23439-DE.20140701.25960-a2261anew.fits'
    name1d = 'spec1d_DE.20140701.23439-DE.20140701.25960-a2261anew.fits'
    textfile = 'spec1d_DE.20140701.23439-DE.20140701.25960-a2261anew.txt'
    target = 'A2261'

# This is the mask that specpro used, idk what it means so I am just using what they have
mask = 'm46'

mscrng = num_msc + 1

hdulist1d = fits.open(name1d)
hdulist2d = fits.open(name2d)

#Creating paths and folders for files
currentpath = str(pathlib.Path(__file__).parent.resolve())
specpath = currentpath + f'/specpro'
imgpath = currentpath + f'/images'


if write_images:
    os.makedirs('images', exist_ok= True)
os.makedirs('specpro', exist_ok= True)



# ---------2D FUNCTIONS----------------
def findslitbound(hdulist, detnum):
    slitsindex = hdulist.index_of(f'MSC0{detnum}-SLITS')

    slits = hdulist[slitsindex].data
    numslits = len(slits)

    rightbound = []
    leftbound = []

    for i in range(numslits):
        currentslit = slits[i]
        leftbound.append(currentslit['left_init'])
        rightbound.append(currentslit['right_init'])
    left = np.reshape(leftbound, np.shape(leftbound))
    right = np.reshape(rightbound, np.shape(rightbound))

    return left, right


def findslitwave(hdulist, detnum):
    left, right = findslitbound(hdulist, detnum)
    fullwave = hdulist[hdulist.index_of(f'MSC0{detnum}-WAVEIMG')].data
    waves = []
    for i in range(len(left)):
        currentwave = fullwave[:, int(left[i, 0]):int(right[i, 0])]
        waves.append(currentwave)
    return waves


def findslitflux(hdulist, detnum):
    left, right = findslitbound(hdulist, detnum)
    fullimg = hdulist[hdulist.index_of(f'MSC0{detnum}-SCIIMG')].data
    fluxs = []
    for i in range(len(left)):
        currentimg = fullimg[:, int(left[i, 0]):int(right[i, 0])]
        fluxs.append(currentimg)
    return fluxs


def findslitivarraw(hdulist, detnum):
    left, right = findslitbound(hdulist, detnum)
    fullivar = hdulist[hdulist.index_of(f'MSC0{detnum}-IVARRAW')].data
    ivars = []
    for i in range(len(left)):
        currentivar = fullivar[:, int(left[i, 0]):int(right[i, 0])]
        ivars.append(currentivar)
    return ivars


def findslitivarmodel(hdulist, detnum):
    left, right = findslitbound(hdulist, detnum)
    fullivar = hdulist[hdulist.index_of(f'MSC0{detnum}-IVARMODEL')].data
    ivars = []
    for i in range(len(left)):
        currentivar = fullivar[:, int(left[i, 0]):int(right[i, 0])]
        ivars.append(currentivar)
    return ivars


def all2dmsc(hdulist):
    leftall = []
    rightall = []
    fluxsall = []
    ivarsall = []
    wavesall = []
    namesall = []

    for detnum in range(1, mscrng):
        slitsindex = hdulist.index_of(f'MSC0{detnum}-SLITS')
        slits = hdulist[slitsindex].data
        spat_list = list(slits['spat_id'])
        for l in range(len(spat_list)):
            namesall.append([spat_list[l], detnum])

    for i in range(1, mscrng):
        leftloop, rightloop = findslitbound(hdulist, i)
        wavesloop = findslitwave(hdulist, i)
        fluxsloop = findslitflux(hdulist, i)
        ivarsloop = findslitivarraw(hdulist, i)

        leftall.extend(leftloop)
        rightall.extend(rightloop)
        wavesall.extend(wavesloop)
        fluxsall.extend(fluxsloop)
        ivarsall.extend(ivarsloop)

    return leftall, rightall, wavesall, fluxsall, ivarsall, namesall


# ----------------1D FUNCTIONS----------------------
def check1das2d(detnum, name, hdu2d=hdulist2d):
    #This function is probably useless but I am too lazy to get rid of it
    slitsindex = hdu2d.index_of(f'MSC0{detnum}-SLITS')
    slits = hdu2d[slitsindex].data

    spat_list = list(slits['spat_id'])
    split1 = name.split('SLIT')[1]
    number = int(split1.split('-MSC0')[0])
    if number in spat_list:
        return True
    else:
        return False


def extract1d(hdulist, detnum):
    indexs = []
    names = []
    for i in range(1, len(hdulist)):
        if ('SPAT' in hdulist[i].header['NAME']) and (f'MSC0{detnum}' in hdulist[i].header['NAME']):
            if check1das2d(detnum, hdulist[i].header['NAME']):
                indexs.append(i)

                name = hdulist[i].header['NAME']
                split1 = name.split('SLIT')[1]
                number = int(split1.split('-MSC0')[0])
                names.append([number, detnum])
    countsall = []
    ivarall = []
    waveall = []
    mask = []


    for i in indexs:
        data = hdulist[i].data
        #print(names[i-1])
        #print(indexs[i-1])
        #print(len(names), len(indexs))
        #print(names[i-1])
        extract_method_string = 'OPT'
        try:
            data['OPT_COUNTS']
            data['OPT_COUNTS_IVAR']

        except KeyError:
            extract_method_string = 'BOX'
            print(f'slit-{i} in MSC0{detnum} OPT extraction not found!!! Using BOX extraction as backup')



        countstring = extract_method_string + '_COUNTS'
        ivarstring = extract_method_string + '_COUNTS_IVAR'
        wavestring = extract_method_string + '_WAVE'
        maskstring = extract_method_string + '_MASK'

        countsall.append(data[countstring])
        ivarall.append(data[ivarstring])
        waveall.append(data[wavestring])
        mask.append(data[maskstring])

    return countsall, ivarall, waveall, names, mask


def all1dmsc(hdulist):
    countsall = []
    ivarall = []
    waveall = []
    namesall = []
    maskall = []
    for detnum in range(1, mscrng):
        currentcount, currentivar, currentwave, names, mask = extract1d(hdulist, detnum)
        countsall.extend(currentcount)
        ivarall.extend(currentivar)
        waveall.extend(currentwave)
        namesall.extend(names)
        maskall.extend(mask)

    return countsall, ivarall, waveall, namesall, maskall


# ---------------INFO FILE---------------------

def find2dradec(hdu2d):  # SLIT RA AND DEC
    ral = []
    decl = []
    for detnum in range(1, mscrng):
        maskindex = hdu2d.index_of(f'MSC0{detnum}-MASKDEF_DESIGNTAB')
        data = hdu2d[maskindex].data
        for i in range(len(data)):
            currentslit = data[i]
            ral.append(currentslit['SLITRA'])
            decl.append(currentslit['SLITDEC'])
    ra = np.reshape(ral, np.shape(ral))
    dec = np.reshape(decl, np.shape(decl))
    return ra, dec


def find1dradec(txtfile):  # RA AND DEC
    file = pd.read_csv(txtfile, sep='|')
    for i in range(len(list(file.keys()))):
        if 'objra' in list(file.keys())[i]:
            raindex = list(file.keys())[i]
        if 'objdec' in list(file.keys())[i]:
            decindex = list(file.keys())[i]

    ra = np.array(file[raindex].to_list())
    dec = np.array(file[decindex].to_list())

    return ra, dec


def findextractpos(txtfile):
    file = pd.read_csv(txtfile, sep='|')
    for i in range(len(list(file.keys()))):
        if 'spat_pixpos' in list(file.keys())[i]:
            spatindex = list(file.keys())[i]
    spat_pixos = np.array(file[spatindex].to_list())
    return spat_pixos


#This function is probably useless now
def findboxwidth(txtfile):
    file = pd.read_csv(txtfile, sep='|')
    for i in range(len(list(file.keys()))):
        if 'box_width' in list(file.keys())[i]:
            boxindex = list(file.keys())[i]
    box_width = np.array(file[boxindex].to_list())
    return box_width


def findboxnpix(hdu1d):
    indexs = []
    names = []
    for detnum in range(1, mscrng):
        for i in range(1, len(hdu1d)):
            if ('SPAT' in hdu1d[i].header['NAME']) and (f'MSC0{detnum}' in hdu1d[i].header['NAME']):
                if check1das2d(detnum, hdu1d[i].header['NAME']):
                    indexs.append(i)

                    name = hdu1d[i].header['NAME']
                    split1 = name.split('SLIT')[1]
                    number = int(split1.split('-MSC0')[0])
                    names.append([number, detnum])
    
    boxnpix = []
    for i in indexs:
        data = hdu1d[i].data
        boxnpix.append(data['BOX_NPIX'])
    return boxnpix


def findslitnums(hdu2d):
    slitlenl = []
    slitwidl = []
    slitpal = []

    for detnum in range(1, mscrng):
        maskindex = hdu2d.index_of(f'MSC0{detnum}-MASKDEF_DESIGNTAB')
        data = hdu2d[maskindex].data
        for i in range(len(data)):
            currentslit = data[i]

            slitlenl.append(currentslit['SLITLEN'])
            slitwidl.append(currentslit['SLITWID'])
            slitpal.append(currentslit['SLITPA'])
    slitlen = np.reshape(slitlenl, np.shape(slitlenl))
    slitwid = np.reshape(slitwidl, np.shape(slitwidl))
    slitpa = np.reshape(slitpal, np.shape(slitpal))
    return slitlen, slitwid, slitpa



# --------------WRITING SPEC DATA------------------

def findindex(namesall1d, namesall2d):
    #leftall2d, rightall2d, wavesall2d, fluxsall2d, ivarsall2d, namesall2d = all2dmsc(hdu2d)
    #countsall1d, ivarsall1d, wavesall1d, namesall1d, mask1d = all1dmsc(hdu1d)

    namesindex = []
    indexs = []
    for i in range(len(namesall2d)):
        name2 = namesall2d[i]
        index2 = i
        index1 = namesall1d.index(name2)
        name1 = namesall1d[index1]
        indexs.append([index1, index2])
        namesindex.append([name1, name2])
    test = []
    for i in range(len(namesall1d)):
        name1 = namesall1d[i]
        index1 = i
        index2 = namesall2d.index(name1)
        name2 = namesall2d[index2]
        #namesindex.append([name1, name2])
        test.append([index1,index2])
    return test, namesindex


def removebadmask(wave1,wave2,ivar1,ivar2,flux1,flux2,mask,boxnpix):
    print(f'Cutting arrays down to {len(wave1[mask])} from {len(mask)} rows')
    assert False not in mask[mask]
    return wave1[mask], wave2[mask,:], ivar1[mask], ivar2[mask,:], flux1[mask], flux2[mask,:], boxnpix[mask]


def compilefiles(hdu1d, hdu2d, textfile):
    leftall2d, rightall2d, wavesall2d, fluxsall2d, ivarsall2d, namesall2d = all2dmsc(hdu2d)
    countsall1d, ivarsall1d, wavesall1d, namesall1d, mask1d = all1dmsc(hdu1d)

    indexs, namesindex = findindex(namesall1d, namesall2d)

    print(f'Found {indexs[-1][0] + 1} 1d spectra and {indexs[-1][1] + 1} 2d spectra')
    print(f'Creating files for {indexs[-1][0] + 1} slits \n')


    for i in range(len(indexs)):
        #---------Info file data-------------------------
        ra1d, dec1d = find1dradec(textfile)
        boxnpix = findboxnpix(hdu1d)
        spat_pos1d = findextractpos(textfile)
        ra2d, dec2d = find2dradec(hdu2d)

        print(ra1d[i], dec1d[i], ra2d[i], dec2d[i])
        exit()
        #---------------FITS data----------------------
        i1d, i2d = indexs[i]

        primary2dhdu = fits.PrimaryHDU()
        slitnum2d = str(np.char.zfill(str(i + 1), 3))
        filename2d = f'spec2d.{mask}.{slitnum2d}.{target}.fits'

        primary1dhdu = fits.PrimaryHDU()
        filename1d = f'spec1d.{mask}.{slitnum2d}.{target}.fits'

        print(f'({i1d+1}/{indexs[-1][0] + 1}), 1d index={i1d}, 2d index={i2d}')
        
        #---------------------------Remove the 0's---------------------------
        cwave1, cwave2, civar1, civar2, cflux1, cflux2, box_width1d= removebadmask(
            wavesall1d[i1d],
            wavesall2d[i2d],
            ivarsall1d[i1d],
            ivarsall2d[i2d],
            countsall1d[i1d],
            fluxsall2d[i2d],
            mask1d[i1d],
            boxnpix[i1d])
        

        # -----------Data Structuring-----------



        try:
            lenlong = len(np.array(cwave2)[:, 1])
            lenshort = len(np.array(cwave2)[1])
        except:
            print(f'Size of array is 0. Skipping slit {i1d + 1}')
            continue

        a1 = np.ones((1, lenshort, lenlong))
        a2 = np.ones((1, lenshort, lenlong))
        a3 = np.ones((1, lenshort, lenlong))

        #a1[0] = np.array(np.array(cflux2).transpose())
        #a2[0] = np.array(np.array(civar2).transpose())
        #a3[0] = np.array(np.array(cwave2).transpose())

        a1[0] = np.array(cflux2.transpose())
        a2[0] = np.array(civar2.transpose())
        a3[0] = np.array(cwave2.transpose())

        length1d = len(civar1)

        b1 = np.ones((1, length1d))
        b2 = np.ones((1, length1d))
        b3 = np.ones((1, length1d))

        #b1[0] = np.array(cflux1)
        #b2[0] = np.array(civar1)
        #b3[0] = np.array(cwave1)

        b1[0] = cflux1
        b2[0] = civar1
        b3[0] = cwave1

        # -------FILE PACKAGING----------------
        dimstring1d = f'{length1d}'
        dimstring2d = f'({lenlong},{lenshort})'

        fluxcolumn2d = fits.Column(name='FLUX', array=a1, format=str(np.size(a1)) + 'D', dim=dimstring2d)
        ivarcolumn2d = fits.Column(name='IVAR', array=a2, format=str(np.size(a2)) + 'D', dim=dimstring2d)
        wavecolumn2d = fits.Column(name='LAMBDA', array=a3, format=str(np.size(a3)) + 'D', dim=dimstring2d)

        tablehdu2d = fits.BinTableHDU.from_columns([fluxcolumn2d, ivarcolumn2d, wavecolumn2d])

        new2dhdu = fits.HDUList([primary2dhdu, tablehdu2d])

        fluxcolumn1d = fits.Column(name='FLUX', array=b1, format=str(np.size(b1)) + 'D', dim=dimstring1d)
        ivarcolumn1d = fits.Column(name='IVAR', array=b2, format=str(np.size(b2)) + 'D', dim=dimstring1d)
        wavecolumn1d = fits.Column(name='LAMBDA', array=b3, format=str(np.size(b3)) + 'D', dim=dimstring1d)

        tablehdu1d = fits.BinTableHDU.from_columns([fluxcolumn1d, wavecolumn1d, ivarcolumn1d])
        new1dhdu = fits.HDUList([primary1dhdu, tablehdu1d])

        print(f'Writing {filename1d} and {filename2d}')

        # -----------INFO FILE-----------------
        #These are declared earlier so can be removed
        #ra1d, dec1d = find1dradec(textfile)
        #box_width1d = findboxwidth(textfile)
        #boxnpix = findboxnpix(hdu1d)
        #spat_pos1d = findextractpos(textfile)
        #ra2d, dec2d = find2dradec(hdu2d)

        objid = slitnum2d  # Maybe change this in the future to an ID unique to the slit from the txtfile
        slitlen2d, slitwid2d, slitpa2d = findslitnums(hdu2d)

        zphot = 0.0
        zpdf = 0.0
        zpdflow = 0.0
        zpdfhigh = 0.0


        txtnames = ['ID', 'RA', 'DEC', 'extractpos', 'extractwidth', 'slitRA', 'slitDEC',
                     'slitlen', 'slitwid', 'slitPA', 'zphot', 'zpdf', 'zpdf-low', 'zpdf-up']
        txtnumbers = [objid, ra1d[i1d], dec1d[i1d], spat_pos1d[i1d] - leftall2d[i2d][0], np.average(box_width1d), 
                      ra2d[i2d], dec2d[i2d], slitlen2d[i2d], slitwid2d[i2d], slitpa2d[i2d], zphot, zpdf, zpdflow,
                      zpdfhigh]
        

        txtfilename = f'info.{mask}.{slitnum2d}.{target}.dat'

        txtcol = np.column_stack((txtnames, txtnumbers))

        print(f'Writing {txtfilename}')
        np.savetxt(txtfilename, txtcol, delimiter=' ', fmt='%s')

        new2dhdu.writeto(filename2d, overwrite=True)
        new1dhdu.writeto(filename1d, overwrite=True)

        new1dhdu.close()
        new2dhdu.close()


        #-----Making the image file----------
        if write_images:
            imgfilename = f'{target}.{slitnum2d}.fits'

            imagePrimaryHDU = fits.PrimaryHDU(a1)

            imagePrimaryHDU.writeto(imgfilename, overwrite = True)
            shutil.move(currentpath + '/' + imgfilename, imgpath+'/'+imgfilename)
            print(f'Writing image file {imgfilename} \n')


        #---Moving Files to new folders---------
        shutil.move(currentpath + '/' + filename1d, specpath+'/'+filename1d)
        shutil.move(currentpath + '/' + filename2d, specpath+'/'+filename2d)
        shutil.move(currentpath + '/' + txtfilename, specpath+'/'+txtfilename)
        
       

    print('Done')

compilefiles(hdulist1d, hdulist2d, textfile)

hdulist2d.close()
hdulist1d.close()