
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval
from astroquery.vizier import Vizier
import subprocess
import os
def process_image(image_filename):

    maxmag = 18
    boxsize = 30
    catNum = 'II/349'

    with fits.open(image_filename) as HDUList:
        image_hdu = next((hdu for hdu in HDUList if hdu.data is not None), None)
        if image_hdu is None:
            raise ValueError(f"No image data found in {image_filename}")
        header = image_hdu.header
        image = image_hdu.data
    w = WCS(header)
    naxis = w.naxis
    if naxis == 2:
        (raImage, decImage) = w.all_pix2world(image.shape[1]/2, image.shape[0]/2, 0)
    elif naxis == 3:
        (raImage, decImage, thirdAxis) = w.all_pix2world(image.shape[2]/2, image.shape[1]/2, image.shape[0]/2, 0)
    else:
        raise ValueError(f"Unexpected number of WCS axes: {naxis}")
    v = Vizier(columns=['*'], column_filters={"gmag":f"<{maxmag}", "Nd":">6", "e_gmag":f"<{1.086/3}"}, row_limit=-1)
    Q = v.query_region(SkyCoord(ra=raImage, dec=decImage, unit=(u.deg, u.deg)), radius=str(boxsize)+'m', catalog=catNum, cache=False)
    good_cat_stars = Q[0]

    catalogName = image_filename + '.cat'
    subprocess.run(['source-extractor', '-c', 'photomCat.sex', image_filename, '-CATALOG_NAME', catalogName, '-PARAMETERS_NAME', 'photomCat.param'], check=True)
    subprocess.run(['psfex', '-c', 'psfex_conf.psfex', catalogName], check=True)

    psfName = image_filename + '.psf'
    psfcatalogName = image_filename + '.psf.cat'
    subprocess.run(['source-extractor', '-c', 'photomCat.sex', image_filename, '-CATALOG_NAME', psfcatalogName, '-PSF_NAME', psfName, '-PARAMETERS_NAME', 'photomPSF.param'], check=True)

    with fits.open(psfcatalogName) as HDU:
        psfsourceTable = Table(HDU[2].data)
    cleanPSFSources = psfsourceTable[(psfsourceTable['FLAGS']==0) & (psfsourceTable['FLAGS_MODEL']==0)  & (psfsourceTable['FWHM_WORLD'] < 2) & (psfsourceTable['XMODEL_IMAGE']<3500) & (psfsourceTable['XMODEL_IMAGE']>500) &(psfsourceTable['YMODEL_IMAGE']<3500) &(psfsourceTable['YMODEL_IMAGE']>500)]

    psfsourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHAWIN_J2000'], dec=cleanPSFSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
    ps1CatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
    photoDistThresh = 0.6
    idx_psfimage, idx_psfps1, d2d, d3d = ps1CatCoords.search_around_sky(psfsourceCatCoords, photoDistThresh*u.arcsec)
    psfoffsets = ma.array(good_cat_stars['gmag'][idx_psfps1] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)

    ra = 210.910674637
    dec = 54.3116510708
    sn2023ixf_coords = SkyCoord(ra=[ra], dec=[dec], frame='icrs', unit='degree')
    idx_sn2023ixf, idx_cleanpsf_sn2023ixf, d2d, d3d = psfsourceCatCoords.search_around_sky(sn2023ixf_coords, photoDistThresh*u.arcsec)
    sn2023ixf_psfinstmag = cleanPSFSources[idx_cleanpsf_sn2023ixf]['MAG_POINTSOURCE'][0]
    sn2023ixf_psfinstmagerr = cleanPSFSources[idx_cleanpsf_sn2023ixf]['MAGERR_POINTSOURCE'][0]
    sn2023ixf_psfmag = zero_psfmed + sn2023ixf_psfinstmag
    sn2023ixf_psfmagerr = np.sqrt(sn2023ixf_psfinstmagerr**2 + zero_psfstd**2)

    return sn2023ixf_psfmag, sn2023ixf_psfmagerr
image_files = [f for f in os.listdir('.') if f.endswith('.fits')]

results = []

for image_file in image_files:
    if image_file.endswith('.fits') and  image_file.startswith('2023'): 
        print("Processing file:", image_file)
        mag, mag_err = process_image(image_file)
        results.append((image_file, mag, mag_err))

for result in results:
    print(f"Image: {result[0]}, Magnitude: {result[1]:.2f} +/- {result[2]:.2f}")

{output}

#Final Result

Image: 20230605170711-835-RA.wcs.proc.fits, Magnitude: 11.30 +/- 0.20
Image: 20230525184701-216-RA.wcs.proc.fits, Magnitude: 11.03 +/- 0.22
Image: 20230521154933-241-RA.wcs.proc.fits, Magnitude: 11.44 +/- 0.22
Image: 20230528152042-139-RA.wcs.proc.fits, Magnitude: 11.06 +/- 0.21
Image: 20230605152158-307-RA.wcs.proc.fits, Magnitude: 11.29 +/- 0.22
Image: 20230526154421-699-RA.wcs.proc.fits, Magnitude: 11.04 +/- 0.21
Image: 20230525195322-030-RA.wcs.proc.fits, Magnitude: 10.94 +/- 0.18
Image: 20230527144518-321-RA.wcs.proc.fits, Magnitude: 11.00 +/- 0.20
Image: 20230530180249-785-RA.wcs.proc.fits, Magnitude: 11.12 +/- 0.19
Image: 20230522153045-729-RA.wcs.proc.fits, Magnitude: 11.04 +/- 0.19
Image: 20230606145107-018-RA.wcs.proc.fits, Magnitude: 11.32 +/- 0.22
Image: 20230524182725-957-RA.wcs.proc.fits, Magnitude: 10.80 +/- 0.17
Image: 20230606151355-919-RA.wcs.proc.fits, Magnitude: 11.32 +/- 0.22
Image: 20230521151417-752-RA.wcs.proc.fits, Magnitude: -- +/- 0.00
Image: 20230520151900-863-RA.wcs.proc.fits, Magnitude: 12.51 +/- 0.21
Image: 20230526152122-558-RA.wcs.proc.fits, Magnitude: 11.01 +/- 0.21
Image: 20230524180955-720-RA.wcs.proc.fits, Magnitude: 10.80 +/- 0.17
Image: 20230604151623-740-RA.wcs.proc.fits, Magnitude: 11.25 +/- 0.22
Image: 20230528153428-317-RA.wcs.proc.fits, Magnitude: 11.02 +/- 0.20
Image: 20230605160040-731-RA.wcs.proc.fits, Magnitude: 11.29 +/- 0.20
Image: 20230605163342-989-RA.wcs.proc.fits, Magnitude: 11.37 +/- 0.20
Image: 20230525190514-248-RA.wcs.proc.fits, Magnitude: 11.00 +/- 0.20
Image: 20230603151714-483-RA.wcs.proc.fits, Magnitude: 11.24 +/- 0.21
Image: 20230527150837-913-RA.wcs.proc.fits, Magnitude: 11.06 +/- 0.21
Image: 20230605172106-597-RA.wcs.proc.fits, Magnitude: 11.29 +/- 0.21
Image: 20230522154624-727-RA.wcs.proc.fits, Magnitude: 11.15 +/- 0.23
Image: 20230603150337-154-RA.wcs.proc.fits, Magnitude: 11.23 +/- 0.20
Image: 20230604153026-612-RA.wcs.proc.fits, Magnitude: 11.25 +/- 0.21
Image: 20230523210842-911-RA.wcs.proc.fits, Magnitude: 10.92 +/- 0.20\
