; NAME:
;
;  GETCV
;
; PURPOSE:
;
;   calculate cosmic variance for given geometry as function of redshift bin and stellar mass
;   driver program for for quickcv including galaxy bias
;
; CATEGORY:
;
;   cosmic variance for WMAP3 cosmology
;   Updated for Planck 2018
;
;
; INPUTS:
;
;   side1: length of 1 side of field on sky, in degrees (rect. geom)
;   side2: length of other side of field on sky, in degrees
;   zarr: array of z's for which cosmic variance will be calculated
;
;
; OUTPUTS:
;  
;  prints table containing the following columns:
;  1) mean redshift
;  2) redshift bin size
;  3) fractional error in a count for dark matter
;  4) fractional error in a count (sigma/N) for galaxies with 8.5 < log(m/M_sun) 9.0
;  5) fractional error in a count (sigma/N) for galaxies with 9.0 < log(m/M_sun) 9.5
;  6) fractional error in a count (sigma/N) for galaxies with 9.5 < log(m/M_sun) 10.0
;  7) fractional error in a count (sigma/N) for galaxies with 10.0 < log(m/M_sun) 10.5
;  8) fractional error in a count (sigma/N) for galaxies with 10.5 < log(m/M_sun) 11.0
;  9) fractional error in a count (sigma/N) for galaxies with 11.0 < log(m/M_sun) 11.5
;  in a volume side1 x side2 degrees x delta_z centered at redshifts
;  zarr.  NOTE THIS IS SIGMA, NOT VARIANCE!!!!
;
; RESTRICTIONS:
;
;  uses power spectrum in pofk.pro (gamma=0.1872, omega0=0.26), distance
;  relation in rz.pro, growth factor in dlin.pro
;
; MODIFICATION HISTORY:
;   released DEC 09
;   added more mass bins
;   changed cosmo parameters to match Planck






; -------- Input parameters ---------------------------------------------------------------------

 s1 = 10.0/60.0 ;side 1, in degrees 
 s2 = 20.0/60.0 ;side 2, in degrees
 zarr = [0, 0.04, 0.07, 0.11, 0.15, 0.19, 0.22, 0.26, 0.31, 0.35, 0.39, 0.43, 0.48, 0.53, 0.57, 0.62, 0.67, 0.73, 0.78, 0.84, 0.89, 0.95, 1.01, 1.08, 1.14, 1.21, 1.28, 1.36, 1.43, 1.51, 1.6, 1.69, 1.78, 1.87, 1.98, 2.08, 2.19, 2.31, 2.43, 2.56, 2.7, 2.84, 2.99, 3.15, 3.32, 3.51, 3.7, 3.9, 4.12, 4.35, 4.6, 4.86, 5.15, 5.45, 5.78, 6.14, 6.52, 6.93, 7.38, 7.87, 8.4, 8.98, 9.61, 10.3, 11.07, 11.92, 12.86, 13.91, 15.07, 16.39, 17.87, 19.54, 21.46]
; -----------------------------------------------------------------------------------------------


nz = n_elements(zarr)-1

; The bias values for 7.0, 7.5, and 8.0 are incorrectly extrapolated and should be removed

b070=0.072
b170=2.50
b270=0.80

b075=0.067
b175=2.62
b275=0.87

b080=0.062
b180=2.74
b280=0.94

b085=0.062
b185=2.59
b285=1.025

b090=0.074
b190=2.58
b290=1.039

b095=0.042
b195=3.17
b295=1.147

b0100=0.053
b1100=3.07
b2100=1.225

b0105=0.069
b1105=3.19
b2105=1.269

b0110=0.173
b1110=2.89
b2110=1.438


  for i=0, nz-1 do begin
      dz      = zarr[i+1]-zarr[i]
      zmid    = zarr[i]+0.5*dz
      s       = quickcv(s1, s2, zmid, dz, sig8=0.81)
      bias70  = b070  * (zmid+1)^b170  + b270
      bias75  = b075  * (zmid+1)^b175  + b275
      bias80  = b080  * (zmid+1)^b180  + b280
      bias85  = b085  * (zmid+1)^b185  + b285
      bias90  = b090  * (zmid+1)^b190  + b290
      bias95  = b095  * (zmid+1)^b195  + b295
      bias100 = b0100 * (zmid+1)^b1100 + b2100
      bias105 = b0105 * (zmid+1)^b1105 + b2105
      bias110 = b0110 * (zmid+1)^b1110 + b2110
      sgg70  = bias70  * s
      sgg75  = bias75  * s
      sgg80  = bias80  * s
      sgg85  = bias85  * s
      sgg90  = bias90  * s
      sgg95  = bias95  * s
      sgg100 = bias100 * s
      sgg105 = bias105 * s
      sgg110 = bias110 * s
      print, zmid, dz, s, sgg70, sgg75, sgg80, sgg85, sgg90, sgg95, sgg100, sgg105, sgg110, FORMAT = '(F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4," ",F7.4)'

;prints the mean redshift, redshift bins size, sigma_dm and the sigma_gg for nine stellar mass bins

  endfor

end



