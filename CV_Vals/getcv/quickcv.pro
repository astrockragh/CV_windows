function quickcv,side1,side2,za,deltaz,sig8=sig8,npoints=npoints
;+
; NAME:
;
;  QUICKCV
;
; PURPOSE:
;
;   calculate cosmic variance for some geometry
;
; CATEGORY:
;
;   cosmic variance
;
;
; CALLING SEQUENCE:
;
;   fracerror=quickcv(side1,side2 [zarr,deltaz,sig8=sig8,npoints=npoints]) 
;
; INPUTS:
;
;   side1: length of 1 side of field on sky, in degrees (rect. geom)
;   side2: length of other side of field on sky, in degrees
;
; OPTIONAL INPUTS:
;   zarr: array of z's for which cosmic variance will be calculated; default:1
;   deltaz: total length of volume in z direction; default: 0.1
;
; KEYWORD PARAMETERS:
;   sig8: desired value of sigma_8; default: 0.77
;   npoints: desired number of integration points for first part of 
;     integration; default: 200
;
;
; OUTPUTS:
;
;  fracerror: array containing the fractional error in a count,
;  sigma/N, due to cosmic variance for objects with bias b=1 
;  in a volume side1 x side2 degrees x delta_z centered at redshifts
;  zarr.  NOTE THIS IS SIGMA, NOT VARIANCE!!!!
;
; RESTRICTIONS:
;
;  uses power spectrum in pofk.pro (gamma=0.1872, omega0=0.26), distance
;  relation in rz.pro 
;
;
; MODIFICATION HISTORY:
; modified by rss july 2006 to include Dlin scaling with z
; modified by bpm sep  2009 to wmap3 cosmology
;-

if n_elements(npoints) eq 0 then npoints=200.

; desired sigma8 
if n_elements(sig8) eq 0 then sig8=0.81

if n_elements(deltaz) eq 0 then deltaz=0.10


if n_elements(za) eq 0 then za=1.
nz=n_elements(za)

nfields=1


; assume LCDM
; uses jonathan baker's old code

rza=rz(za,omegam=0.31,omegal=0.69)
rzmin=rz(za-deltaz/2,omegam=0.31,omegal=0.69)
rzmax=rz(za+deltaz/2,omegam=0.31,omegal=0.69)
cvi=fltarr(nz)


x1a=3000./2.*(rza * side1)/!radeg
x2a=3000./2.*(rza * side2)/!radeg
x3a=3000.*(rzmax-rzmin)/2.

for i=0,nz-1 do  $
  cvi(i)=intpk4(x1=x1a(i),x2=x2a(i),x3=x3a(i),nsteps=npoints)

;cvi is the fractional variance in a count in a redshift bin
;cvar is the fractional error (sigma)

cvi=cvi*sig8^2/1.033^2
  d = Dlin(0.26, za)
  ;print, 'dlin: ', d
  cvar=sqrt(cvi)*d

;print,'Fractional sigma: ',cvar

  return,cvar

end
