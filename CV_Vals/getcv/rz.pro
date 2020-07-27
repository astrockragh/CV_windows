
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 
;; Routines for computing the Alcock-Paczynski distortion
;;
;; 1998/10  Jonathan Baker <jbaker@astro.berkeley.edu>
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; H(z) in units where H(0) = 1

FUNCTION hubble, z

   COMMON cosmol, OmegaM, OmegaL, OmegaQ, wQ

   x = 1. + z
   OmegaK = 1. - OmegaM - OmegaL - OmegaQ
   h2 = OmegaM * x^3 + OmegaK * x^2 + OmegaL + OmegaQ * x^(3.+3.*wQ)
   
   RETURN, sqrt(h2)

END


; Compute d(chi)/dz where d(chi)^2 = dr^2 / (1-kr^2)

FUNCTION dchi_dz, z

   RETURN, 1./hubble(z)

END


; Compute coordinate distance r(z) for FRW metric

FUNCTION rz, z, $
            OmegaMatter = OmegaM, $
            OmegaLambda = OmegaL, $
            OmegaQ = OmegaQ, $
            wQ = wQ, $
            eps = eps

   COMMON cosmol, cOmegaM, cOmegaL, cOmegaQ, cwQ
   
   IF NOT keyword_set(OmegaM) THEN OmegaM = 1.d0
   IF NOT keyword_set(OmegaL) THEN OmegaL = 0.d0
   IF NOT keyword_set(OmegaQ) THEN OmegaQ = 0.d0
   IF NOT keyword_set(wQ) THEN wQ = 0.d0
   
   cOmegaM = OmegaM
   cOmegaL = OmegaL
   cOmegaQ = OmegaQ
   cwQ = wQ
   
   kurv = OmegaM + OmegaL + OmegaQ - 1.d0
   
   nz = n_elements(z)
   dchi = dblarr(nz)
   chi = dblarr(nz)
   r = dblarr(nz)
   
   z2 = z[0]
   IF z2 EQ 0. THEN $
     dchi[0] = 0.d0 $
   ELSE $
     dchi[0] = qromb('dchi_dz', 0., z2, /double, eps=eps)
   
   FOR i = 1, nz-1 do begin
      z1 = z[i-1]
      z2 = z[i]
      dchi[i] = qromb('dchi_dz', z1, z2, /double, eps=eps)
   ENDFOR
   
   chi[0] = dchi[0]
   FOR i = 1, nz-1 do $
     chi[i] = chi[i-1] + dchi[i]
   
   IF abs(kurv) LT 1.e-4 THEN $ ; flat
     r = chi $
   ELSE IF kurv GT 0.d0 THEN $  ; closed
     r = sin(chi*sqrt(kurv))/sqrt(kurv) $
   ELSE $                       ; open
     r = sinh(chi*sqrt(-kurv))/sqrt(-kurv)
   
   RETURN, r
   
END


; Make Fig. 13.5 from Peebles (1993) book: 
; angular size / physical size vs. z

PRO alcock_test1, $
                  zMin = zMin, $
                  zMax = zMax, $
                  nz = nz
   
   IF NOT keyword_set(zMin) THEN zMin = 0.001d0
   IF NOT keyword_set(zMax) THEN zMax = 10.d0
   IF NOT keyword_set(nz) THEN nz = 200
   
   dz = (zMax-zMin) / (nz-1.d0)
   z = zMin + dz*findgen(nz)
   
   lincolr
   
; models with no curvature
   
   window, 0
   r = rz(z) & plot, z, (1+z)/r, yr=[0,10], /ys, color=1
   r = rz(z, OmegaM=2, OmegaL=-1) & oplot, z, (1+z)/r, color=2
   r = rz(z, OmegaM=0.5, OmegaL=0.5) & oplot, z, (1+z)/r, color=3
   r = rz(z, OmegaM=0.2, OmegaL=0.8) & oplot, z, (1+z)/r, color=4
   r = rz(z, OmegaM=0.1, OmegaL=0.9) & oplot, z, (1+z)/r, color=5
   r = rz(z, OmegaM=0.05, OmegaL=0.95) & oplot, z, (1+z)/r, color=6
   
; models with no lambda
   
   window, 1
   r = rz(z) & plot, z, (1+z)/r, yr=[0,10], /ys, color=1
   r = rz(z, OmegaM=2) & oplot, z, (1+z)/r, color=2
   r = rz(z, OmegaM=0.5) & oplot, z, (1+z)/r, color=3
   r = rz(z, OmegaM=0.2) & oplot, z, (1+z)/r, color=4
   r = rz(z, OmegaM=0.1) & oplot, z, (1+z)/r, color=5
   r = rz(z, OmegaM=0.05) & oplot, z, (1+z)/r, color=6
   
END


; Plot comoving size / angular size vs. z

PRO alcock_plot1, $
                  zMin = zMin, $
                  zMax = zMax, $
                  nz = nz, $
                  _extra = extra
   
   IF NOT keyword_set(zMin) THEN zMin = 0.d0
   IF NOT keyword_set(zMax) THEN zMax = 4.d0
   IF NOT keyword_set(nz) THEN nz = 200
   
   dz = (zMax-zMin) / (nz-1.d0)
   z = zMin + dz*findgen(nz)
   
   lincolr
   !p.multi = [0,1,2]
   
   r = rz(z) 
   mplot, z, r, _extra=extra, $
     mxtitle="z", mytitle=textoidl("H_0 r(z) / c"), $
     mtitle="Comoving size / Angular size", $
     charsize=1.5
   
   r = rz(z, OmegaM=0.2, OmegaL=0.8) & oplot, z, r, linesty=4, color=3
   r = rz(z, OmegaM=0.4, OmegaL=0.6) & oplot, z, r, linesty=5, color=3
   r = rz(z, OmegaM=0.2) & oplot, z, r, linesty=2, color=2
   r = rz(z, OmegaM=0.4) & oplot, z, r, linesty=3, color=2
   oplot, [0,10], [0,10], linesty=1
   
   xr = [0.55, 0.65]
   yr = [0.4, 0.32, 0.24, 0.16]
   reg2dev, xr, yr[0]*[1,1], x, y
   plots, x, y, linesty=4, color=3, /device
   reg2dev, xr, yr[1]*[1,1], x, y
   plots, x, y, linesty=5, color=3, /device
   reg2dev, xr, yr[2]*[1,1], x, y
   plots, x, y, linesty=2, color=2, /device
   reg2dev, xr, yr[3]*[1,1], x, y
   plots, x, y, linesty=3, color=2, /device

   xr = [0.67]
   yr = yr - 0.02
   reg2dev, xr, [yr[0]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.2, \Omega_\Lambda=0.8"), /device, $
     color=3, charsize=1.4
   reg2dev, xr, [yr[1]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.4, \Omega_\Lambda=0.6"), /device, $
     color=3, charsize=1.4
   reg2dev, xr, [yr[2]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.2, \Omega_\Lambda=0"), /device, $
     color=2, charsize=1.4
   reg2dev, xr, [yr[3]], x, y
   xyouts, x, y, textoidl("\Omega_m=0.4, \Omega_\Lambda=0"), /device, $
     color=2, charsize=1.4
   xyouts, 2.3, 0.95, "E-dS", orientation=15, chars=1.2
   xyouts, 0.8, 0.9, "Euclidean", orientation=55, chars=1.2

   r = rz(z)
   mplot, z, r, _extra=extra, $
     charsize=1.5

   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-1.) 
   oplot, z, r, color=3
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.6667) 
   oplot, z, r, color=2, lines=2
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.5) 
   oplot, z, r, color=4, lines=3
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.3333) 
   oplot, z, r, color=7, lines=4
   r = rz(z, OmegaM=0.4, OmegaQ=0.6, wQ=-0.1667) 
   oplot, z, r, color=6, lines=5
   oplot, [0,10], [0,10], linesty=1

   xr = [0.55, 0.65]
   yr = [0.4, 0.32, 0.24, 0.16]
   reg2dev, xr, yr[0]*[1,1], x, y
   plots, x, y, linesty=2, color=2, /device
   reg2dev, xr, yr[1]*[1,1], x, y
   plots, x, y, linesty=3, color=4, /device
   reg2dev, xr, yr[2]*[1,1], x, y
   plots, x, y, linesty=4, color=7, /device
   reg2dev, xr, yr[3]*[1,1], x, y
   plots, x, y, linesty=5, color=6, /device
   
   xr = [0.67]
   yr = yr - 0.02
   reg2dev, xr, [yr[0]], x, y
   xyouts, x, y, textoidl("w=-2/3"), /device, color=2, charsize=1.4
   reg2dev, xr, [yr[1]], x, y
   xyouts, x, y, textoidl("w=-1/2"), /device, color=4, charsize=1.4
   reg2dev, xr, [yr[2]], x, y
   xyouts, x, y, textoidl("w=-1/3"), /device, color=7, charsize=1.4
   reg2dev, xr, [yr[3]], x, y
   xyouts, x, y, textoidl("w=-1/6"), /device, color=6, charsize=1.4
   xyouts, 2.3, 0.8, "E-dS (w=0)", orientation=15, chars=1.2
   xyouts, 0.8, 0.9, "Euclidean", orientation=55, chars=1.2
   xyouts, 0.3, 1.3, textoidl("\Omega_m=0.4"), chars=1.5
   xyouts, 0.3, 1.2, textoidl("\Omega_Q=0.6"), chars=1.5
   xyouts, 1.7, 1.08, textoidl("\Lambda (w=-1)"), orient=25, chars=1.2, color=3
   
END


FUNCTION squish_factor, z, $
                        _extra=extra
   
   f = z / (rz(z, _extra=extra) * hubble(z))

   RETURN, f
END


