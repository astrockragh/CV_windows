function pofk,karr,gamma=gamma, n=n


	h=0.678
	omega0=0.308
	omegaB=0.022/h/h
	sigma8 = 0.82

if n_elements(n) eq 0 then n = 0.96 ;after Planck
if n_elements(gamma) eq 0 then gamma=omega0*h*exp(-omegaB*(1.+sqrt(2*h)/omega0))

; determine sigma8 from model
; following Peacock, p. 528

	keff=(0.172+0.011*(alog(gamma/0.34))^2)*1.
	q=keff/(h*gamma)
	tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
	sigma8squn=keff^(3.0+n)*tk^2

	q=karr/(h*gamma)   ; k/h Gamma
	tk=alog(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)^2+(5.46*q)^3+(6.71*q)^4)^(-0.25)
	delsq=karr^n*tk^2       ; took out factor of k^3
	delsq=delsq*(sigma8)^2/sigma8squn        ;normalize to Borgani et al. sigma 8, omegam=0.31,omegal=0.69

	pofk=2*!Pi^2*delsq

return,pofk
end

