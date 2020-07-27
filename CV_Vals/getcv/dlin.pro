function Dlin, Omega_0, z
;assumes a flat cosmology
  Omega = Omega_0*(1+z)^3/(1-Omega_0 + (1+z)^3*Omega_0)
  gz = 2.5*Omega/(1./70. + 209.0*Omega/140. + Omega^(4./7.))
  g0 = 2.5*Omega_0/(1./70. + 209.0*Omega_0/140. + Omega_0^(4./7.))
  D = gz/(g0*(1+z))
  return, D
end

