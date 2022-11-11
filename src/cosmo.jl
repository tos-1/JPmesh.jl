#############################################################################
#cosmology, D1 and EisHu from https://github.com/komatsu5147/MatterPower.jl
#############################################################################

"Growth function linear order, defined in MatterPower"
function D1(a, Om_m, Om_l)
    sol = setup_growth(Om_m, Om_l)
    sol(a)[1]
end

"returns Pk a function of k (h/Mpc) at chosen z(a_ini) using Eisenstein and Hu fitting formula, 
 already normalized to have the chosen s8. No wiggle power spectrum, defined in MatterPower"
function EisHu(simbox)
   pk(kovh) = simbox.A *#* D1(simbox.a_i,simbox.Om_m, simbox.Om_l)^2 * 
   (kovh * simbox.h / simbox.k_p)^(simbox.ns - 1) *
   (2 * kovh^2 * 2998^2 / 5 / simbox.Om_m)^2 *
   t_nowiggle(kovh * simbox.h, simbox.Om_m * simbox.h^2, simbox.Om_b/simbox.Om_m)^2 *
   2 * π^2 / kovh^3
   s8_A = sigma2(pk, 8)
   pks8(kovh) = pk(kovh) * (simbox.s8/s8_A) * ( D1(simbox.a_i,simbox.Om_m, simbox.Om_l) / D1(1.0,simbox.Om_m, simbox.Om_l) )^2 
   return pks8
end

"closure equation"
function Om_k(Om_m, Om_l)
    return 1-Om_m-Om_l
end

"adimensional Hubble E = H/H0"
function E(a,Om_l,Om_m)
    sqrt( Om_l + Om_m * a^-3 + Om_k(Om_m,Om_l) * a^-2 )
end

"Omega as a function of scale factor"
function Om0z(a, Om_l, Om_m)
    Om_m/a^3/E(a,Om_l,Om_m)^2.
end

"Omega_Lambda as a function of scale factor"
function OmLz(a, Om_l, Om_m)
    Om_l/E(a,Om_l,Om_m)^2.
end

"kernel of the leap-frog integrator"
function Fa(a,Om_l, Om_m)
    sqrt(a/(Om_m+Om_k(Om_m,Om_l)*a+Om_l*a^3))
end

"Growth function from Carrol+92. At early times d1~a"
function D1fit(a, Om_m, Om_l)
    5/2.0 * Om0z(a, Om_l, Om_m)*a / ( Om0z(a, Om_l, Om_m)^(4/7.) - OmLz(a, Om_l, Om_m) + (1.0 + Om0z(a, Om_l, Om_m)/2.)*(1. + OmLz(a, Om_l, Om_m)/70.) )
end
