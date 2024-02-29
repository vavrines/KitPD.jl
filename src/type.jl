"""
Peridynamic matter properties

"""
@with_kw mutable struct PDMater{T} <: KitBase.AbstractProperty
    δ::T
    rh::T = 3.015 * δ
    emod::T = 2e11 # modulus Pa
    ν::T = 0.333 # poison ratio
    kbk::T = emod / (2 - 2 * ν) # bulk modulus
    ksh::T = emod / (2 + 2 * ν) # shear modulus
    c::T = 9.0 * emod / (π * (rh^3) * δ) # micro-modulus in PD
    kc::T = 51.9 # thermal conductivity W/mK
    kp::T = 6.0 * kc / (π * (rh^3) * δ) # micro-conductivity in PD
    aph::T = 11.5e-6 # thermal expansion K-1
    cᵥ::T = 472.0 # specific heat capacity J/kgK
    ρ::T = 7870.0 # kg/m^3
    gc::T = 42320.0 # critial energy release rate
    sc::T = sqrt(gc / (rh * (ksh * 6.0 / π + 16.0 / (9.0 * (π^2)) * (kbk - 2.0 * ksh))))
end
