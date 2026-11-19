module ActiveDNAWaveguide

export ActiveMotorFrame, DNAViscoelasticMedium, RemodelingFlow, refract

struct ActiveMotorFrame
    l_sys_mm::Float64
    zeta::Float64
    lam_align::Float64
end

struct DNAViscoelasticMedium
    lc_um::Float64
    g_p::Float64
    tau_p::Float64
end

struct RemodelingFlow
    gamma_q::Float64
    d_q::Float64
    ell_v::Float64
end

function refract(sub::ActiveMotorFrame, obj::DNAViscoelasticMedium, act::RemodelingFlow)
    k = pi / sub.l_sys_mm
    visc = 1.0 + (act.ell_v * k)^2
    big_gamma_q = act.gamma_q + act.d_q * k^2
    gamma_act = (sub.lam_align * sub.zeta) / visc - big_gamma_q
    gamma_ve = (1.0 / obj.tau_p) + (obj.g_p / visc)

    tr = gamma_act - gamma_ve
    det_j = big_gamma_q * (obj.g_p / visc) - gamma_act / obj.tau_p
    de_nem = obj.tau_p * big_gamma_q
    de_crit = obj.g_p > 0.0 ? (1.0 + visc / (obj.g_p * obj.tau_p)) : Inf

    if det_j <= 0.0
        return (+1, UInt8(237), de_nem, de_crit)
    elseif tr > 0.0 && tr^2 < 4.0 * det_j
        return (0, UInt8(121), de_nem, de_crit)
    else
        return (-1, UInt8(0), de_nem, de_crit)
    end
end

end # module
