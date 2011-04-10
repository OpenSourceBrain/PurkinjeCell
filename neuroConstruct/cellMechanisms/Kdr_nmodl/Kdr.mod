TITLE Delayed rectifire
 
COMMENT
  from "An Active Membrane Model of the Cerebellar Purkinje Cell
        1. Simulation of Current Clamp in Slice"

Taken from De Schutter model conversion from GENESIS by Jenny Davie, Arnd Roth,
   Volker Steuber, Erik De Schutter & Michael Hausser 28.8.2004

ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX %Name%
        USEION k WRITE ik
        RANGE  gmax, gk, minf, hinf, mexp, hexp, ik, alpha2, beta2, gamma2, zeta2
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gmax = %Max Conductance Density% (mho/cm2)
        ek   = -85 (mV)

}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp alpha2 beta2 gamma2 zeta2
}
 
BREAKPOINT {
        SOLVE states
        gk = gmax *m*m*h
        ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
    trates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  :Computes state variables m,h
        trates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}

PROCEDURE trates(v) {
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        rates(v)
}


PROCEDURE rates(vm) {  :Computes rate and other constants at current v.
                       :Call once from HOC to initialize inf at resting v.
                      
        LOCAL  q10, tinc, tauh, alpha, beta, gamma, zeta, taum
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10

                :"m" potassium activation system
    alpha = trap0(vm,8,23.5e-3,12)
    beta = 5*exp((vm+127)/(-30))
    minf = alpha/(alpha+beta)
    
    gamma = trap0(vm,-12,23.5e-3,12)
    zeta = 5*exp((vm+147)/(-30))
    taum = 1/(gamma + zeta)
    
    
    mexp = 1 - exp(tinc/taum)

                :"h" potassium activation system
        hinf = 1.0 / (1+exp((vm+25)/4))
        if(vm<-25) {
    tauh = 1200
    }else{
    tauh = 10
    }
    hexp = 1 - exp(tinc/tauh)
    alpha2 = alpha
    beta2 = beta
    gamma2 = gamma
    zeta2 = zeta

}

FUNCTION trap0(vn,th,a,q) {
    if (fabs(vn-th) > 1e-6) {
            trap0 = a * (vn - th) / (1 - exp(-(vn - th)/q))
    } else {
            trap0 = a * q
    }
}

 
UNITSON

