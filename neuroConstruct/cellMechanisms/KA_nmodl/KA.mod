TITLE A current
 
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
        USEION k READ ek WRITE ik             :////PG///  Note!!! originally no 'READ ek' here, so it used -85 (below), ignoring external ek
        RANGE  gmax, gk, minf, hinf, mexp, hexp, ik
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gmax    = %Max Conductance Density% (mho/cm2)
        ek  = -85 (mV)

}
 
STATE {
        m h
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gk = gmax *m*m*m* m*h 
        ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
}

PROCEDURE states() {  :Computes state variables m, h
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
                      
        LOCAL  q10, tinc, alpha, beta, sum
        TABLE minf, mexp, hinf, hexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
        
        
                :"m" potassium activation system
        alpha = 1.4/(1+exp((v+27)/(-12)))
        beta =  0.49/(1+exp((v+30)/4))
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)  
        
        
                :"h" potassium inactivation system
        alpha = 0.0175/(1+exp((v+50)/8))
        beta = 1.3/(1+exp((v+13)/(-10)))
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON

