TITLE Persistent potassium current
 
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
        USEION k  WRITE ik
        RANGE  gmax, gk, minf,  mexp, ik
} 

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gmax   = %Max Conductance Density% (mho/cm2)
        ek      = -85 (mV)

}
 
STATE {
        m 
}
 
ASSIGNED {
        ik (mA/cm2)
        gk minf  mexp 
}
 
BREAKPOINT {
        SOLVE states
        gk = gmax *m
        ik = gk* (v-ek)
}
 
UNITSOFF
 
INITIAL {
        rates(v)
        m = minf

}

PROCEDURE states() {  :Computes state variables m
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)
}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
                      
        LOCAL  q10, tinc, sum
        TABLE minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
        
        
                :"m" potassium activation system
:sum = 3.3*(exp((v+35)/40)+exp(-(v+35)/20))/200
:sum = 3.3*(exp((v+35)/20)+exp(-(v+35)/20))/200
:changed possition of brackets to mimick genesis
        sum = (3.3*exp((v + 35)/20) + exp(-(v+35)/20))/200
        minf = 1.0 / (1+exp(-(v+35)/10))
        mexp = 1 - exp(tinc*sum)
               
}

 
UNITSON
