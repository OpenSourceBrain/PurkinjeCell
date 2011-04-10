TITLE Fast sodium current
 
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
    USEION na READ ena WRITE ina             :////PG///  Note!!! originally no 'READ ena' here, so it used 45 (below), ignoring external ena
        RANGE  gmax, gna, minf, hinf, mexp, hexp, ina
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        dt (ms)
        gmax	= %Max Conductance Density% (mho/cm2)
        ena	= 45 (mV)

}
 
STATE {
        m h
}

ASSIGNED {
        ina (mA/cm2)
        gna minf hinf mexp hexp 
}
 
BREAKPOINT {
        SOLVE states
        gna = gmax *m*m* m*h 
	ina = gna* (v-ena)
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
                :"m" sodium activation system
        alpha = 35/exp((v+5)/(-10))
        beta =  7/exp((v+65)/20)
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
                :"h" sodium inactivation system
        alpha = 0.225/(1+exp((v+80)/10))
        beta = 7.5/exp((v-3)/(-18))
        sum = alpha + beta
        hinf = alpha/sum
        hexp = 1 - exp(tinc*sum)
}

 
UNITSON

