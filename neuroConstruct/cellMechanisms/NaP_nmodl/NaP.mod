TITLE Persistent sodium current
 
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
	RANGE gmax, gna, minf, ina
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	v		(mV)
	celsius	= 37	(degC)
	ena	= 45	(mV)
	gmax	= %Max Conductance Density% (mho/cm2)
    dt
}
 
STATE {
        m 
}
 
ASSIGNED {
        ina (mA/cm2)
 	    minf
        mexp
        gna
}
 
BREAKPOINT {
        SOLVE states
       
        gna = gmax*m*m*m
        ina = gna*(v - ena)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
}

PROCEDURE states() {  :Computes states variable m
        rates(v)      :             at the current v and dt.
        m = m + mexp*(minf-m)

}
 
PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  q10, tinc, alpha, beta, sum

        TABLE minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200

        q10 = 3^((celsius - 37)/10)
        tinc = -dt * q10
                :"m" sodium activation system
        alpha = 200/(1+exp((v-18)/(-16)))
        beta =  25/(1+exp((v+58)/8))
        sum = alpha + beta
        minf = alpha/sum
        mexp = 1 - exp(tinc*sum)
}
 

 
UNITSON

