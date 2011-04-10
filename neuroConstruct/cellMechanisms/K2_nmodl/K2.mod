TITLE K2 calcium-activated potassium current
: Calcium activated K channel.
COMMENT
  from "An Active Membrane Model of the Cerebellar Purkinje Cell
        1. Simulation of Current Clamp in Slice"
ENDCOMMENT

UNITS {
    (molar) = (1/liter)
}

UNITS {
    (mV) =  (millivolt)
    (mA) =  (milliamp)
    (mM) =  (millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX %Name%
    USEION ca READ cai
    USEION k WRITE ik
    RANGE gmax,gk,zinf,ik
}


PARAMETER {
    celsius=37  (degC)
    v       (mV)
    gmax=%Max Conductance Density%    (mho/cm2)   : Maximum Permeability
    cai     (mM)
    ek  = -85   (mV)
    dt      (ms)
}


ASSIGNED {
    ik      (mA/cm2)
    minf
    mexp
    zinf
    zexp
    gk
}

STATE { m z }       : fraction of open channels

BREAKPOINT {
    SOLVE state
    gk = gmax*m*z*z
    ik = gk*(v - ek)
}
:UNITSOFF
:LOCAL fac

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.
PROCEDURE state() { : exact when v held constant; integrates over dt step

    ::::::::::::ratezinf(floor((cai - 4e-5)/9.9986667e-5)*9.9986667e-5 + 4e-5)
    ratezinf(cai)
    rate(v)
    m = m + mexp*(minf - m)
    z = z + zexp*(zinf - z)
    VERBATIM
    return 0;
    ENDVERBATIM
}

INITIAL {
    ratezinf(cai)
    rate(v)
    m = minf
    z = zinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { :callable from hoc
    alp = 0.2/(ca*1000)
}

FUNCTION bet(v (mV)) (1/ms) { :callable from hoc
    bet = 0.075/exp((v+25)/6)
}

PROCEDURE rate(v (mV)) { :callable from hoc
    LOCAL a,b
    TABLE zexp, minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200
    zexp = (1 - exp(-dt/10))
    b = bet(v)
    minf = 25/(25+b)
    mexp = (1 - exp(-dt*(25+b)))
}
PROCEDURE ratezinf(ca (mM)) { :callable from hoc
    LOCAL a,b
    a = alp(0,ca)
    zinf = 1/(1+a)
}
:UNITSON
