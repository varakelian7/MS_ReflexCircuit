NEURON {
    POINT_PROCESS NMDA
    USEION ca READ cai, cao
    NONSPECIFIC_CURRENT i
    RANGE gbar, e, tau_r, tau_d, mg, factor
}

PARAMETER {
    gbar   = 0.001    (uS)
    e      = 0        (mV)
    tau_r  = 2        (ms)
    tau_d  = 100      (ms)
    mg     = 1        (mM)     : extracellular Mg2+
}

ASSIGNED {
    v       (mV)
    i       (nA)
    g       (uS)
    factor  (1)
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    A = B = 0
    tp = (tau_r*tau_d)/(tau_d - tau_r)*log(tau_d/tau_r)
    factor = 1/(-exp(-tp/tau_r) + exp(-tp/tau_d))
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = (B - A)
    i = g*(v - e) / (1 + mg * exp(-0.062*v)/3.57)
}

DERIVATIVE state {
    A' = -A/tau_r
    B' = -B/tau_d
}

NET_RECEIVE(weight) {
    A = A + weight*factor
    B = B + weight*factor
}
