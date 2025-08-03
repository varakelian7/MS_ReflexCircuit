NEURON {
    POINT_PROCESS GABAB
    NONSPECIFIC_CURRENT i
    RANGE gbar, e, tau_r, tau_d, factor
}

PARAMETER {
    gbar  = 0.0005   (uS)
    e     = -95      (mV)
    tau_r = 50       (ms)
    tau_d = 300      (ms)
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
    i = g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau_r
    B' = -B/tau_d
}

NET_RECEIVE(weight) {
    A = A + weight*factor
    B = B + weight*factor
}
