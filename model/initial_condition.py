from .name2idx import variables as V

def initial_values():
    y0 = [0]*V.len_f_vars

    y0[V.E1] = 274.0
    y0[V.E2] = 158.0
    y0[V.E3] = 294.0
    y0[V.E4] = 399.0
    y0[V.G] = 82.4
    y0[V.S] = 11.5
    y0[V.I] = 46.4
    y0[V.R] = 93.6
    y0[V.T] = 500.0
    y0[V.O] = 82.3
    y0[V.A] = 43.1
    y0[V.P2] = 197
    y0[V.Akt] = 444.0
    y0[V.RsD] = 95.7
    y0[V.Raf] = 743.0
    y0[V.MEK] = 772.0
    y0[V.ERK] = 750.0
    y0[V.ERKPpase] = 35.0

    return y0