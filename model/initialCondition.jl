function initialValues()
    u0::Vector{Float64} = zeros(V.len_f_vars);

    u0[V.E1] = 274.0;
    u0[V.E2] = 158.0;
    u0[V.E3] = 294.0;
    u0[V.E4] = 399.0;
    u0[V.G] = 82.4;
    u0[V.S] = 11.5;
    u0[V.I] = 46.4;
    u0[V.R] = 93.6;
    u0[V.T] = 500.0;
    u0[V.O] = 82.3;
    u0[V.A] = 43.1;
    u0[V.P2] = 197;
    u0[V.Akt] = 444.0;
    u0[V.RsD] = 95.7;
    u0[V.Raf] = 743.0;
    u0[V.MEK] = 772.0;
    u0[V.ERK] = 750.0;
    u0[V.ERKPpase] = 35.0;

    return u0
end