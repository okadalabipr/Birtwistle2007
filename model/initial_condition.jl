function initial_values()
    u0::Vector{Float64} = zeros(length(F_V));

    u0[E1] = 274.0;
    u0[E2] = 158.0;
    u0[E3] = 294.0;
    u0[E4] = 399.0;
    u0[G] = 82.4;
    u0[S] = 11.5;
    u0[I] = 46.4;
    u0[R] = 93.6;
    u0[T] = 500.0;
    u0[O] = 82.3;
    u0[A] = 43.1;
    u0[P2] = 197;
    u0[Akt] = 444.0;
    u0[RsD] = 95.7;
    u0[Raf] = 743.0;
    u0[MEK] = 772.0;
    u0[ERK] = 750.0;
    u0[ERKPpase] = 35.0;

    return u0
end