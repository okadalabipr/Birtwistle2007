using ODE;

include("./model/setParamConst.jl");
include("./model/setVarEnum.jl");
include("./model/diffeq.jl");
include("./model/initialValues.jl");

p = setParamConst();
u0 = initialValues();

t = 0.0:1.0:1800.0;
condition = 8;

ERK_act = zeros(length(t),condition);
Akt_act = zeros(length(t),condition);

for i=1:condition
    if i==1
        u0[E] = 0.0;
        u0[H] = 0.5;
    elseif i==2
        u0[E] = 0.0;
        u0[H] = 10.0;
    elseif i==3
        u0[E] = 0.5;
        u0[H] = 0.0;
    elseif i==4
        u0[E] = 0.5;
        u0[H] = 0.5;
    elseif i==5
        u0[E] = 0.5;
        u0[H] = 10.0;
    elseif i==6
        u0[E] = 10.0;
        u0[H] = 0.0;
    elseif i==7
        u0[E] = 10.0;
        u0[H] = 0.5;
    elseif i==8
        u0[E] = 10.0;
        u0[H] = 10.0;
    end

    (sol_t,sol_u) = ode23s(diffeq,u0,t;points=:specified);

    for j=1:length(t)
        ERK_act[j,i] = sol_u[j][ERKstar] + sol_u[j][pERK_ERKPpase];
        Akt_act[j,i] = sol_u[j][Aktstar];
    end

end
