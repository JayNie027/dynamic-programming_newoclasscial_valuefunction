using Printf

sigma = 1;

g = 0;

n = 0;

alpha = 1/3;

delta = 0.1;

beta = 0.96 ;

#predetermined steady state kstar

ks = ( (1/alpha)*((1+n)*(1+g)/beta-1+delta) )^(1/(alpha-1))

N=150


K=collect(LinRange(.85*ks,1.05*ks, N)) 
# K=collect(LinRange(.1*ks,1.5*ks, N))    # State space

KC = repeat(K,1,N)

KT = KC';

Y = KC.^alpha;

C = Y + (1-delta)*KC - KT*(1+n)*(1+g);

# C[C .<= 0] .= 0.0001;

UU = log.(C);

UU[C .<= 0] .= -10000;

# vv = zeros(size(K,1),1)    # initial gusss for value function 150 * 1
vv = zeros(N,1)  
new_vv = zeros(N,1)
ind = zeros(N,1)
polk = zeros(N,1)

maxiter = 1000;

iter = 0;
eps = 10^(-6)

for iter in 1:maxiter 

    iter += 0

    vv = copy(new_vv)

    for i in 1:N

        new_vv[i] = maximum( UU[i, :]' + beta .* vv' )
    
        ind[i] = argmax(UU[i, :]' + beta .* vv')[2]
    
        polk[i] = K[convert.(Int, ind[i])]
    
    end

    diff = sqrt( sum( (vv - new_vv).^2 ))

    if diff < eps
        break 
    end

    # vv = copy(new_vv)

    @printf "%d interations, gap: %f \n" iter diff

end

new_vv


using Plots

plot(KC[:, 1], polk)
Plots.abline!(1,0)
vline!([ks])

plot(KC[:, 1], new_vv)