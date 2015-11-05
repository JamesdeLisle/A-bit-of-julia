using PyCall
@pyimport matplotlib.pyplot as plt

type parameterspace
    
    mu_disc::Int
    Delta_disc::Int
    p_disc::Int
    mu::LinSpace{Float64}
    Delta::LinSpace{Float64}
    p::LinSpace{Float64}

    function parameterspace{T<:Number}(mu_min::T,mu_max::T,mu_disc::Int,Delta_min::T,Delta_max::T,Delta_disc::Int,p_min::T,p_max::T,p_disc::Int)
        mu_space = linspace(mu_min,mu_max,mu_disc)
        Delta_space = linspace(Delta_min,Delta_max,Delta_disc)
        p_space = linspace(p_min,p_max,p_disc)
        new(mu_disc,Delta_disc,p_disc,mu_space,Delta_space,p_space)
    end

end

function hamiltonian(mu::Float64,Delta::Float64,t::Float64,p::Float64)
    
    H_00 = mu + cos(p)  
    H_01 = 1im * Delta * sin(p)
    
    ham = [H_00 H_01
    conj(H_01) -H_00]

    return ham

end

function getGap(mu::Float64,Delta::Float64,space::parameterspace)
    
    gap = 1.0e4 
    for (p) in space.p
        evals = sort(eigvals(hamiltonian(mu,Delta,1.0,p)))
        if evals[2] < gap
            gap = evals[2]
        end
    end
    
    return gap

end

function scan(space::parameterspace)
    
    eigenArray = zeros(space.mu_disc,space.Delta_disc)
    
    for (mu_idx,mu) in enumerate(space.mu)
        for (Delta_idx,Delta) in enumerate(space.Delta)
            #println((float(mu_idx)*float(space.Delta_disc)+float(Delta_idx))/(float(space.mu_disc)*float(space.Delta_disc)))
            eigenArray[mu_idx,Delta_idx] = getGap(mu,Delta,space)
        end
    end
    
    return eigenArray

end

space = parameterspace(-2.0,2.0,1000,-4.0,4.0,1000,0.0,2*3.14159,100)

gap = scan(space)
plt.pcolor(gap)
plt.show()

