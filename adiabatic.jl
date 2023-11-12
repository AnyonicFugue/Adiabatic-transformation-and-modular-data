include("utils.jl")

function adiabatic(A0::Array{Float64,2},A1::Array{Float64,2},nstep::Float64)

    # A_0 is the initial Hamiltonian, and the final Hamiltonian is A1.
    # nstep is number of steps. We cannot perform a conditinuous adiabatic evolution, so we need to discretize it.

    dA=(A1-A0)/nstep # The increment of A at each step.
    res=0.0

    A0_occupied_modes=obtain_occupied_modes(A0)
    last_occupied_modes=A0_occupied_modes # The occupied modes at the last step.

    overlap_with_A0=zeros(ComplexF64,nstep)
    phase_arr=zeros(Float64,nstep)

    for i in range(1,nstep)
        A=A0+i*dA
        A_occupied_modes=obtain_occupied_modes(A)

        # Use the formula $\begin{aligned} \mathrm{d} \Theta= & \arg \langle\Psi(\lambda) \mid \Psi(\lambda+\mathrm{d} \lambda)\rangle \\ & -\arg \left\langle\Psi(\lambda) \mid \Psi_{\text {ref }}\right\rangle+\arg \left\langle\Psi(\lambda+\mathrm{d} \lambda) \mid \Psi_{\text {ref }}\right\rangle .\end{aligned}$

        overlap_with_A0[i]=Pfaffian(A_occupied_modes'*A0_occupied_modes)

        if(i==1)
            res=res+0.0
        else
            phase_arr[i]=angle(Pfaffian(last_occupied_modes'*A_occupied_modes))-angle(overlap_with_A0[i-1])+angle(overlap_with_A0[i])
            res=res+phase_arr[i]
        end

        last_occupied_modes=A_occupied_modes
    end

    return res,phase_arr,overlap_with_A0 # Return the Berry phase together with the overlaps with the initial state.

end