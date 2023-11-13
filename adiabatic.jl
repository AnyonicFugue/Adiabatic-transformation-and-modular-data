include("utils.jl")
import LinearAlgebra.det as det

function adiabatic(A0::Matrix,A1::Matrix,nstep::Int,debug::Bool=false)

    # A_0 is the initial Hamiltonian, and the final Hamiltonian is A1.
    # nstep is number of steps. We cannot perform a conditinuous adiabatic evolution, so we need to discretize it.

    dA=(A1-A0)/nstep # The increment of A at each step.
    res=0.0

    A0_occupied_modes=obtain_occupied_modes(A0)
    last_occupied_modes=A0_occupied_modes # The occupied modes at the last step.

    overlap_with_A0=zeros(ComplexF64,nstep)
    phase_arr=zeros(Float64,nstep)

    for i in range(1,nstep)
        if(debug)
            println("step",i)
        end

        A=A0+i*dA
        A_occupied_modes=obtain_occupied_modes(A)

        if(debug)
            println("Calculated occupied modes")
        end

        # Brute-force calculation of Pfaffian is doomed. Make use of the determinant to clever calculate it.
        # Note that the Pfaffian is the square root of the determinant. We only need to determine the sign.



        if(debug)
            println("Calculated overlap")
        end

        if(i==1)
            res=res+0.0
            tmp=sqrt(det(A0_occupied_modes'*A_occupied_modes))
            
            # Take the value closer to 1, since A only differs from A0 by a small amount.
            if(abs(tmp-1.0)<abs(tmp+1.0))
                overlap_with_A0[i]=tmp
            else
                overlap_with_A0[i]=-tmp
            end
            
        else
            # Calculate the overlap with A0
            tmp=sqrt(det(A0_occupied_modes'*A_occupied_modes))
            # Take the value closer to the last overlap.
            if(abs(tmp-overlap_with_A0[i-1])<abs(tmp+overlap_with_A0[i-1]))
                overlap_with_A0[i]=tmp
            else
                overlap_with_A0[i]=-tmp
            end

            # Calculate the overlap with the previous A

            tmp=sqrt(det(last_occupied_modes'*A_occupied_modes))
            # Take the value closer to 1.
            if(abs(tmp-1)<abs(tmp+1))
                tmp=tmp
            else
                tmp=-tmp
            end

            phase_arr[i]=angle(tmp)-angle(overlap_with_A0[i-1])+angle(overlap_with_A0[i])
            res=res+phase_arr[i]
        end

        last_occupied_modes=A_occupied_modes
    end

    return res,phase_arr,overlap_with_A0 # Return the Berry phase together with the overlaps with the initial state.

end

function adiabatic_ver1(A0::Matrix,A1::Matrix,nstep::Int,debug::Bool=false)

    # A_0 is the initial Hamiltonian, and the final Hamiltonian is A1.
    # nstep is number of steps. We cannot perform a conditinuous adiabatic evolution, so we need to discretize it.

    dA=(A1-A0)/nstep # The increment of A at each step.
    res=0.0

    A0_occupied_modes=obtain_occupied_modes(A0)
    last_occupied_modes=A0_occupied_modes # The occupied modes at the last step.

    overlap_with_A0=zeros(ComplexF64,nstep)
    phase_arr=zeros(Float64,nstep)

    for i in range(1,nstep)
        if(debug)
            println("step",i)
        end

        A=A0+i*dA
        A_occupied_modes=obtain_occupied_modes(A)

        if(debug)
            println("Calculated occupied modes")
        end

        # Brute-force calculation of Pfaffian is doomed. Make use of the determinant to clever calculate it.
        # Note that the Pfaffian is the square root of the determinant. We only need to determine the sign.

        if(debug)
            println("Calculated overlap")
        end

        if(i==1)
            res=res+0.0

            tmpmat=A0_occupied_modes'*A_occupied_modes
            tmpmat=(tmpmat-tmpmat')/2.0 # Make it skew-symmetric.
            tmp=sqrt(det(tmpmat))
            
            # Take the value closer to 1, since A only differs from A0 by a small amount.
            if(abs(tmp-1.0)<abs(tmp+1.0))
                overlap_with_A0[i]=tmp
            else
                overlap_with_A0[i]=-tmp
            end
            
        else
            # Calculate the overlap with A0
            tmpmat=A0_occupied_modes'*A_occupied_modes
            tmpmat=(tmpmat-tmpmat')/2.0 # Make it skew-symmetric.
            tmp=sqrt(det(tmpmat))
            
            # Take the value closer to the last overlap.
            if(abs(tmp-overlap_with_A0[i-1])<abs(tmp+overlap_with_A0[i-1]))
                overlap_with_A0[i]=tmp
            else
                overlap_with_A0[i]=-tmp
            end

            # Calculate the overlap with the previous A
            tmpmat=last_occupied_modes'*A_occupied_modes
            tmpmat=(tmpmat-tmpmat')/2.0 # Make it skew-symmetric.
            tmp=sqrt(det(tmpmat))
            
            # Take the value closer to 1.
            if(abs(tmp-1)<abs(tmp+1))
                tmp=tmp
            else
                tmp=-tmp
            end

            phase_arr[i]=angle(tmp)-angle(overlap_with_A0[i-1])+angle(overlap_with_A0[i])
            res=res+phase_arr[i]
        end

        last_occupied_modes=A_occupied_modes
    end

    return res,phase_arr,overlap_with_A0 # Return the Berry phase together with the overlaps with the initial state.

end

