import Combinatorics
import LinearAlgebra
import CurveFit
import PyPlot

function construct_A_honeycomb(nx::Int,ny::Int,Jx::Float64,Jy::Float64,Jz::Float64,xPBC::Bool=true,yPBC::Bool=true,xShift::Int=0,yShift::Int=0)
    # xShift, yShift are the shifts of the lattice.

    Lx=2*nx
    Ly=2*ny
    # The lattice size is Lx \times Ly.

    N=Lx*Ly # The total number of sites
    
    A=zeros(Float64,N,N) # The A matrix encoding interactions.

    for x in 1:Lx
        for y in 1:Ly
            if((x+(x-1)*yShift+y+(y-1)*xShift)%2==1)
                # The site belongs to the odd sublattice. We should only iterate over even ones here.
                continue
            end

            m=x+(y-1)*Lx # Index of the current even vertex.

            x_sign=1
            y_sign=1
            # The signs to take into account the boudary conditions.

            x_neighbor=0
            y_neighbor=0
            # The coordinates of the neighboring vertex.



            # The interaction with the left-side vertex. yShift needs to be considered.
            if(x>1) 
                # There is a site to the left.
                x_sign=1
                x_neighbor=x-1
            else
                x_sign=2*xPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                x_neighbor=Lx # There is a site to the left, but it is on the other side of the lattice.
            end

            if(y-yShift>0) # The y coordinate of the left-side vertex should be y-yShift mod Ly.
                y_sign=1
                y_neighbor=y-yShift
            else
                y_sign=2*yPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                y_neighbor=y-yShift+Ly # There is a site to the left, but it is on the other side of the lattice.
            end
                
            n_left=x_neighbor+(y_neighbor-1)*Lx
            A[m,n_left]=Jy*x_sign*y_sign
            A[n_left,m]=-Jy*x_sign*y_sign # A is antisymmetric.
            

            
            # The interaction with the right-side vertex. yShift needs to be considered.
            if(x<Lx) 
                # There is a site to the right.
                x_sign=1
                x_neighbor=x+1
            else
                x_sign=2*xPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                x_neighbor=1 # There is a site to the right, but it is on the other side of the lattice.
            end

            if(y+yShift<=Ly) # The y coordinate of the right-side vertex should be y+yShift mod Ly.
                y_sign=1
                y_neighbor=y+yShift
            else
                y_sign=2*yPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                y_neighbor=y+yShift-Ly # There is a site to the left, but it is on the other side of the lattice.
            end
                
            n_left=x_neighbor+(y_neighbor-1)*Lx
            A[m,n_left]=Jx*x_sign*y_sign
            A[n_left,m]=-Jx*x_sign*y_sign # A is antisymmetric.
            


            # The interaction with the down-side vertex. xShift needs to be considered.

            if(y<Ly) 
                # There is a site on the downside.
                y_sign=1
                y_neighbor=y+1
            else
                y_sign=2*yPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                y_neighbor=1 # There is a site on the downside, but it is on the other side of the lattice.
            end

            if(x+xShift<=Lx) # The x coordinate of the down-side vertex should be x+xShift mod Lx.
                x_sign=1
                x_neighbor=x+xShift
            else
                x_sign=2*xPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                x_neighbor=x+xShift-Lx # There is a site to the left, but it is on the other side of the lattice.
            end
                
            n_left=x_neighbor+(y_neighbor-1)*Lx
            A[m,n_left]=Jz*x_sign*y_sign
            A[n_left,m]=-Jz*x_sign*y_sign # A is antisymmetric.
        end
    end

    return A
end

function construct_A_square(nx::Int,ny::Int,mu::Float64,t::Float64,delta::Float64,xPBC::Bool=true,yPBC::Bool=true,xShift::Int=0,yShift::Int=0)
    # xShift, yShift are the shifts of the lattice.

    #=
    mu=-2.0
    t=1.0
    delta=1.0
    =#

    Lx=nx
    Ly=ny

    # The lattice size is Lx \times Ly.

    N_majorana=2*Lx*Ly # The total number of majorana fermions. Note that there are two on each site.
    
    A=zeros(ComplexF64,N_majorana,N_majorana) # The A matrix encoding interactions.

    for x in 1:Lx
        for y in 1:Ly
            
            # Here we do not only iterate over even sites, but also odd sites. But we only process the interactions with the right-side vertex and the down-side vertex.

            m=x+(y-1)*Lx # Index of the current vertex.


            x_sign=1
            y_sign=1
            # The signs to take into account the boudary conditions.

            x_neighbor=0
            y_neighbor=0
            # The coordinates of the neighboring vertex.


            # Self-interaction term
            A[2*m-1,2*m]=-mu/2
            A[2*m,2*m-1]=mu/2



           
            # The interaction with the right-side vertex. yShift needs to be considered.
            if(x<Lx) 
                # There is a site to the right.
                x_sign=1
                x_neighbor=x+1
            else
                x_sign=2*xPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                x_neighbor=1 # There is a site to the right, but it is on the other side of the lattice.
            end

            if(y+yShift<=Ly) # The y coordinate of the right-side vertex should be y+yShift mod Ly.
                y_sign=1
                y_neighbor=y+yShift
            else
                y_sign=2*yPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                y_neighbor=y+yShift-Ly # There is a site to the left, but it is on the other side of the lattice.
            end
                
            n_right=x_neighbor+(y_neighbor-1)*Lx
            A[2*m-1,2*n_right]=x_sign*y_sign*(-t/2-delta/2)
            A[2*n_right,2*m-1]=-x_sign*y_sign*(-t/2-delta/2) # A is antisymmetric.

            A[2*m,2*n_right-1]=x_sign*y_sign*(t/2-delta/2)
            A[2*n_right-1,2*m]=-x_sign*y_sign*(t/2-delta/2) # A is antisymmetric.
            


            # The interaction with the down-side vertex. xShift needs to be considered.

            if(y<Ly) 
                # There is a site on the downside.
                y_sign=1
                y_neighbor=y+1
            else
                y_sign=2*yPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                y_neighbor=1 # There is a site on the downside, but it is on the other side of the lattice.
            end

            if(x+xShift<=Lx) # The x coordinate of the down-side vertex should be x+xShift mod Lx.
                x_sign=1
                x_neighbor=x+xShift
            else
                x_sign=2*xPBC-1 # With PBC the sign is +1, while with APBC the sign is -1.
                x_neighbor=x+xShift-Lx # There is a site to the left, but it is on the other side of the lattice.
            end
                
            n_down=x_neighbor+(y_neighbor-1)*Lx
            A[2*m-1,2*n_down]=x_sign*y_sign*(-t/2-im*delta/2)
            A[2*n_down,2*m-1]=x_sign*y_sign*(t/2-im*delta/2) # A is anti-hermitian.

            A[2*m,2*n_down-1]=x_sign*y_sign*(t/2-im*delta/2)
            A[2*n_down-1,2*m]=x_sign*y_sign*(-t/2-im*delta/2) # A is anti-hermitian.

        end
    end

    return A
end

function obtain_occupied_modes(Amat::Matrix)
    # A is constructed in the above function. We should calculate the eigenvalues and eigenvectors of im*A, which is Hermitian.
    # Exactly half the eigenvalues are positive and half are negative.

    # Since eigen() returns the eigenvalues in ascending order, we can just take the first half of them.
    eigenvectors=eigen(im*Amat).vectors # Note that eigenvectors are vertical columns, not horizontal rows.
    # Eigenstates should correspond to the eigenvectors of iA.
    
    # return first half of the eigenvectors
    return eigenvectors[:,1:Int(size(eigenvectors,2)/2)]

end

function Pfaffian(U::Matrix)
    # Calculate the Pfaffian of the matrix U, which is assumed to be a square complex matrix

    l=size(U,1)
    n=Int(l/2)

    res=0.0+0.0im # The result is complex
    perm_iter=Combinatorics.permutations(1:l)

    for p in perm_iter
        parity=Combinatorics.parity(p)

        prod=1.0+0.0im # The product is complex
        for i in 1:n
            prod*=U[p[2*i-1],p[2*i]]
        end

        if(parity==1)
            # The permutation is even. We should add the corresponding term to the Pfaffian.
            res+=prod
        else
            # The permutation is odd. We should subtract the corresponding term from the Pfaffian.
            res-=prod
        end
    end

    res=res/(2^n*factorial(n))
    return -res
end


function Fit_and_plot(phase_arr::Array{Float64},N_arr::Array{Int})
    # Fit the phase_arr with respect to N_arr. The intercept is the desired phase.
    # Use functions in CurveFit.jl to perform linear fit and those in PyPlot to plot the results.
    a,k=CurveFit.linear_fit(N_arr,phase_arr)
    println("intercept:",a)
    println("slope:",k)

    # Plot phase vs N
    PyPlot.plot(N_arr,phase_arr,"o")
    PyPlot.plot(N_arr,a*N_arr.+k)
    PyPlot.show()
end


function plot_A(A::Matrix)
    N=size(A,1) # The total number of sites.
    Lx=Ly=Int64(sqrt(N)) # The length of the lattice in the x direction. The lattice is assumed to be square.
    
    x_arr=zeros(Float64,N)
    y_arr=zeros(Float64,N)
    # The coordinates of the sites.
    
    for x in 1:Lx
        for y in 1:Ly
            m=x+(y-1)*Lx
            x_arr[m]=x
            y_arr[m]=y
        end
    end
    
    # Plot the interactions.
    for m in 1:N
        for n in 1:N
            if(A[m,n]!=0)
                PyPlot.plot([x_arr[m],x_arr[n]],[y_arr[m],y_arr[n]],color="black")
            end
        end
    end
    
    PyPlot.scatter(x_arr,y_arr,color="red")
    
    PyPlot.xlim(0,Lx+1)
    PyPlot.ylim(0,Ly+1)
    
    PyPlot.show()
end