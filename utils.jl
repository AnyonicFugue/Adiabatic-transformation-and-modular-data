import Combinatorics

function construct_A(nx::Int,ny::Int,Jx::Float64,Jy::Float64,Jz::Float64,xPBC::Bool=true,yPBC::Bool=true)
    Lx=2*nx
    Ly=2*ny
    # The lattice size is Lx \times Ly.

    N=Lx*Ly # The total number of sites
    
    A=zeros(Float64,N,N) # The A matrix encoding interactions.

    for x in 1:Lx
        for y in 1:Ly
            if((x+y)%2==1)
                # The site belongs to the odd sublattice. We should only iterate over even ones here.
                continue
            end

            m=x+(y-1)*Lx
            
            if(x>1)
                # There is a site to the left.
                n_left=x-1+(y-1)*Lx
                A[m,n_left]=Jy
                A[n_left,m]=-Jy # A is antisymmetric.
            else
                n_left=Lx+(y-1)*Lx # There is a site to the left, but it is on the other side of the lattice.
                if(xPBC)
                    A[m,n_left]=Jy # Periodic boundary condition
                    A[n_left,m]=-Jy # A is antisymmetric.
                else
                    A[m,n_left]=-Jy # Anti-periodic boundary condition
                    A[n_left,m]=Jy # A is antisymmetric.
                end
            end

            if(x<Lx)
                # There is a site to the right.
                n_right=x+1+(y-1)*Lx
                A[m,n_right]=Jx
                A[n_right,m]=-Jx # A is antisymmetric.
            else
                n_right=1+(y-1)*Lx # There is a site to the right, but it is on the other side of the lattice.
                if(xPBC)
                    A[m,n_right]=Jx # Periodic boundary condition
                    A[n_right,m]=-Jx # A is antisymmetric.
                else
                    A[m,n_right]=-Jx # Anti-periodic boundary condition
                    A[n_right,m]=Jx # A is antisymmetric.
                end
            end

            if(y<Ly)
                # There is a site below.
                n_down=x+y*Lx
                A[m,n_down]=Jz
                A[n_down,m]=-Jz # A is antisymmetric.
            else
                n_down=x # There is a site above, but it is on the other side of the lattice.
                if(yPBC)
                    A[m,n_down]=Jz # Periodic boundary condition
                    A[n_down,m]=-Jz # A is antisymmetric.
                else
                    A[m,n_down]=-Jz # Anti-periodic boundary condition
                    A[n_down,m]=Jz # A is antisymmetric.
                end
            end
        end
    end

    return A
end

function obtain_occupied_modes(Amat::Array{Float64,2})
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