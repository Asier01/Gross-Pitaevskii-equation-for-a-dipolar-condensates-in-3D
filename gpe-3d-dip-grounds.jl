
import Pkg
Pkg.add("FFTW")
using FFTW

#PARAMETERS:

const n1 = 64
const n2 = 64
const n3 = 64

const dim = [n1,n2,n3]

const ci = 1im
pi

#PHYSICAL constants                 # UNITS

const hbar   = 1.054572E-34         # J*s
const amu    = 1.66053873E-27       # Kg (atomic mass unit)
const a_bohr = 0.5291772083E-10     # m 
const muB    = 9.274009994E-24      #
const mu0    = 4*pi*1E-7            #
const mass1  = 161.9267984*amu      # Kg (162Dy)

#SCALES OF ENERGY/LENGTH               # USE

const l0    = 1E-6                     # To express coordinates in micron
const e0    = hbar^(2/(mass1*l0^2))    # Energy scale (ho correspondence)
const w0    = hbar/(mass1*l0^2)        # Frequency scale (ho correspondence)


#SCATTERING AND ATOMS ------------------------------------------
const nat    = 3.0
const a_init = 96.0                 # a_init is the biggest value on the
const a_fin  = 91.0                 # number of atoms in the system
const steps  = 1000

#a , g11 = 0.0 # Initialization as floats

#DIPOLAR INTERACTION -------------------------------------------

const dipol = true                  # true if there are dipolar interactions 
const dm    = 9.93*muB              # The magnetic momentum of the particles
const Cdd   = mu0*dm^(2)/e0/l0^3    # The dipolar interaction weight

#LHY -----------------------------------------------------------

const lhy_power = 2.5               # Exponent for LHY interactions

# ANGLES (DIPOLE ORIENTATION)

const conversion_angle = pi/180
const phi              = 0.0 * conversion_angle #???
const theta            = 0.0 * conversion_angle #???

# GRID SIZE ---------------------------------------------------- 
                    # In microns
const size_x = 18.0
const size_y = 9.0
const size_z = 9.0
const xmax   = [size_x, size_y, size_z]
const xmin   = -xmax

# TOLERANCE OF GS(Ground State) AND INITIAL TRY CONDITION ------

const rtol = -1E-8
displaced  = false

#ENERGY PARAMETERS  --------------------------------------------

#phi0 = Array{Float64}(undef, n1,n2,n3)
phi0         = zeros(n1,n2,n3)
dphi0        = zeros(n1,n2,n3)

dx           = zeros(3)
dp           = zeros(3)

const nparam = 11
param        = zeros(nparam)
#lhy_coeff = 0.0
#edd = 0.0

#END PARAMETERS:

#MAIN PROGRAM VARIABLES

# POTENTIAL
omega  = 2*pi*[20,67,102]/w0         # frequency of the trap on each direction
uext   = zeros(n1,n2,n3)
psq    = zeros(n1,n2,n3)
potdd  = zeros(n1,n2,n3)

# WAVE FUNCTION
psi0   = zeros(n1,n2,n3)
psi    = zeros(n1,n2,n3)
df     = zeros(n1,n2,n3)
dphi   = zeros(n1,n2,n3)

# GRID
x1 = zeros(n1)
x2 = zeros(n2)
x3 = zeros(n3)

p1 = zeros(n1)
p2 = zeros(n2)
p3 = zeros(n3)


# ENERGY
ener = zeros(6)
ener0 = zeros(6)
#=
mu0 = 0.0
rho = 0.0
rho0 = 0.0
p1max = 0.0
p2max = 0.0
=#
file = ""
fileT = ""
filenumber = ""

i1 = 0
i2 = 0
i3 = 0

iter=0
icount=0 
i = 0

#END MAIN PROGRAM VARIABLES

#SETTING UP FOR FFT

plan_f = Int32(0)
plan_b = Int32(0)
iret = 6
nthreads = iret

println("NUMBER OF THREADS FOR FFT : " * string(nthreads))

    #SET NUMBER OF NTHREADS FOR USE ON ON FFT
FFTW.set_num_threads(nthreads)

    #PLANS FOR FOURIER TRANSFORM AND INVERSE FOURIER TRANSFORM
fft_plan = FFTW.plan_fft(dim , flags = FFTW.MEASURE) 
ifft_plan = FFTW.plan_ifft(dim, flags = FFTW.MEASURE)

#END SETTING UP FOR FFT



# FUNCTIONS / SUBROUTINES

function init_sys(a::Float64) # A :: Float64

    g11=4*pi*(a*a_bohr/l0)*nat
    edd = Cdd/(12*pi*(a*a_bohr/l0))
    lhy_coeff = (256*sqrt(pi)/15)*(a*a_bohr/l0)^2.5*(1 + 1.5*edd^2)*nat^1.5

    return lhy_coeff # ?????????

end

function init_stats()

    p_vec = zeros(3)
    u_B = zeros(3)

    dx = (xmax-xmin)./dim ##/dim?????
    dp = 2*pi ./ (xmax-xmin) #Element wise??
    
 #   WIP
 #   forall(i1=1:n1) x1(i1) = xmin(1) + (i1-1)*dx(1)
 #   forall(i1=1:n1/2) p1(i1) = dp(1)*(i1-1)
 #   forall(i1=n1/2+1:n1) p1(i1) = dp(1)*(i1-1-n1)
#
 #   forall(i2=1:n2) x2(i2) = xmin(2) + (i2-1)*dx(2)
 #   forall(i2=1:n2/2) p2(i2) = dp(2)*(i2-1)
 #   forall(i2=n2/2+1:n2) p2(i2) = dp(2)*(i2-1-n2)
#
 #   forall(i3=1:n3) x3(i3) = xmin(3) + (i3-1)*dx(3)
 #   forall(i3=1:n3/2) p3(i3) = dp(3)*(i3-1)
 #   forall(i3=n3/2+1:n3) p3(i3) = dp(3)*(i3-1-n3)

 
 #    forall(i1=1:n1,i2=1:n2,i3=1:n3)
 #    uext(i1,i2,i3) = 0.5d0*((omega(1)*x1(i1))**2. + (omega(2)*x2(i2))**2. + (omega(3)*x3(i3))**2.)
 #   end forall
 #   
 #   forall(i1=1:n1,i2=1:n2,i3=1:n3) psq(i1,i2,i3) = p1(i1)**2 + p2(i2)**2 + p3(i3)**2

 #  ! dipolar potential in momentum space
 #  u_B(:) = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
#   
 #  allocate(ca(n1,n2,n3))
#   
 #  do i1 = 1, n1
 #     do i2 = 1, n2
 #        do i3 = 1, n3
 #           p_vec(:) = [p1(i1),p2(i2),p3(i3)]
 #           if(norm2(p_vec).eq.0) then
 #              potdd(i1,i2,i3) = -1.d0/3.d0
 #           else
 #              ca(i1,i2,i3) = dot_product(u_B,p_vec)/sqrt(psq(i1,i2,i3))
 #              potdd(i1,i2,i3) = ca(i1,i2,i3)**2. - 1.d0/3.d0
 #           end if
 #        end do
 #     end do
 #  end do
#   
 #  !regularization in 3D
 #  where(psq.gt.0)
 #     potdd = (1.d0 + 3.d0*cos(R*sqrt(psq))/(R**2.*psq) - 3.d0*sin(R*sqrt(psq))/(R**3.*psq**1.5))*potdd
 #  end where
#   
 #  deallocate(ca)
#   

end

function dipolar(phi, u_dipole)
    #what is phi?
    phi = rand(ComplexF64, 64, 64, 64) 
    phi[1,2,2]
    #

    df = abs(phi)
    df = FFTW.fft(df)
    df = potdd * df
    
    u_dipole = FFTW.ifft(df)
    
end


    # END FUNCTIONS / SUBROUTINES

#MAIN PROGRAM



#END MAIN PROGRAM


  