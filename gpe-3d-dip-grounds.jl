
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
const phi              = 0.0 * conversion_angle
const theta            = 0.0 * conversion_angle

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


#MAIN PROGRAM


#END MAIN PROGRAM


  