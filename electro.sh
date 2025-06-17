#!/bin/bash
#===============================================================================
# VASP Electronic Structure Workflow - ZnO
#===============================================================================

# Configuration
KPOINTS_MESH="5 5 1"
#MAGMOM='2*2.0 2*2.0'


BASE_INCAR="#============================
LWAVE = .TRUE.                          # Write WAVECAR file (wavefunction)
LCHARG = .TRUE.                         # Write CHGCAR file (charge density)
LREAL = Auto                            # Real space projection (Auto for >20 atoms)
EDIFF = 1E-5                           # Electronic convergence criterion (eV)
EDIFFG = -0.01                         # Ionic convergence criterion (eV/Å)
SIGMA = 0.05                           # Smearing width (eV)
POTIM = 0.4                            # Time step for ionic motion (fs)
PREC = Accurate                        # Precision level (Normal/Accurate)
ALGO = Fast                            # Electronic minimization algorithm
#NCORE = 2                              # Number of cores per orbital
LMAXMIX = 4                            # Maximum l-quantum number for charge density mixing
#MAGMOM = ${MAGMOM} ;ISPIN=2                    # Initial magnetic moments per atom
IDIPOL = 3          # Apply correction along surface normal (z-axis)
LDIPOL = .TRUE.     # Enable dipole potential correction
#DIPOL = 0.5 0.5 0.5 # Center of cell in fractional coordinates
 "

KPATH="HEX (hexagonal) G-M-K-G path from BandLines
40   ! Total number of k-points (1 + 50 + 50 + 50)
Line-mode
reciprocal
   0.0000   0.0000   0.0000   ! \Gamma
   0.5000   0.0000   0.0000   ! M
 
   0.5000   0.0000   0.0000   ! M
   0.3333   0.3333   0.0000   ! K
 
   0.3333   0.3333   0.0000   ! K
   0.0000   0.0000   0.0000   ! \Gamma
"

# Helper functions
create_kpoints() {
    cat > KPOINTS << EOF
K-Points
0
Monkhorst Pack
${KPOINTS_MESH}
0 0 0
EOF
}

copy_scf_files() {
    cp ../2.Self_consistent/{CONTCAR,CHGCAR,WAVECAR} .
    cp CONTCAR POSCAR
    cp ../POTCAR .
}

#===============================================================================
# 1. Geometry Optimization
#===============================================================================
mkdir 1.Geometry_optimization && cd 1.Geometry_optimization
cp ../{POSCAR,POTCAR} .
create_kpoints

cat > INCAR << EOF
ISTART = 0                             # Start from scratch (no WAVECAR)
ICHARG = 2                             # Calculate charge density from atomic densities
ISMEAR = 0
NSW = 1500                             # Maximum number of ionic steps
IBRION = 2                             # Conjugate gradient ionic relaxation algorithm
ISIF = 2                               # Relax ions and cell volume/shape
${BASE_INCAR}
EOF

gerun vasp_std

#===============================================================================
# 2. Self-Consistent Field Calculation
#===============================================================================
cd .. && mkdir 2.Self_consistent && cd 2.Self_consistent
cp ../1.Geometry_optimization/{CONTCAR,CHGCAR,WAVECAR} .
cp CONTCAR POSCAR
cp ../POTCAR .
 cat > KPOINTS << EOF
K-Points
0
Monkhorst Pack
6 6 1               # please change it for other system
0 0 0
EOF

cat > INCAR << EOF
ISTART = 1                             # Continue from existing WAVECAR
ICHARG = 1                             # Read charge density from CHGCAR and extrapolate
ISMEAR = 0
NSW = 0                                # No ionic relaxation (static calculation)
IBRION = -1                            # No ionic updates
ISIF = 2                               # Relax ions and cell volume/shape
ENCUT = 425                            # Plane wave cutoff energy (eV)
NELM = 200                             # Maximum number of electronic steps
${BASE_INCAR}
EOF

gerun vasp_std

#===============================================================================
# 3. Density of States Calculation
#===============================================================================
cd .. && mkdir 3.DOS_calculation && cd 3.DOS_calculation
copy_scf_files
   
 cat > KPOINTS << EOF
K-Points
0
Monkhorst Pack
6 6 1               # please change it for other system
0 0 0
EOF

cat > INCAR << EOF
ICHARG = 11                           # Non-SCF calculation using existing charge density
ENCUT = 425                            # Plane wave cutoff energy (eV)
ISTART = 1                             # Continue from existing WAVECAR
NSW = 0                                # No ionic relaxation
IBRION = -1                            # No ionic updates
ISIF = 2                               # Relax ions and cell volume/shape
ISMEAR = -5
#NELM = 1                               # Single electronic step (non-SCF)
NEDOS = 3000                           # Number of grid points for DOS
LORBIT = 11                            # Calculate PDOS with phase factors
${BASE_INCAR}
EOF

gerun vasp_std
echo -e "11\n113" | vaspkit

#===============================================================================
# 4. Band Structure Calculation
#===============================================================================
cd .. && mkdir 4.Band_calculation && cd 4.Band_calculation
copy_scf_files

cat > KPOINTS << EOF
${KPATH}
EOF

cat > INCAR << EOF
ICHARG = 11                            # Non-SCF calculation using existing charge density
ISTART = 1                             # Continue from existing WAVECAR
NSW = 0                                # No ionic relaxation
ISIF = 2                               # Relax ions and cell volume/shape
ISMEAR  = -5    # Must use Gamma-centered KPOINTs
ENCUT = 425                            # Plane wave cutoff energy (eV)
IBRION = -1                            # No ionic updates
#NELM = 1                               # Single electronic step (non-SCF)
LORBIT = 11                            # Calculate PDOS with phase factors
${BASE_INCAR}
EOF

gerun vasp_std
echo -e "21\n211" | vaspkit
