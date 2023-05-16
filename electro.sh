#!/bin/bash
##################################### EDIT THE INFORMATION ##############
KPOINTS_unitcell="5 5 1" 
MAGMOM='2*2.0 2*2.0 2*2.0' 

incar="LWAVE = .TRUE.
LCHARG = .TRUE.
LREAL = Auto
EDIFF = 1E-6
EDIFFG = -0.01
ISMEAR = 0
SIGMA = 0.05
POTIM = 0.4
ENCUT = 425
PREC = Accurate
ALGO = Fast
NCORE = 8
ISIF = 2
LORBIT = 11
LMAXMIX = 4
ISPIN = 2"

KPATH="KPOINTS information
  40
Line-Mode
Reciprocal
   0.0000000000   0.0000000000   0.0000000000     GAMMA          
   0.5000000000   0.0000000000   0.0000000000     X              
 
   0.5000000000   0.0000000000   0.0000000000     X              
   0.5000000000   0.5000000000   0.0000000000     M              
 
   0.5000000000   0.5000000000   0.0000000000     M              
   0.00000000       0.500000000  0.00000000       Y  
   
   0.00000000       0.500000000   0.00000000      Y  
   0.0000000000   0.0000000000   0.0000000000    GAMMA"

###################################END OF EDITING##########################################


mkdir 1.relax_volume && cd 1.relax_volume && cp ../POSCAR . && cp ../POTCAR .
cat>KPOINTS <<**
K-Points
 0
Monkhorst Pack         
$KPOINTS_unitcell
 0  0  0    
**
cat>INCAR <<**
ICHARG = 2
ISTART = 0
ISMEAR = 0
NSW = 1500
IBRION = 2
ISIF = 4 
$incar
**
gerun vasp_std
gnuplot --persist ~/usr/check.g
#####################################################################
cd ../ && mkdir 2.relax_only_atom && cd 2.relax_only_atom && cp ../1.relax_volume/CONTCAR . && cp CONTCAR POSCAR && cp ../POTCAR .

cat>KPOINTS <<**
K-Points
 0
Monkhorst Pack          
$KPOINTS_unitcell
 0  0  0    
**
cat>INCAR <<**
ICHARG = 2
ISTART = 0
ISMEAR = 0
NSW = 1500
IBRION = 2
ISIF = 2  
$incar
**
gerun vasp_std
gnuplot --persist ~/usr/check.g
###############################3.Self_consistent##########
cd ../ && mkdir 3.Self_consistent && cd 3.Self_consistent && cp ../2.relax_only_atom/CONTCAR . && cp ../2.relax_only_atom/CHGCAR . && cp ../2.relax_only_atom/WAVECAR . && cp CONTCAR POSCAR && cp ../POTCAR .
cat>INCAR <<**
ICHARG = 1
LREAL = Auto
ISTART = 1
ISMEAR = 0
NSW = 0
IBRION = -1
NELM = 200 
$incar
**
cat>KPOINTS <<**
K-Points
 0
Monkhorst Pack          #change it
$KPOINTS_unitcell
 0  0  0    
**
gerun vasp_std
gnuplot --persist ~/usr/check.g
#####################################################N0n-Self_consistent
cd ../ && mkdir 3.DOS && cd 3.DOS && cp ../3.Self_consistent/CONTCAR . && cp ../3.Self_consistent/CHGCAR . && cp ../3.Self_consistent/WAVECAR . && cp CONTCAR POSCAR && cp ../POTCAR .

cat>INCAR <<**
ICHARG = 11
ISTART = 1
ISMEAR = 0
NSW = 0
IBRION = -1
ISIF = 2
NELM = 200  
NEDOS = 3000
$incar
**
gerun vasp_std
gnuplot --persist ~/usr/check.g

echo "11" | vaspkit
echo "113" | vaspkit

cd ../ && mkdir 4.Bands && cd 4.Bands && cp ../3.Self_consistent/CONTCAR . && cp ../3.Self_consistent/CHGCAR . && cp ../3.Self_consistent/WAVECAR . && cp CONTCAR POSCAR && cp ../POTCAR .

cat>INCAR <<**
ICHARG = 11
ISTART = 1
ISMEAR = 0
NSW = 0
IBRION = -1
ISIF = 2
LORBIT = 11
NEDOS = 3000
$incar
**
cat>KPOINTS <<**
$KPATH
**
gerun vasp_std


echo "21" | vaspkit
echo "211" | vaspkit
gnuplot --persist ~/usr/band.g
done

