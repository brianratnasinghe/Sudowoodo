#!/bin/bash
set -euo pipefail

# Energy minimization
gmx grompp -f EM.mdp -c afm_system.gro -p afm_system.top -o EM.tpr
gmx mdrun -deffnm EM -ddcheck -ntmpi 1 -ntomp 24 -dlb no

# Equilibration
gmx grompp -f EQ.mdp -c EM.gro -p afm_system.top -o EQ.tpr -maxwarn 2
gmx mdrun -deffnm EQ -ddcheck -ntmpi 1 -ntomp 24 -dlb no -v

# Production
gmx grompp -f production.mdp -c EQ.gro -p afm_system.top -o production.tpr
gmx mdrun -deffnm production -ddcheck -ntmpi 1 -ntomp 24 -dlb no -v