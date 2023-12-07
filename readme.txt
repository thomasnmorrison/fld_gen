readme.txt

Readme file for fld_gen.
This is a program to generate fields using a combination of (Fourier space) transfer functions and pointwise nonlinear functions. There is not currently support for general nonlinear functions of the field.

Files
------------------
fld_gen.f90          - Top program
io_mod.f90
lin_tran.f90
nlin_tran.f90
spec_init.f90
fftw_mod.f90

eom_pert_cosmic.f90  - Equations of motion for perturbations
eom_bg_cosmic.f90    - Equations of motion for background
pert_macros.h        - macros

newton_root.f90      - Root finding for GL integrator
gl_integrator.f90    - GL integrator
butcher_table.f90
gl_macros.h

params.f90           - Has parameters for the lattice
macros.h
vars.f90             - Field variables and initializers

Makefile

Outputs
------------------
gfld*.out            - Binary file of initial Gaussian field
ngfld*.out           - Binary file after initial Gaussian field has gone through a nonlinear transform
ngtfld*.out          - Binary file after nonGaussian field has gone through a linear transform
