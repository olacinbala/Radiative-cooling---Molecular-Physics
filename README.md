# Radiative cooling (Molecular Physics)


# About the simulation

This code simulates the spontaneous radiative cooling of vibrationally excited molecules in isolated conditions. 
The output files are
_ the emission spectra of the molecule resulting from this spontaneous radiative decay.
_ Temporal decay of internal (vibrational)energy of the molecule

This kind of photophysics is suspected to occur in interstellar and circumstellar environments for molecules containing dozens of carbon atoms.

# Utility of this code

This code has been the starting point for various simulations which leaded to three research publications related to the recurrent fluorescence
of carbon clusters in the interstellar medium context. 
These publications are listed below:
O. Lacinbala et al, J. Chem. Phys, 156, 144305 (2022) https://doi.org/10.1063/5.0080494
O. Lacinbala et al, Astron. Astrophys., 671, A86 (2023) https://doi.org/10.1051/0004-6361/202245421 
O. Lacinbala et al, Phys. Rev. A, 107, 062808 (2023) https://doi.org/10.1103/PhysRevA.107.062808

The simulation is performed by using a kinetic Monte Carlo method. The radiative cooling takes place via vibrational emission and recurrent fluorescence.

# Method

This code simulates the radiative cooling of highly vibrationally excited molecules in the harmonic approximation of vibrational normal modes.
This simulation requires several inputs about:
- vibrational structure (vibrational normal modes, IR intensities) of the molecule.
- electronic structure (electronic energies and electronic oscillator strength) of the molecule.

The vibrational density of states (VDOS) are computed with the Beyer-Swinehart algorithm. The VDOS allows to compute the
occupation probability of vibronic states. Then, the vibrational emission and recurrent fluorescence rates are computed for each internal energy (energy bin= 1 cm^-1).

Once all rates are computed at each internal energy and the initial vibrational energy of the molecule is defined, the kinetic Monte Carlo simulation is run to simulate the radiative cascade.

# Output files

Several output files are generated at the end of this simulation.

decroissance_E.txt: Temporal evolution of the internal energy

emissionIR_niveaux_elec.txt: time- and internal energy-dependence of total emitted photon from all vibrational modes for each electronic states.

emissionIR_par_modes.txt: time- and internal energy-dependence of emitted photon from each vibrational normal modes of the ground electronic states.

emissionIR_par_modes_global.txt: total emitted photon from each vibrational normal modes resulting from the radiative relaxation.

emissionUV.txt: time- and internal energy-dependence of total emitted photon via recurrent fluorescence from each electronic states.

spectreUV.txt: Total emitted photon number via recurrent fluorescence from each electronic states.
