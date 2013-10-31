#!/usr/bin/python
"""
 Constants.py: Constants for simulation

"""

# Misc units
Clight    = 2.99792458e8   # speed of light in m/s
Clight_au = 1.3703574503876348e2    # in au
Kboltz    = 3.16683136e-6  # Boltzmann constant
e2 = 14.399             # Coulomb's law coeff if R in \AA and resulting E in eV
planck = 6.6260755e-34  # Planck's constant, in Js
Planck    = 6.6260755e-34           # Planck's constant, in Js

pi = 3.1415926535897931e0

# Distance units
bohr2ang = 0.529177249          # Conversion of length from bohr to angstrom
ang2bohr = 1.0/bohr2ang


# Energy units
kelvin2hartree = 3.16683136e-6        # Kelvin to Hartree
hartree2kelvin = 1.0/kelvin2hartree   # Hartree to Kelvin

cm2hartree = 4.55635e-6               # wave number to Hartree
hartree2cm = 1.0/cm2hartree           # Hartree to wave number

kj2hartree = 3.80880e-4               # kJ/mol to Hartree
hartree2kj = 1.0/kj2hartree           # Hartree to kJ/mol

kcal2hartree = 1.59360e-3             # kcal/mol to Hartree
hartree2kcal = 1.0/kcal2hartree       # Hartree to kcal/mol

ev2hartree = 3.67493e-2               # eV to Hartree
hartree2ev = 1.0/ev2hartree           # Hartree to eV

# Energy units part 2

kelvin2cm = 0.695039                # Kelvin to wave number
cm2kelvin = 1.0/kelvin2cm

kelvin2kj = 8.31451e-3              # Kelvin to kJ/mol
kj2kelvin = 1.0/kelvin2kj           

kelvin2kcal = 1.98722e-3            # Kelvin to kcal/mol
kcal2kelvin = 1.0/kelvin2kcal

kelvin2ev = 8.61739e-5              # Kelvin to eV
ev2kelvin = 1.0/kelvin2ev

cm2kj = 1.19627e-2                  # wave number to kJ/mol
kj2cm = 1.0/cm2kj

cm2kcal = 2.85914e-3                # wave number to kcal/mol
kcal2cm = 1.0/cm2kcal

cm2ev = 1.23985e-4                  # wave number to eV
ev2cm = 1.0/cm2ev

kj2kcal = 2.39006e-1                # kJ/mol to kcal/mol
kcal2kj = 1.0/kj2kcal

kj2ev = 1.03643e-2                  # kJ/mol to eV
ev2kj = 1.0/kj2ev

kcal2ev = 4.33642e-2                # kcal/mol to eV
ev2kcal = 1.0/kcal2ev

hartree2joule = 4.3597482e-18       # Hartree to Joule
joule2hartree = 1/hartree2joule     # Joule to Hartree


# Mass units
amu2me = 1822.888                   # mass in amu to mass in au (m_e)
me2amu = 1/amu2me


# Time units
ps2tau = 41341.447                  # time in ps to time in au
tau2ps = 1/ps2tau

fs2tau = ps2tau/1000.0
tau2fs = 1.0/fs2tau 

# Derived quantities
Rgas = Kboltz*hartree2kcal*1000.0   # gas constant R = 1.98722 cal/mole/K

# Frequency unit
cm2au = 2.e0*pi*Clight*1.e2*1.e-12/ps2tau  #  cm-1 (frequency) to au
au2cm = 1.e0/cm2au

debye2au = 0.393456e0
