"""
Atom object.

An atom object holds the following information: chemical symbol,
position, atomic number, mass, momentum, velocity, magnetic moment and
a tag.  The chemical symbol and the atomic number are kept syncronized
- the same goes for the momentum and the velocity."""

import numpy as np

from Constants import amu2me, me2amu
from Elements import symbols, masses, names, numbers


class Atom:
    """Atom object."""
    def __init__(self, symbol=None, position=[0, 0, 0],
                 Z=None, mass=None,
                 momentum=None, magmom=0.0, tag=0):
        """Atom(symbol, position, ...) -> atom object.

        The chemical symbol or the atomic number must be given
        (``symbol`` or ``Z``).  The rest of the arguments (``mass``,
        ``tag``, ``momentum``, ``velocity`` and ``magmom``) have
        default values."""
        
        if symbol is None:
            if Z is None:
                raise ValueError, 'Missing symbol or atomic number!'
            symbol = symbols[Z]
        else:
            symbol = symbol.lower().capitalize() 
            if Z is None:
                Z = numbers[symbol]
            else:
                if symbols[Z] != symbol:
                    raise ValueError, 'Incompatible atomic number and symbol'

        self.Z = Z
        self.symbol = symbol
        self.name = names[Z]
        self.force = np.array([0,0,0],dtype="float") 
        self.SetCartesianPosition(position)
        if mass is None:
            mass = masses[Z] * amu2me
        self.SetMass(mass)
        if momentum is None:
	    self.SetCartesianMomentum([0, 0, 0])
        else:
	    self.SetCartesianMomentum(momentum)
        
        self.SetMagneticMoment(magmom)
        self.SetTag(tag)

    def __repr__(self):
        return "%s %s" % (self.symbol, self.position[:])

    def ndims(self):
        return len(self.position)

    def repr_trj(self):
        return "%s  % 10.8f  % 10.8f  % 10.8f" % (self.symbol,
                self.position[0],self.position[1],self.position[2])

    def __cmp__(self, other):
        return cmp(self.Z, other.Z)

    def Copy(self):
        return Atom(position=self.position, Z=self.Z, mass=self.mass,
                    tag=self.tag, momentum=self.momentum,
                    magmom=self.magmom)

    def SetCartesianPosition(self, position):
        """Set the position."""
        self.position = np.array(position, dtype=float)

    def GetCartesianPosition(self):
        """Get the position as a numarray."""
        return self.position

    def SetCartesianMomentum(self, momentum):
        """Set the momentum."""
        self.momentum = np.array(momentum, dtype=float)

    def GetCartesianMomentum(self):
        """Get the momentum as a numarray."""
        return self.momentum

    def SetCartesianVelocity(self, velocity):
        """Set the momentum from the velocity as a numarray."""
        self.momentum = np.array(velocity, dtype=float) * self.mass

    def GetCartesianVelocity(self):
        """Get the velocity as a numarray."""
        return self.momentum / self.mass

    def SetCartesianForce(self, force):
        """Set the force."""
        self.force = np.array(force, dtype=float)

    def GetCartesianForce(self):
        """Get the force as a numarray."""
        return self.force

    def SetCartesianAcceleration(self, acceleration):
        """Set the Acceleration."""
        self.force = np.array(acceleration, dtype=float) * self.mass 
    
    def GetCartesianAcceleration(self):
        """Get the acceleration as a numarray."""
        return self.force / self.mass

    def SetMass(self, mass):
        """Set the atoms mass."""
        self.mass = mass

    def GetAtomicMass(self):
        """Get the atomic mass."""
        return masses[self.Z]

    def GetMass(self):
        """Get the atomic mass."""
        return self.mass

    def SetMagneticMoment(self, magmom):
        """Set the atoms magnetic moment."""
        self.magmom = float(magmom)

    def GetMagneticMoment(self):
        """Get the atomic magnetic moment."""
        return self.magmom

    def SetTag(self, tag):
        """Set integer tag."""
        self.tag = int(tag)

    def GetTag(self):
        """Get integer tag."""
        return self.tag

    def SetAtomicNumber(self, Z):
        """Set the atomic number."""
        self.symbol = symbols[Z]
        self.Z = Z

    def GetAtomicNumber(self):
        """Get the atomic number."""
        return self.Z

    def SetChemicalSymbol(self, symbol):
        """Set the chemical symbol."""
        self.Z = numbers[symbol]
        self.symbol = symbol

    def GetChemicalSymbol(self):
        """Get the chemical symbol."""
        return self.symbol

    def GetElementName(self):
        """Get the element name."""
        return self.name

    def GetKineticEnergy(self):
        """Get the kinetic energy."""
        return 0.5 * np.dot(self.momentum, self.momentum) / self.mass

def test():
#    at1 = Atom(symbol='Ne',position=[0,0,0])
    at1 = Atom(symbol='Ne', position=[0,0,0])
    print at1
    print at1.GetAtomicNumber()
    print at1.GetMass()
    print at1.GetAtomicMass()
    print at1.GetElementName()
    at1.SetCartesianVelocity([2,2,2])
    print at1.GetCartesianMomentum()
    print at1.GetKineticEnergy()

def test2():
    at = Atom(symbol='D', position=[0,0,0])
    print at
    print at.GetAtomicNumber()
    print at.GetMass()
    print at.GetAtomicMass()
    print at.GetElementName()

def test3():
    at = Atom(symbol='Ge', position=[0,0,0])
    print at
    print at.GetAtomicNumber()
    print at.GetMass()
    print at.GetAtomicMass()
    print at.GetElementName()

if __name__ == '__main__': test()
