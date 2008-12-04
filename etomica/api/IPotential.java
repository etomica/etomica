package etomica.api;

import etomica.units.Dimension;

public interface IPotential {

    /**
     * Returns the range over which the potential applies.  IAtoms with a
     * greater separation do not interact.
     */
    public double getRange();

    /**
     * Informs the potential of the box on which it acts so that it can
     * properly consider the boundaries.
     */
    public void setBox(IBox box);

    /**
     * The number of atoms on which the potential depends.
     */
    public int nBody();

    /**
     * Returns the interaction energy between the given atoms.  There might be
     * 0, 1, 2 or more atoms in the AtomSet.
     */
    public abstract double energy(IAtomList atoms);

}