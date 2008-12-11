package etomica.api;


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
}