package etomica.atom;

import etomica.Atom;
import etomica.DataSource;
import etomica.DataTranslator;
import etomica.Space;
import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.data.DataTranslatorVector;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Calculates the geometric center over a set of atoms. The position of all
 * atoms passed to the actionPerformed are accumulated and used to compute their
 * center (unweighted by mass). Calculated center is obtained via the getData
 * method, which returns an array with the elements of the calculated center
 * vector, or via the getCenter method, which returns a vector with the
 * calculated center. Calculation is zeroed via the reset method. <br>
 * A typical use of this class would have it passed to the allAtoms method of an
 * iterator, or wrapped in an AtomGroupAction to calculate the geometric center
 * of the atoms in an atom group.
 * 
 * @author David Kofke
 */

/*
 * History Created on Jan 26, 2005 by kofke
 */
public class AtomPositionGeometricCenter extends AtomActionAdapter implements DataSource,
        AtomPositionDefinition {

    public AtomPositionGeometricCenter(Space space) {
        this.space = space;
        vectorSum = space.makeVector();
        center = space.makeVector();
        dataTranslator = new DataTranslatorVector(space);
        groupWrapper = new AtomGroupAction(new MyAction());
        setLabel("Geometric center");
    }

    /*
     * (non-Javadoc)
     * 
     * @see etomica.action.AtomAction#actionPerformed(etomica.Atom)
     */
    public void actionPerformed(Atom a) {
        groupWrapper.actionPerformed(a);
    }

    /**
     * Returns the center of mass (calculated so far) as an array. Does not
     * reset COM sums. Implementation of DataSource interface.
     */
    public double[] getData() {
        return dataTranslator.toArray(getCenter());
    }

    /**
     * Returns the lenght of the data returned, which equals the dimension of
     * the space.
     */
    public int getDataLength() {
        return space.D();
    }

    /**
     * Returns the center of mass (calculated so far) as a vector. Does not
     * reset center-position sums.
     */
    public Vector getCenter() {
        center.Ea1Tv1(1.0 / natoms, vectorSum);
        return center;
    }

    //AtomPositionDefinition interface implementation
    public Vector position(Atom atom) {
        reset();
        actionPerformed(atom);
        return getCenter();
    }

    /**
     * Returns Dimension.LENGTH
     */
    public Dimension getDimension() {
        return Dimension.LENGTH;
    }

    /**
     * Returns a new DataTranslatorVector instance
     */
    public DataTranslator getTranslator() {
        return new DataTranslatorVector(space);
    }

    /**
     * Sets all accumulation sums to zero, readying for new center-of-mass
     * calculation
     */
    public void reset() {
        vectorSum.E(0.0);
        natoms = 0;
    }

    private class MyAction extends AtomActionAdapter {

        public void actionPerformed(Atom a) {
            vectorSum.PE(a.coord.position());
            natoms++;
        }
    }

    private final Space space;
    private final Vector vectorSum, center;
    private final DataTranslatorVector dataTranslator;
    private AtomGroupAction groupWrapper;
    private int natoms;
}