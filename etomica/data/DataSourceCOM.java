package etomica.data;

import etomica.Atom;
import etomica.DataSource;
import etomica.DataTranslator;
import etomica.Space;
import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.atom.AtomPositionDefinition;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * Calculates the center of mass (COM) over a set of atoms. The mass and
 * position of all atoms passed to the actionPerformed are accumulated and used
 * to compute their center of mass. Calculated COM is obtained via the getData
 * method, which returns an array with the elements of the calculated COM
 * vector, or via the getCOM method, which returns a vector with the calculated
 * COM. Calculation is zeroed via the reset method. <br>
 * A typical use of this class would have it passed to the allAtoms method of an
 * iterator, or wrapped in an AtomGroupAction to calculate the COM of the atoms
 * in an atom group.
 * 
 * @author David Kofke
 */

/*
 * History Created on Jan 26, 2005 by kofke
 */
public class DataSourceCOM extends AtomActionAdapter implements DataSource, AtomPositionDefinition {

    public DataSourceCOM(Space space) {
        this.space = space;
        vectorSum = space.makeVector();
        com = space.makeVector();
        dataTranslator = new DataTranslatorVector(space);
        groupWrapper = new AtomGroupAction(new MyAction());
        setLabel("COM");
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
     * Returns the center of mass (calculated so far) as an array.
     * Does not reset COM sums.  Implementation of DataSource interface.
     */
    public double[] getData() {
        return dataTranslator.toArray(getCOM());
    }

    /**
     * Returns the center of mass (calculated so far) as a vector.
     * Does not reset COM sums.
     */
    public Vector getCOM() {
        com.Ea1Tv1(1.0 / massSum, vectorSum);
        return com;
    }
    
    //AtomPositionDefinition interface implementation
    public Vector position(Atom atom) {
        reset();
        actionPerformed(atom);
        return getCOM();
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
     * Sets all accumulation sums to zero, readying 
     * for new center-of-mass calculation
     */
    public void reset() {
        vectorSum.E(0.0);
        massSum = 0.0;
    }

    private class MyAction extends AtomActionAdapter {
        public void actionPerformed(Atom a) {
            double mass = a.type.getMass();
            vectorSum.PEa1Tv1(mass, a.coord.position());
            massSum += mass;
        }
    }
    
    private final Space space;
    private final Vector vectorSum, com;
    private final DataTranslatorVector dataTranslator;
    private double massSum;
    private AtomGroupAction groupWrapper;
}