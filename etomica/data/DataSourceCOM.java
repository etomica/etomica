package etomica.data;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionDefinition;
import etomica.atom.AtomTypeLeaf;
import etomica.data.types.DataVector;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Length;

/**
 * Calculates the center of mass (COM) over a set of atoms. The mass and
 * position of all atoms passed to the actionPerformed are accumulated and used
 * to compute their center of mass. Calculated COM is obtained via the getData
 * method, which returns an array with the elements of the calculated COM
 * vector, or via the getCOM method, which returns a vector with the calculated
 * COM. Calculation is zeroed via the reset method.
 * <p>
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
        vectorSum = space.makeVector();
        data = new DataVector(space, "Center of Mass", Length.DIMENSION);
        myAction = new MyAction(vectorSum);
        groupWrapper = new AtomGroupAction(myAction);
    }
    
    public DataInfo getDataInfo() {
        return data.getDataInfo();
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
    public Data getData() {
        data.x.Ea1Tv1(1.0 / myAction.massSum, vectorSum);
        return data;
    }
   
    /**
     * Returns the center of mass (calculated so far) as a vector.
     * Does not reset COM sums.
     */
    public Vector getCOM() {
        getData();
        return data.x;
    }
    
    //AtomPositionDefinition interface implementation
    public Vector position(Atom atom) {
        reset();
        actionPerformed(atom);
        return getCOM();
    }

    /**
     * Sets all accumulation sums to zero, readying 
     * for new center-of-mass calculation
     */
    public void reset() {
        vectorSum.E(0.0);
        myAction.massSum = 0.0;
    }

    private static class MyAction extends AtomActionAdapter {
        public MyAction(Vector sum) {
            vectorSum = sum;
        }
        public void actionPerformed(Atom a) {
            double mass = ((AtomTypeLeaf)a.type).getMass();
            vectorSum.PEa1Tv1(mass, ((AtomLeaf)a).coord.position());
            massSum += mass;
        }
        private Vector vectorSum;
        public double massSum;
    }
    
    private final Vector vectorSum;
    private final DataVector data;
    private AtomGroupAction groupWrapper;
    private MyAction myAction;
}