package etomica.data;

import etomica.Atom;
import etomica.Data;
import etomica.DataInfo;
import etomica.DataSource;
import etomica.Space;
import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.data.types.DataVector;
import etomica.space.ICoordinateKinetic;
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
public class DataSourceVelocityAverage extends AtomActionAdapter implements DataSource {

    public DataSourceVelocityAverage(Space space) {
        this.space = space;
        data = new DataVector(space, new DataInfo("Average Velocity", Dimension.UNDEFINED));
        vectorSum = space.makeVector();
        groupWrapper = new AtomGroupAction(new MyAction());
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
        data.x.Ea1Tv1(1.0 / numAtoms, vectorSum);
        return data;
    }
    
    /**
     * Returns the center of mass (calculated so far) as a vector.
     * Does not reset COM sums.
     */
    public Vector getVelocityAverage() {
        getData();
        return data.x;
    }

    /**
     * Sets all accumulation sums to zero, readying 
     * for new calculation
     */
    public void reset() {
        vectorSum.E(0.0);
        numAtoms = 0;
    }

    private class MyAction extends AtomActionAdapter {
        public void actionPerformed(Atom a) {
            vectorSum.PE(((ICoordinateKinetic)a.coord).velocity());
            numAtoms++;
        }
    }

    private final Space space;
    private final Vector vectorSum;
    private final DataVector data;
    private double massSum;
    private int numAtoms;
    private AtomGroupAction groupWrapper;
}