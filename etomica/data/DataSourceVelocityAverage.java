package etomica.data;

import etomica.Atom;
import etomica.AtomTypeLeaf;
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
 * Calculates the mass average velocity over a set of atoms. The velocity
 * and mass of all atoms passed to the actionPerformed are accumulated and used
 * to compute the mass-average velocity (total momentum divided by the total mass).
 * Calculated average is obtained via the getData
 * method, which returns a DataVector with the average.
 * Calculation is zeroed via the reset method. <br>
 * 
 * @author David Kofke
 */

/*
 * History Created on Jan 26, 2005 by kofke
 */
public class DataSourceVelocityAverage extends AtomActionAdapter implements DataSource {

    public DataSourceVelocityAverage(Space space) {
        data = new DataVector(space, "Average Velocity", Dimension.UNDEFINED);
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
     * Does not reset sums.  Implementation of DataSource interface.
     */
    public Data getData() {
        data.x.Ea1Tv1(1.0 / massSum, vectorSum);
        return data;
    }
    
    /**
     * Returns the average velocity (calculated so far) as a vector.
     * Does not reset sums.
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
        massSum = 0;
    }

    private class MyAction extends AtomActionAdapter {
        public void actionPerformed(Atom a) {
            vectorSum.PE(((ICoordinateKinetic)a.coord).velocity());
            massSum += ((AtomTypeLeaf)a.type).getMass();
        }
    }

    private final Vector vectorSum;
    private final DataVector data;
    private double massSum;
    private AtomGroupAction groupWrapper;
}