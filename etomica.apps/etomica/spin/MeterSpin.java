package etomica.spin;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Undefined;


/**
 * Meter that provides the x-component of the vector average of
 * spin values (which is represented by the atom's position
 * vector).
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 25, 2005 by kofke
 */
public class MeterSpin extends DataSourceScalar implements Meter {

    /**
     * 
     */
    public MeterSpin(Space space) {
        super("Spin",Undefined.DIMENSION);
        sum = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.MeterScalar#getDataAsScalar(etomica.Phase)
     */
    public double getDataAsScalar() {
        sum.E(0.0);
        int count = 0;
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
            sum.PE(atom.coord.position());
            count++;
        }
        return sum.x(0)/count;
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private final Vector sum;
}
