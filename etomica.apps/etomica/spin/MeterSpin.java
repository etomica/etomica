package etomica.spin;

import etomica.Atom;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.space.Vector;
import etomica.units.Dimension;


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
        super("Spin",Dimension.UNDEFINED);
        sum = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.MeterScalar#getDataAsScalar(etomica.Phase)
     */
    public double getDataAsScalar() {
        sum.E(0.0);
        int count = 0;
        iterator.setList(phase.getSpeciesMaster().atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
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
    private final AtomIteratorListTabbed iterator = new AtomIteratorListTabbed();
    private final Vector sum;
}
