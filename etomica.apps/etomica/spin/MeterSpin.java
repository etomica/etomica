package etomica.spin;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Undefined;


/**
 * Meter that provides the x-component of the vector average of
 * spin values (which is represented by the atom's position
 * vector).
 *
 * @author David Kofke
 *
 */
public class MeterSpin extends DataSourceScalar implements DataSource {

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
            sum.PE(atom.getPosition());
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

    private static final long serialVersionUID = 1L;
    private Phase phase;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private final IVector sum;
}
