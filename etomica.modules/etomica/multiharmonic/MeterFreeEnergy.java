package etomica.multiharmonic;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.potential.P1Harmonic;
import etomica.units.Dimension;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Oct 28, 2005 by kofke
 */
public class MeterFreeEnergy extends DataSourceScalar implements Meter {

    /**
     * @param label
     * @param dimension
     */
    public MeterFreeEnergy(P1Harmonic reference, P1Harmonic target) {
        super("Free energy", Dimension.ENERGY);
        this.reference = reference;
        this.target = target;
    }

    /* (non-Javadoc)
     * @see etomica.data.DataSourceScalar#getDataAsScalar()
     */
    public double getDataAsScalar() {
        iterator.reset();
        double sum = 0.0;
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            sum += target.energy(a) - reference.energy(a);
        }
        return Math.exp(-sum);
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.Meter#setPhase(etomica.phase.Phase)
     */
    public void setPhase(Phase phase) {
        iterator.setPhase(phase);
        this.phase = phase;
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.Meter#getPhase()
     */
    public Phase getPhase() {
        return phase;
    }

    AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    Phase phase;
    P1Harmonic reference, target;
}
