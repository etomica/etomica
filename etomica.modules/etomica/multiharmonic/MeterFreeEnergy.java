package etomica.multiharmonic;

import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.potential.P1Harmonic;
import etomica.units.Dimension;


/**
 * Computes the free-enery difference through free-energy perturbation.
 * Written specifically for harmonic 1-body potential, but wouldn't be hard
 * to modify for more general cases.
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
//        if(sum < 0) System.out.println("deltaU: "+sum);
//       System.out.println(target.getSpringConstant()+" "+target.getX0().x(0)+" "+reference.getSpringConstant()+" "+reference.getX0().x(0));
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
