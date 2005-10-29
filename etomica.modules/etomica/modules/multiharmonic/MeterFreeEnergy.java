package etomica.modules.multiharmonic;

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
    
    public MeterFreeEnergy(P1Harmonic reference, P1Harmonic target) {
        super("Free energy", Dimension.ENERGY);
        this.reference = reference;
        this.target = target;
    }
    
    public double getDataAsScalar() {
        iterator.reset();
        double sum = 0.0;
        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            sum += target.energy(a) - reference.energy(a);
        }
        return Math.exp(-sum);
    }
    
    public void setPhase(Phase phase) {
        iterator.setPhase(phase);
        this.phase = phase;
    }
    
    public Phase getPhase() {
        return phase;
    }

    AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    Phase phase;
    P1Harmonic reference, target;
}
