package etomica.modules.multiharmonic;

import etomica.api.IBox;
import etomica.atom.AtomSet;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.potential.P1Harmonic;
import etomica.units.Energy;


/**
 * Computes the free-enery difference through free-energy perturbation.
 * Written specifically for harmonic 1-body potential, but wouldn't be hard
 * to modify for more general cases.
 *
 * @author David Kofke
 *
 */
public class MeterFreeEnergy extends DataSourceScalar implements DataSource {
    
    public MeterFreeEnergy(P1Harmonic reference, P1Harmonic target) {
        super("Free energy", Energy.DIMENSION);
        this.reference = reference;
        this.target = target;
    }
    
    public double getDataAsScalar() {
        iterator.reset();
        double sum = 0.0;
        for (AtomSet a = iterator.next(); a != null; a = iterator.next()) {
            sum += target.energy(a) - reference.energy(a);
        }
        return Math.exp(-sum);
    }
    
    public void setBox(IBox box) {
        iterator.setBox(box);
        this.box = box;
    }
    
    public IBox getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    IBox box;
    P1Harmonic reference, target;
}
