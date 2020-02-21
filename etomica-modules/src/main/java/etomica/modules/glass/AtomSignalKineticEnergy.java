package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.data.meter.MeterStructureFactor;

public class AtomSignalKineticEnergy implements MeterStructureFactor.AtomSignalSource {

    protected double avg;

    public void setDoSubtractAvg(double avg) {
        this.avg = avg;
    }

    public double signal(IAtom atom) {
        double m = atom.getType().getIndex();
        if (m == 0) return 0;
        double KE = 0.5 * m * ((IAtomKinetic) atom).getVelocity().squared();
        return KE - avg;
    }
}
