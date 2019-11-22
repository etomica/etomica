package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.meter.MeterStructureFactor;

public class AtomSignalStress implements MeterStructureFactor.AtomSignalSource {

    protected final AtomStressSource stressSource;
    protected final int comp;

    public AtomSignalStress(AtomStressSource stressSource, int comp) {
        this.stressSource = stressSource;
        this.comp = comp;
    }

    public double signal(IAtom atom) {
        return stressSource.getStress()[atom.getLeafIndex()][comp];
    }
}
