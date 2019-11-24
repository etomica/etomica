package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.meter.MeterStructureFactor;

public class AtomSignalStress implements MeterStructureFactor.AtomSignalSource {

    protected final AtomStressSource stressSource;
    protected final int comp1, comp2;

    public AtomSignalStress(AtomStressSource stressSource, int comp1, int comp2) {
        this.stressSource = stressSource;
        this.comp1 = comp1;
        this.comp2 = comp2;
    }

    public double signal(IAtom atom) {
        return stressSource.getStress()[atom.getLeafIndex()][comp1][comp2];
    }
}
