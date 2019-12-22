package etomica.modules.glass;

import etomica.atom.IAtom;
import etomica.data.meter.MeterStructureFactor;

public class AtomSignalStress implements MeterStructureFactor.AtomSignalSource {

    protected final AtomStressSource stressSource;
    protected final int[][] comps;
    protected double avg;

    public AtomSignalStress(AtomStressSource stressSource, int[][] comps) {
        this.stressSource = stressSource;
        this.comps = comps;
    }

    public boolean ready() {
        double[][][] stress = stressSource.getStress();
        avg = 0;
        for (int i=0; i<stress.length; i++) {
            double sum = 0;
            for (int j=0; j<comps.length; j++) {
                sum += stress[i][comps[j][0]][comps[j][1]];
            }
            avg += (sum - avg)/(i+1);
        }
        return true;
    }

    public double signal(IAtom atom) {
        double sum = 0;
        double[][] stress = stressSource.getStress()[atom.getLeafIndex()];
        for (int i=0; i<comps.length; i++) {
            sum += stress[comps[i][0]][comps[i][1]];
        }
        return sum - avg;
    }
}
