package etomica.modules.entropylottery;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.units.Null;

/**
 * Calculates the entropy of the distribution of Atoms in a Phase
 * using Stirling's approximation.
 * @author Andrew Schultz
 */
public class MeterEntropy extends DataSourceScalar implements DataSource {

    public MeterEntropy() {
        super("entropy", Null.DIMENSION);
        atomCount = new int[0];
        atomIterator = new AtomIteratorLeafAtoms();
    }

    public double getDataAsScalar() {
        IVector dimensions = phase.getBoundary().getDimensions();
        if (atomCount.length != (int)Math.round(dimensions.x(0))) {
            atomCount = new int[(int)Math.round(dimensions.x(0))];
        }
        for (int i=0; i<atomCount.length; i++) {
            atomCount[i] = 0;
        }
        atomIterator.setPhase(phase);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            AtomLeaf a = (AtomLeaf)atomIterator.nextAtom();
            int x = (int)Math.round(a.getPosition().x(0)+dimensions.x(0)*0.5-0.5);
            atomCount[x]++;
        }
        
        double sum = 0;
        double atomTotal = atomIterator.size();
        for (int i=0; i<atomCount.length; i++) {
            if (atomCount[i] > 0) {
                sum += atomCount[i] * Math.log(atomCount[i]/atomTotal);
            }
        }
        return -sum;
    }

    public Phase getPhase() {
        return phase;
    }

    public void setPhase(Phase newPhase) {
        phase = newPhase;

    }

    private static final long serialVersionUID = 1L;
    protected Phase phase;
    protected int[] atomCount;
    protected final AtomIteratorLeafAtoms atomIterator;
}
