package etomica.modules.entropylottery;

import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IVector;

import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.units.Null;

/**
 * Calculates the entropy of the distribution of Atoms in a Box
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
        IVector dimensions = box.getBoundary().getDimensions();
        if (atomCount.length != (int)Math.round(dimensions.x(0))) {
            atomCount = new int[(int)Math.round(dimensions.x(0))];
        }
        for (int i=0; i<atomCount.length; i++) {
            atomCount[i] = 0;
        }
        atomIterator.setBox(box);
        atomIterator.reset();
        for (IAtomPositioned a = (IAtomPositioned)atomIterator.nextAtom(); a != null;
             a = (IAtomPositioned)atomIterator.nextAtom()) {
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

    public IBox getBox() {
        return box;
    }

    public void setBox(IBox newBox) {
        box = newBox;

    }

    private static final long serialVersionUID = 1L;
    protected IBox box;
    protected int[] atomCount;
    protected final AtomIteratorLeafAtoms atomIterator;
}
