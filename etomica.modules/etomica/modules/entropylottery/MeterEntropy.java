package etomica.modules.entropylottery;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.data.DataSourceScalar;
import etomica.data.IEtomicaDataSource;
import etomica.units.Null;

/**
 * Calculates the entropy of the distribution of Atoms in a Box
 * using Stirling's approximation.
 * @author Andrew Schultz
 */
public class MeterEntropy extends DataSourceScalar implements IEtomicaDataSource {

    public MeterEntropy() {
        super("entropy", Null.DIMENSION);
        atomCount = new int[0];
    }

    public double getDataAsScalar() {
        IVector dimensions = box.getBoundary().getDimensions();
        if (atomCount.length != (int)Math.round(dimensions.getX(0))) {
            atomCount = new int[(int)Math.round(dimensions.getX(0))];
        }
        for (int i=0; i<atomCount.length; i++) {
            atomCount[i] = 0;
        }
        IAtomList leafList = box.getLeafList();
        for (int i=0; i<leafList.getAtomCount(); i++) {
            IAtomPositioned a = (IAtomPositioned)leafList.getAtom(i);
            int x = (int)Math.round(a.getPosition().getX(0)+dimensions.getX(0)*0.5-0.5);
            atomCount[x]++;
        }
        
        double sum = 0;
        double atomTotal = box.getLeafList().getAtomCount();
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
}
