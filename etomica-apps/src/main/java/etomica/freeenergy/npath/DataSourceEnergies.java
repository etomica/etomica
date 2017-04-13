package etomica.freeenergy.npath;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialMaster;
import etomica.units.Energy;

/**
 * Created by andrew on 4/12/17.
 */
public class DataSourceEnergies implements IEtomicaDataSource {

    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final PotentialCalculationEnergies pc;
    protected final PotentialMaster potentialMaster;
    protected IBox box;
    protected final IteratorDirective id;

    public DataSourceEnergies(PotentialMaster potentialMaster) {
        this.potentialMaster = potentialMaster;
        data = new DataDoubleArray(2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("energies!", Energy.DIMENSION, new int[]{2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        pc = new PotentialCalculationEnergies();
        id = new IteratorDirective();
    }

    public void setBox(IBox box) {
        this.box = box;
    }

    @Override
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        double[] x = data.getData();
        x[0] = pc.getSum1();
        x[1] = pc.getSum2();
        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public static class PotentialCalculationEnergies implements PotentialCalculation {

        protected double sum1, sum2;

        public double getSum1() {
            return sum1;
        }

        public double getSum2() {
            return sum2;
        }

        public void reset() {
            sum1 = sum2 = 0;
        }

        @Override
        public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
            double u = potential.energy(atoms);
            if (potential.nBody() == 1) {
                sum1 += u;
            }
            else {
                sum2 += u;
            }
        }
    }
}
