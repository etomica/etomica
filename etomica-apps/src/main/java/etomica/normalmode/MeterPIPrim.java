package etomica.normalmode;

import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Null;

public class MeterPIPrim implements IDataSource {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected int nBeads;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected double beta;
    protected int numAtoms;
    protected Box box;
    protected double EnShift;
    protected double volume;

    public MeterPIPrim(PotentialMasterBonding pmBonding, PotentialCompute pcP1, int nBeads, double temperature, Box box) {
        int nData = 2;
        this.volume = box.getBoundary().volume();
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        this.box = box;
        this.EnShift = 0;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }


    @Override
    public IData getData() {
        double[] x = data.getData();

        pmBonding.computeAll(false);
        pcP1.computeAll(true);
        x[0] = dim*numAtoms*nBeads/2.0/beta + pcP1.getLastEnergy() - pmBonding.getLastEnergy() - EnShift; //En
        x[1] = dim*numAtoms*nBeads/2.0/beta/beta - 2*pmBonding.getLastEnergy()/beta + x[0]*x[0];
//        x[2] = nBeads*numAtoms/volume/beta - 2.0/dim/volume*pmBonding.getLastEnergy() - pcP1.getLastVirial()/dim/volume;
        return data;
    }

    public void setEnShift(double E) { this.EnShift = E; }
}