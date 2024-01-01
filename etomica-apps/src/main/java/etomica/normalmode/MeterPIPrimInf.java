package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Null;

public class MeterPIPrimInf implements IDataSource {
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
    protected double omega, R, cothR, hbar;

    public MeterPIPrimInf(PotentialMasterBonding pmBonding, PotentialCompute pcP1, int nBeads, double temperature, Box box, double omega, double hbar) {
        int nData = 2;
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
        this.omega = omega;
        this.hbar = hbar;
        this.R = beta*hbar*omega/nBeads;
        this.cothR = Math.cosh(R)/Math.sinh(R);
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
        pcP1.computeAll(false);

        if (numAtoms==1){
            x[0] = dim*hbar*omega/2*cothR - R*cothR*pmBonding.getLastEnergy() + pcP1.getLastEnergy();
        } else {
            x[0] = dim*(numAtoms*nBeads-1)/2.0/beta + pcP1.getLastEnergy() - pmBonding.getLastEnergy(); //En
        }
        x[1] = nBeads/2.0/beta/beta - 2*pmBonding.getLastEnergy()/beta + x[0]*x[0];
        return data;
    }
}

