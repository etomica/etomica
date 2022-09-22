package etomica.normalmode;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.units.dimensions.Null;

public class MeterPIPrim implements IDataSource {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double beta, betaN;
    protected int nBeads;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected double EnShift;

    public MeterPIPrim(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, int nBeads, double betaN) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        this.beta = this.betaN*nBeads;
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
        pcP1.computeAll(false);
        x[0] = 1.0/2.0/betaN + pcP1.getLastEnergy() - pmBonding.getLastEnergy() - EnShift; //En
        x[1] = nBeads/2.0/beta/beta - 2*pmBonding.getLastEnergy()/beta;  //Cvn/kb^2, without Var
        return data;
    }

    public void setShift(double EnShift){
        this.EnShift = EnShift;
    }

}
