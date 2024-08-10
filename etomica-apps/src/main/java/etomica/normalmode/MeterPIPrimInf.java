package etomica.normalmode;

import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.units.dimensions.Null;

public class MeterPIPrimInf implements IDataSource {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected int nBeads;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected double beta;
    protected int numAtoms;
    protected Box box;
    protected double hbar, omega, R, sinhR, tanhR, cothR;
    protected double fac1, fac2, fac3;
    public MeterPIPrimInf(PotentialMasterBonding pmBonding, PotentialCompute pcP1harm,PotentialCompute pcP1ah, int nBeads, double temperature, Box box, double omega, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1harm = pcP1harm;
        this.pcP1ah = pcP1ah;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        this.omega = omega;
        this.hbar = hbar;
        this.R = beta*hbar*omega/nBeads;
        this.sinhR = Math.sinh(R);
        this.tanhR = Math.tanh(R);
        this.cothR = 1/tanhR;
        this.fac1 = dim*numAtoms*nBeads/2.0/beta/beta*R*R/sinhR/sinhR;
        this.fac2 = R*R/beta*(1-2*cothR*cothR);
        this.fac3 = R*R/2/beta/Math.cosh(R/2)/Math.cosh(R/2);

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
        pcP1harm.computeAll(false);
        pcP1ah.computeAll(false);

//        x[0] = dim*numAtoms*hbar*omega/2*cothR - R*cothR*pmBonding.getLastEnergy() + R/sinhR*pcP1harm.getLastEnergy()+pcP1ah.getLastEnergy()
//        - dim*numAtoms*hbar*omega/2/Math.tanh(beta*hbar*omega/2);
        x[0] = dim*numAtoms*hbar*omega/2*cothR - R*cothR*pmBonding.getLastEnergy() + R/sinhR*pcP1harm.getLastEnergy()+pcP1ah.getLastEnergy();
        x[1] = fac1 + fac2*pmBonding.getLastEnergy() + fac3*pcP1harm.getLastEnergy() + x[0]*x[0];
        return data;
    }
}