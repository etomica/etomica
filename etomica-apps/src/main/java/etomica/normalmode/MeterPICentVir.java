package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPICentVir implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialComputeField pcP1;
    protected double betaN, beta;
    protected int nBeads;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector rc;
    protected double EnShift;


    public MeterPICentVir(PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        beta = this.betaN*this.nBeads;
        this.box = box;
        rc = box.getSpace().makeVector();
        this.EnShift = 0;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        double vir = 0;
        rc.E(0);
        Vector ri;
        for (int i = 0; i < nBeads; i++){
            ri = box.getLeafList().get(i).getPosition();
            rc.PE(ri);
        }
        rc.TE(1.0/nBeads);

        pcP1.computeAll(true, this);//it needs rc
        Vector[] forces = pcP1.getForces();
        for (int i = 0; i < nBeads; i++){
            ri = box.getLeafList().get(i).getPosition();
            vir -= forces[i].dot(ri);
            vir += forces[i].dot(rc);
        }
        x[0] = 1.0/2.0/beta + pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift; //En
        x[1] = 1.0/2.0/beta/beta + 1.0/4.0/beta*(-3.0*vir - rHr); //Cvn/kb^2, without Var
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
        Vector tmpV = box.getSpace().makeVector();
        tmpV.Ev1Mv2(rj, rc);
        Hij.transform(tmpV);
        rHr += ri.dot(tmpV);
        rHr -= rc.dot(tmpV);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public boolean wantsHessian() {
        return true;
    }

    public void setShift(double EnShift){
        this.EnShift = EnShift;
    }

}
