package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVir implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialComputeField pcP1;
    protected double betaN, beta;
    protected int nBeads;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector ri, rj;


    public MeterPIVir(PotentialComputeField pcP1, double betaN, int nBeads, Box box) {
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
        ri = box.getSpace().makeVector();
        rj = box.getSpace().makeVector();
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        pcP1.computeAll(true, this);
        Vector[] forces = pcP1.getForces();
        double vir = 0;
        for (int i = 0; i < nBeads; i++){
            ri = box.getLeafList().get(i).getPosition();
            vir -= forces[i].dot(ri);
        }
        x[0] = pcP1.getLastEnergy() + 1.0/2.0*vir;//En
        x[1] = 1.0/4.0/beta*(-3.0*vir - rHr); //Cvn/kb^2, without Var
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        IAtom ai = box.getLeafList().get(i);
        IAtom aj = box.getLeafList().get(j);
        ri.E(ai.getPosition());
        rj.E(aj.getPosition());
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(rj);
        Hij.transform(tmpV);
        rHr += ri.dot(tmpV); //rHr
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

}
