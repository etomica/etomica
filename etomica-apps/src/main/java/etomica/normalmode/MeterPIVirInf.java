package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVirInf implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected double beta;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms;
    protected Vector[] rc;
    protected double omega, hbar, R, cothR, tanhR2, fac1, fac2, fac3;

    public MeterPIVirInf(PotentialCompute pcP1harm, PotentialCompute pcP1ah, double temperature, Box box, int nBeads, double omega, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pcP1harm = pcP1harm;
        this.pcP1ah = pcP1ah;
        this.omega = omega;
        this.hbar = hbar;
        this.beta = 1/temperature;
        this.box = box;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        this.R = beta*hbar*omega/nBeads;
        this.tanhR2 = Math.tanh(R/2);
        this.cothR = 1/tanhR2;
        this.fac1 = -R*R/2/beta*(2+1/Math.sinh(R/2)/Math.sinh(R/2));
        this.fac2 = R*R/4/beta*(1-1/Math.sin(R)/Math.sin(R)) + R*cothR/beta;
        this.fac3 = -R*R*cothR*cothR/4/beta;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        pcP1harm.computeAll(false);
//        pcP1ah.computeAll(true); //no Cv
        pcP1ah.computeAll(true, this); //with Cv
        Vector[] forces = pcP1ah.getForces();
        double vir = 0;

        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                vir -= forces[atom.getLeafIndex()].dot(ri);
            }
        }

        x[0] =   pcP1ah.getLastEnergy() + R/tanhR2*pcP1harm.getLastEnergy() + R*cothR/2.0*vir;//En
        x[1] =  fac1*pcP1harm.getLastEnergy() - fac2*vir + fac3*rHr + x[0]*x[0];
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
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
