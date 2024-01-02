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

public class MeterPIVirHarmStagingInf implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1ah;
    protected double beta;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms, nBeads;
    protected Vector[] rc;
    protected double omega, R, cothR, hbar;
    double[] D;

    public MeterPIVirHarmStagingInf(PotentialCompute pcP1ah, double temperature, Box box, int nBeads, double omega, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pcP1ah = pcP1ah;
        this.omega = omega;
        this.hbar = hbar;
        this.beta = 1/temperature;
        this.box = box;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        this.nBeads = nBeads;
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        this.R = beta*hbar*omega/nBeads;
        this.cothR = 1.0/Math.tanh(R);
        D = new double[nBeads];
        double Ai = 0;
        for (int i = nBeads - 1; i > 0; i--){
            Ai = Ai*Math.sinh(i*R)/Math.sinh((i+1)*R) + Math.sinh(R)/Math.sinh((i+1)*R);
            D[i] = Ai + Math.sinh(i*R)/Math.sinh(beta*hbar*omega);
        }
        D[0] = 1;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        pcP1ah.computeAll(true); //no Cv
//        pcP1ah.computeAll(true, this); //with Cv
        Vector[] forces = pcP1ah.getForces();
        double vir = 0;
        double vir2 = 0;

        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atomi : molecule.getChildList()) {
                Vector ri = atomi.getPosition();
                vir -= forces[atomi.getLeafIndex()].dot(ri);
                for (IAtom atomj : molecule.getChildList()) {
                    Vector rj = atomj.getPosition();
                    int ij = atomi.getIndex() - atomj.getIndex();
                    double DD;
                    if (ij >= 0) {
                        DD = D[ij];
                    } else {
                        DD = D[nBeads+ij];
                    }
                    vir2 -= DD * forces[atomi.getLeafIndex()].dot(rj);
                }
            }
        }
        x[0] =  1.0/2.0*hbar*omega/Math.tanh(beta*hbar*omega/2) + 1.0/2.0*R*cothR*vir - 1.0/2.0*R/Math.tanh(beta*hbar*omega/2)*vir2 + pcP1ah.getLastEnergy();
        x[1] = 1.0/4.0/beta*(-3.0*vir - rHr) + x[0]*x[0];
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
