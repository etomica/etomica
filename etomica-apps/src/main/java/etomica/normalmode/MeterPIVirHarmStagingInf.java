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
    protected double omega, R, cothR, hbar, tanhR_2, fac1, fac2, fac3, fac4;
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
        D = new double[nBeads];
        double Ai = 0;
        for (int i = nBeads - 1; i > 0; i--){
            Ai = Ai*Math.sinh(i*R)/Math.sinh((i+1)*R) + Math.sinh(R)/Math.sinh((i+1)*R);
            D[i] = Ai + Math.sinh(i*R)/Math.sinh(beta*hbar*omega);
        }
        D[0] = 1;

        this.cothR = 1.0/Math.tanh(R);
        this.tanhR_2 = Math.tanh(R/2);
        this.fac1 = -R*R/2/beta*Math.cosh(R)/Math.sinh(R/2)/Math.sinh(R/2);
        this.fac2 = R*R/4/beta*(1-1/Math.sinh(R)/Math.sinh(R)) + R*cothR/beta;
        this.fac3 = -R*R*cothR*cothR/4/beta;
        this.fac4 = fac1*tanhR_2/R;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
//        pcP1ah.computeAll(true); //no Cv
        pcP1ah.computeAll(true, this); //with Cv
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
        double Einf = hbar*omega/2/Math.tanh(beta*hbar*omega/2);
        x[0] =  R/2*cothR*vir + Einf*(1 - beta/nBeads*vir2) + pcP1ah.getLastEnergy();
        x[1] = -fac2*vir + fac4*Einf*(1 - beta/nBeads*vir2) + fac3*rHr + x[0]*x[0];

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
