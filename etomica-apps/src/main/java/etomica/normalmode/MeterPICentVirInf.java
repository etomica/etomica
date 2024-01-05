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

public class MeterPICentVirInf implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected double beta;
    protected int nBeads;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rc;
    protected int dim;
    protected int numAtoms;
    protected double omega, R, cothR, tanhR_2, hbar;
    protected double fac0, fac1, fac2, fac3, fac4;

    public MeterPICentVirInf(PotentialCompute pcP1harm, PotentialCompute pcP1ah, double temperature, Box box, int nBeads, double omega, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        this.pcP1harm = pcP1harm;
        this.pcP1ah = pcP1ah;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        this.omega = omega;
        this.hbar = hbar;
        this.box = box;
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        this.R = beta*hbar*omega/nBeads;
        this.tanhR_2 = Math.tanh(R/2);
        this.cothR = Math.cosh(R)/Math.sinh(R);
        double coshR = Math.cosh(R);
        double sinhR = Math.sinh(R);
        double massRing=1;
        this.fac0 = dim*numAtoms*R*R/2/beta/beta/sinhR/sinhR;
        this.fac1 = R/beta*massRing*omega*omega*tanhR_2*(1+2*coshR/sinhR/sinhR);
        this.fac2 = -R*R/2/beta*(2+1/Math.sinh(R/2)/Math.sinh(R/2));
        this.fac3 = R*R/4/beta*(1-1/sinhR/sinhR) + R*cothR/2/beta;
        this.fac4 = -R*R*cothR*cothR/4/beta;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        Vector rirc = box.getSpace().makeVector();
        rHr = 0;
        double vir = 0;
        pcP1harm.computeAll(false);
//        pcP1ah.computeAll(true); //no Cv
        pcP1ah.computeAll(true, this); //with Cv
        Vector[] forces = pcP1ah.getForces();
        double sumRc2 = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            sumRc2 += rc[molecule.getIndex()].squared();
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                rirc.Ev1Mv2(ri, rc[molecule.getIndex()]);
                box.getBoundary().nearestImage(rirc);
                vir -= forces[atom.getLeafIndex()].dot(rirc);
            }
        }
        sumRc2 /=numAtoms;
        double mass = 1;
        if(numAtoms==1){
            x[0] = dim/2.0/beta*R*cothR + pcP1ah.getLastEnergy() +  R/tanhR_2 *pcP1harm.getLastEnergy() + R*cothR/2.0*vir
                    - mass*omega*omega*cothR* tanhR_2 *sumRc2;
        } else {
//            x[0] = -dim*(nBeads-1)/(2.0*nBeads)/beta - dim / 2.0 / beta + dim * (nBeads - 1.0) / 2.0 / nBeads / beta + dim * numAtoms / 2.0 / beta + pcP1.getLastEnergy() + 1.0 / 2.0 * vir; //En
//            x[0] = dim * (numAtoms-1) / 2.0 / beta + pcP1.getLastEnergy() + 1.0 / 2.0 * vir; //En
        }
        x[1] = fac0 + fac1*sumRc2 + fac2*pcP1harm.getLastEnergy()
                - fac3*vir + fac4*rHr + x[0]*x[0];
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
        Vector tmpVecI = box.getSpace().makeVector();
        Vector tmpVecJ = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(ri, CenterOfMass.position(box, box.getLeafList().get(i).getParentGroup()));
        tmpVecJ.Ev1Mv2(rj, CenterOfMass.position(box, box.getLeafList().get(j).getParentGroup()));
        box.getBoundary().nearestImage(tmpVecI);
        box.getBoundary().nearestImage(tmpVecJ);
        Hij.transform(tmpVecI);
        rHr += tmpVecI.dot(tmpVecJ);
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
