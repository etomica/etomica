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

public class MeterPIHMAcInf implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected double beta;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rc;
    protected int dim;
    protected int numAtoms, nBeads;
    protected double omega, R, cothR, tanhR_2, hbar;
    protected double fac0, fac1, fac2, fac3, fac4, fac5;
    protected Vector[] latticePositions;

    public MeterPIHMAcInf(PotentialCompute pcP1harm, PotentialCompute pcP1ah, double temperature, Box box, int nBeads, double omega, double hbar) {
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
        this.cothR = 1/Math.tanh(R);
        double coshR = Math.cosh(R);
        double sinhR = Math.sinh(R);
        double massRing=1;
        this.fac0 = dim*numAtoms*R*R/beta/beta/sinhR/sinhR;
        this.fac1 = 2*R/beta*massRing*omega*omega*tanhR_2/sinhR/sinhR*(2*coshR-1);
        this.fac2 = -R*R/2/beta*coshR/Math.sinh(R/2)/Math.sinh(R/2);
        this.fac3 = R*R/4/beta*(1-1/sinhR/sinhR) + R*cothR/beta;
        this.fac4 = -R*R*cothR*cothR/4/beta;
        this.fac5 = R*R*cothR*cothR/2/beta;

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        Vector rcNI = box.getSpace().makeVector();
        rHr = 0;
        double vir = 0;
        double virc = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();

        pcP1harm.computeAll(false);
//        pcP1ah.computeAll(true); //no Cv
        pcP1ah.computeAll(true, this); //with Cv
        Vector[] forces = pcP1ah.getForces();


        double sumRc2 = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                dri.Ev1Mv2(ri, latticePositions[molecule.getIndex()]);
//                if (atom.getLeafIndex() == 0) {
//                    dr0Ref.E(ri);
//                }
//                dri.ME(dr0Ref);
                drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
//                drc.ME(dr0Ref);
                box.getBoundary().nearestImage(dri);
                box.getBoundary().nearestImage(drc);

                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
                virc -= forces[atom.getLeafIndex()].dot(drc);
            }


            rcNI.E(rc[molecule.getIndex()]);
            box.getBoundary().nearestImage(rcNI);
            sumRc2 += rcNI.squared();


        }

        double mass = 1;
        x[0] = dim/beta*R*cothR + pcP1ah.getLastEnergy() +  R/tanhR_2 *pcP1harm.getLastEnergy() + R*cothR/2*vir
                - 2*mass*omega*omega*cothR*tanhR_2*sumRc2;
        x[1] = fac0 + fac1*sumRc2 + fac2*pcP1harm.getLastEnergy() - fac3*vir - fac5*virc + fac4*rHr + x[0]*x[0];
        return data;
    }
    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        int moleculeIndexI = box.getLeafList().get(i).getParentGroup().getIndex();
        int moleculeIndexJ = box.getLeafList().get(j).getParentGroup().getIndex();
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
        Vector dri = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        dri.Ev1Mv2(ri, latticePositions[moleculeIndexI]);
        drj.Ev1Mv2(rj, latticePositions[moleculeIndexJ]);
        box.getBoundary().nearestImage(dri);
        box.getBoundary().nearestImage(drj);
        Vector drcI = box.getSpace().makeVector();
        Vector drcJ = box.getSpace().makeVector();
        drcI.Ev1Mv2(CenterOfMass.position(box, box.getLeafList().get(i).getParentGroup()), latticePositions[moleculeIndexI]);
        drcJ.Ev1Mv2(CenterOfMass.position(box, box.getLeafList().get(j).getParentGroup()), latticePositions[moleculeIndexJ]);
        box.getBoundary().nearestImage(drcI);
        box.getBoundary().nearestImage(drcJ);
        Vector tmpVecI = box.getSpace().makeVector();
        Vector tmpVecJ = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(dri, drcI);
        tmpVecJ.Ev1Mv2(drj, drcJ);
        tmpVecI.ME(drcI);
        tmpVecJ.ME(drcJ);
        Hij.transform(tmpVecJ);
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
