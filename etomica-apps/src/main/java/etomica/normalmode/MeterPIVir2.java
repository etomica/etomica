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

public class MeterPIVir2 implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms;
    protected Vector[] rc;
    protected double EnShift, sdot;
    protected Vector drShift;
    protected final Vector[] latticePositions;

    public MeterPIVir2(PotentialCompute pcP1, double temperature, Box box) {
        int nData = 7;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pcP1 = pcP1;
        this.beta = 1/temperature;
        this.box = box;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        this.EnShift = 0;
        drShift = box.getSpace().makeVector();

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (IMolecule m : box.getMoleculeList()) {
            latticePositions[m.getIndex()].E(m.getChildList().get(0).getPosition());
        }
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        if (box.getMoleculeList().size() > 1) drShift = computeShift();
        rHr = 0;
        pcP1.computeAll(true, this);
        Vector[] forces = pcP1.getForces();
        double vir = 0;
        double sumr2 = 0;
        Vector drj = box.getSpace().makeVector();

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                drj.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                drj.PE(drShift);
                vir -= forces[atom.getLeafIndex()].dot(drj);
                sumr2 += drj.squared();
            }
        }
        double U = pcP1.getLastEnergy();
        double Fx = -vir;
        x[0] = Fx;
        x[1] = rHr;
        x[2] = Fx*Fx;
        x[3] = U;
        x[4] = dim*(numAtoms)/2.0/beta + U + 0.5*Fx;
        x[5] = -dim*(numAtoms)*sdot + U - sdot*beta*Fx;

//        double sdot2 = -1/beta/beta/(rHr+1/beta);
        double sdot2 = -1.0/beta*Fx/(Fx-rHr);
//        System.out.println(sdot2);
//        x[6] = -sdot2 + U - sdot2*beta*Fx;
        x[6] = 2*rHr*Fx*Fx/(Fx-rHr)/(Fx-rHr) + U - sdot2*beta*Fx;
//        System.out.println(x[3] + " " + x[4] + " "+ x[6]);
        return data;
    }


    public void fieldComputeHessian(int i, Tensor Hii) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(ri);
        Hii.transform(tmpV);
        rHr += ri.dot(tmpV); //rHr
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        int moleculeIndexI = box.getLeafList().get(i).getParentGroup().getIndex();
        int moleculeIndexJ = box.getLeafList().get(j).getParentGroup().getIndex();
        Vector dri = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        dri.Ev1Mv2(box.getLeafList().get(i).getPosition(), latticePositions[moleculeIndexI]);
        drj.Ev1Mv2(box.getLeafList().get(j).getPosition(), latticePositions[moleculeIndexJ]);
        dri.PE(drShift);
        drj.PE(drShift);
        Vector drij = box.getSpace().makeVector();
        drij.Ev1Mv2(dri, drj);
        Hij.transform(drij);

        rHr += drij.dot(dri)  - drij.dot(drj) ;
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

    public void setEnShift(double E) { this.EnShift = E; }

    public void setSdot(double sdot) { this.sdot = sdot; }

    protected Vector computeShift() {
        drShift.E(0);
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                Vector drTmp = box.getSpace().makeVector();
                drTmp.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                drShift.PE(drTmp);
            }
        }
        drShift.TE(-1.0/box.getLeafList().size());
        return drShift;
    }

}
