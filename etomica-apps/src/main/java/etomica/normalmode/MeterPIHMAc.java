package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
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
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 Only works for PIMD where atoms are not wrapped back to box. Needs fixed for PIMC.
 */
public class MeterPIHMAc implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected int nBeads;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rc;
    protected int dim;
    protected int numAtoms;
    protected Vector[] latticePositions;
    protected Vector drShift;
    protected double EnShift;

    public MeterPIHMAc(PotentialCompute pcP1, double temperature, int nBeads, Box box) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pcP1 = pcP1;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        this.box = box;
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }
        this.drShift = box.getSpace().makeVector();
        this.EnShift = 0;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        if (numAtoms > 1) {
            drShift = computeShift();
        }

        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()].E(0);
            for (IAtom atom : molecule.getChildList()) {
                rc[molecule.getIndex()].PE(atom.getPosition());
            }
            rc[molecule.getIndex()].TE(1.0 / nBeads);
        }

        rHr = 0;
        double vir = 0;
        double virc = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();

        pcP1.computeAll(true, this); //with Cv
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
            drc.PE(drShift);
            box.getBoundary().nearestImage(drc);

            for (IAtom atom : molecule.getChildList()) {
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                box.getBoundary().nearestImage(dri);
                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
                virc -= forces[atom.getLeafIndex()].dot(drc);
            }
        }
        if (numAtoms == 1) {
            x[0] = dim/beta + pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift;
            x[1] = dim/beta/beta + 1.0/4.0/beta*(-3.0*vir - 2.0*virc - rHr)  + x[0]*x[0];
        } else {
            x[0] = -dim/2.0/beta + dim*numAtoms/beta + pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift;
            x[1] = -dim/2.0/beta/beta + dim*numAtoms/beta/beta + 1.0/4.0/beta*(-3.0*vir - 2.0*virc - rHr)  + x[0]*x[0];
        }

        return data;
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

    public void fieldComputeHessian(int i, Tensor Hii) { // in general potential, Hij is the Hessian between same beads of atom i and j
        int moleculeIndexI = box.getLeafList().get(i).getParentGroup().getIndex();
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector dri = box.getSpace().makeVector();
        dri.Ev1Mv2(ri, latticePositions[moleculeIndexI]);
        Vector drcI = box.getSpace().makeVector();
        drcI.Ev1Mv2(rc[box.getLeafList().get(i).getParentGroup().getIndex()], latticePositions[moleculeIndexI]);
        Vector tmpVecI = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(dri, drcI);
        tmpVecI.ME(drcI);
        Hii.transform(tmpVecI);
        rHr += tmpVecI.dot(dri) - 2*tmpVecI.dot(drcI);
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
        dri.PE(drShift);
        drj.PE(drShift);
        box.getBoundary().nearestImage(dri);
        box.getBoundary().nearestImage(drj);

        Vector drcI = box.getSpace().makeVector();
        Vector drcJ = box.getSpace().makeVector();
        drcI.Ev1Mv2(rc[box.getLeafList().get(i).getParentGroup().getIndex()], latticePositions[moleculeIndexI]);
        drcJ.Ev1Mv2(rc[box.getLeafList().get(j).getParentGroup().getIndex()], latticePositions[moleculeIndexJ]);
        drcI.PE(drShift);
        drcJ.PE(drShift);
        box.getBoundary().nearestImage(drcI);
        box.getBoundary().nearestImage(drcJ);

        Vector tmpVecI = box.getSpace().makeVector();
        Vector tmpVecJ = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(dri, drcI);
        tmpVecJ.Ev1Mv2(drj, drcJ);
        tmpVecI.ME(drcI);
        tmpVecJ.ME(drcJ);

        Vector tmpVecIJ = box.getSpace().makeVector();
        tmpVecIJ.Ev1Mv2(tmpVecI, tmpVecJ);
        Hij.transform(tmpVecIJ);
        rHr += tmpVecIJ.dot(tmpVecI) - tmpVecIJ.dot(tmpVecJ);
    }


    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) return box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        Vector drShift = box.getSpace().makeVector();
        drShift.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        drShift.TE(-1);
        return drShift;
    }
}