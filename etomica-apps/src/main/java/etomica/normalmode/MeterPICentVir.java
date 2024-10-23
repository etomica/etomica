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

/**
 Only works for PIMD where atoms are not wrapped back to box. Needs fixed for PIMC.
 */
public class MeterPICentVir implements IDataSource, PotentialCallback {
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
    protected double EnShift;
    protected double volume;

    public MeterPICentVir(PotentialCompute pcP1, double temperature, int nBeads, Box box) {
        int nData = 2;
        this.volume = box.getBoundary().volume();
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
        this.EnShift = 0;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        Vector rirc = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()].E(0);
            for (IAtom atom : molecule.getChildList()) {
                rc[molecule.getIndex()].PE(atom.getPosition());
            }
            rc[molecule.getIndex()].TE(1.0/nBeads);
        }
        rHr = 0;
        double vir = 0;
        pcP1.computeAll(true, this); //with Cv

        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                rirc.Ev1Mv2(ri, rc[molecule.getIndex()]);
                vir -= forces[atom.getLeafIndex()].dot(rirc);
            }
        }

        // same for both bound and unbound systems
        x[0] = dim*numAtoms/2.0/beta + pcP1.getLastEnergy() + 1.0 / 2.0 * vir - EnShift;
        x[1] = dim*numAtoms/2.0/beta/beta + 1.0/4.0/beta*(-3.0*vir - rHr) +  x[0]*x[0];
//        System.out.println("cv: " + x[0]);
//        x[2] = numAtoms/volume/beta - pcP1.getLastVirial()/dim/volume + vir/dim/volume; //Pressure
        return data;
    }

    public void fieldComputeHessian(int i, Tensor Hii) {
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector tmpVecI = box.getSpace().makeVector();
        Vector tmpVecI2 = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(ri, rc[box.getLeafList().get(i).getParentGroup().getIndex()]);
        tmpVecI2.E(tmpVecI);
        Hii.transform(tmpVecI);
        rHr += tmpVecI.dot(tmpVecI2);
    }


    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
        Vector tmpVecI = box.getSpace().makeVector();
        Vector tmpVecJ = box.getSpace().makeVector();
        tmpVecI.Ev1Mv2(ri, rc[box.getLeafList().get(i).getParentGroup().getIndex()]);
        tmpVecJ.Ev1Mv2(rj, rc[box.getLeafList().get(j).getParentGroup().getIndex()]);
        Vector tmpVecIJ = box.getSpace().makeVector();
        tmpVecIJ.Ev1Mv2(tmpVecI, tmpVecJ);
        Hij.transform(tmpVecIJ);
        rHr += tmpVecIJ.dot(tmpVecI) - tmpVecIJ.dot(tmpVecJ);
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
}
