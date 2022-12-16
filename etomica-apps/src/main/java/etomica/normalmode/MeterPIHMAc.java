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

public class MeterPIHMAc implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1;
    protected double betaN, beta;
    protected int nBeads;
    protected double rHr;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rc;
    protected double EnShift;
    protected double dim;
    protected int numAtoms;
    protected Vector[] latticePositions;

    public MeterPIHMAc(PotentialCompute pcP1, double betaN, int nBeads, Box box) {
        int nData = 1;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        this.pcP1 = pcP1;
        this.betaN = betaN;
        this.nBeads = nBeads;
        beta = this.betaN*this.nBeads;
        this.box = box;
        this.EnShift = 0;
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();

        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        double vir = 0;
        double virc = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();

        pcP1.computeAll(true);
//        pcP1.computeAll(true, this);//it needs rc
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                dri.Ev1Mv2(ri, latticePositions[molecule.getIndex()]);
                drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
                box.getBoundary().nearestImage(dri);
                box.getBoundary().nearestImage(drc);

                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
                virc -= forces[atom.getLeafIndex()].dot(drc);
            }
        }

        x[0] = dim*numAtoms/beta + pcP1.getLastEnergy() + 1.0/2.0*vir;// - EnShift; //En
//        x[1] = 1.0/beta/beta + 1.0/4.0/beta*(3.0*vir - 2.0*virc - rHr); //Cvn/kb^2, without Var
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector ri = box.getLeafList().get(i).getPosition();
        Vector rj = box.getLeafList().get(j).getPosition();
        Vector tmpV = box.getSpace().makeVector();
        int moleculeIndex = box.getLeafList().get(i).getParentGroup().getIndex();
        tmpV.Ev1Mv2(rj, rc[moleculeIndex]);
        tmpV.ME(rc[moleculeIndex]);
        Hij.transform(tmpV);
        rHr += ri.dot(tmpV) - 2.0*rc[moleculeIndex].dot(tmpV);
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
