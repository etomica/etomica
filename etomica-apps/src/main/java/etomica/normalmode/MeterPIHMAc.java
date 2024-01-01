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

    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
        double vir = 0;
        double virc = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();
//        Vector dr0Ref = box.getSpace().makeVector();

//        pcP1.computeAll(true, null); // no Cv (rHr=0)
        pcP1.computeAll(true, this); //with Cv
        Vector[] forces = pcP1.getForces();
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
        }
        if(numAtoms == 1) {
            x[0] = dim/beta + pcP1.getLastEnergy() + 1.0 / 2.0 * vir; //En in 1D
        } else {
            x[0] = -dim/2.0/beta + dim*numAtoms/beta + pcP1.getLastEnergy() + 1.0/2.0*vir; //En in 3D
        }
        x[1] = dim/beta/beta + 1.0/4.0/beta*(-3.0*vir - 2.0*virc - rHr)  + x[0]*x[0];
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
