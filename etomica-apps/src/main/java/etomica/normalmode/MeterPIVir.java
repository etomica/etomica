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

public class MeterPIVir implements IDataSource, PotentialCallback {
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

    public MeterPIVir(PotentialCompute pcP1, double temperature, Box box) {
        int nData = 2;
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
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        rHr = 0;
//        pcP1.computeAll(true); // no Cv (rHr=0)
        pcP1.computeAll(true, this); //with Cv
        Vector[] forces = pcP1.getForces();
        double vir = 0;

        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                vir -= forces[atom.getLeafIndex()].dot(ri);
            }
        }

        if (numAtoms == 1) { // e.g., HO or EC
            x[0] =   pcP1.getLastEnergy() + 1.0/2.0*vir;
            x[1] = 1.0/4.0/beta*(-3.0*vir - rHr) + x[0]*x[0];
        } else {
            x[0] =   dim/2.0/beta + pcP1.getLastEnergy() + 1.0/2.0*vir;
            x[1] = dim/2.0/beta/beta + 1.0/4.0/beta*(-3.0*vir - rHr) + x[0]*x[0];
        }
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
