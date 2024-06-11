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

public class MeterPICentVirFD implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected int nBeads;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected Vector[] rc;
    protected int dim;
    protected int numAtoms;
    protected double EnShift;
    protected double volume;
    protected double dbeta;
    protected Vector[] rOrig;

    public MeterPICentVirFD(PotentialCompute pcP1, double temperature, int nBeads, Box box, double dbeta) {
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
        this.dbeta = dbeta;
        rOrig = box.getSpace().makeVectorArray(numAtoms*nBeads);
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
        }
        double vir = 0;
        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        Vector rirc = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
                rirc.Ev1Mv2(rOrig[atom.getLeafIndex()], rc[molecule.getIndex()]);
                box.getBoundary().nearestImage(rirc);
                vir -= forces[atom.getLeafIndex()].dot(rirc);
            }
        }
        x[0] = dim*numAtoms/2.0/beta + pcP1.getLastEnergy() + 1.0/2.0 * vir - EnShift;

        double dbUdb_p = dbUdb(beta+dbeta);
        double dbUdb_m = dbUdb(beta-dbeta);
        double d2bUdb2 = (dbUdb_p - dbUdb_m)/(2*dbeta);

        x[1] = dim*numAtoms/2.0/beta/beta - d2bUdb2 +  x[0]*x[0];

//        System.out.println("FD-CV: "+ (-d2bUdb2));
//        System.out.println();

        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setEnShift(double E) { this.EnShift = E; }

    public double dbUdb(double beta_new) {
        Vector rirc = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                rirc.Ev1Mv2(atom.getPosition(), rc[molecule.getIndex()]);
                box.getBoundary().nearestImage(rirc);
                rirc.TE(Math.sqrt(beta_new/beta));
                atom.getPosition().Ev1Pv2(rc[molecule.getIndex()], rirc);
            }
        }

        double vir = 0;
        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                rirc.Ev1Mv2(ri, rc[molecule.getIndex()]);//no need to recompute rc, because it isn't changed!
                box.getBoundary().nearestImage(rirc);
                vir -= forces[atom.getLeafIndex()].dot(rirc);
            }
        }

        double dbUdb = pcP1.getLastEnergy() + 1.0 / 2.0 * vir;

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                atom.getPosition().E(rOrig[atom.getLeafIndex()]);
            }
        }

        return dbUdb;
    }
}
