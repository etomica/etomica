package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIVirFD implements IDataSource, PotentialCallback {
    protected Box box;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int dim;
    protected int numAtoms;
    protected double EnShift;
    protected double volume;
    protected double dbeta;
    protected Vector[] rOrig;

    public MeterPIVirFD(PotentialCompute pcP1, double temperature, Box box, double dbeta) {
        int nData = 2;
        this.volume = box.getBoundary().volume();
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pcP1 = pcP1;
        this.beta = 1/temperature;
        this.box = box;
        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();
        this.EnShift = 0;
        this.dbeta = dbeta;
        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        double vir = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
                vir -= forces[atom.getLeafIndex()].dot(rOrig[atom.getLeafIndex()]);
            }
        }

        if (numAtoms == 1) {
            x[0] = pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift;
        } else {
            x[0] = dim/2.0/beta + pcP1.getLastEnergy() + 1.0/2.0*vir- EnShift;
        }

        double dbUdb_p = dbUdb(beta+dbeta);
        double dbUdb_m = dbUdb(beta-dbeta);
        double d2bUdb2 = (dbUdb_p - dbUdb_m)/(2*dbeta);


        if (numAtoms == 1) {
            x[1] = - d2bUdb2 +  x[0]*x[0];
        } else {
            x[1] = dim/2.0/beta/beta - d2bUdb2 +  x[0]*x[0];
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

    public double dbUdb(double beta_new) {
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                atom.getPosition().TE(Math.sqrt(beta_new/beta));
            }
        }

        double vir = 0;
        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                vir -= forces[atom.getLeafIndex()].dot(ri);
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
