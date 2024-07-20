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
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 Only works for PIMD where atoms are not wrapped back to box. Needs fixed for PIMC.
 */
public class MeterPIHMAcFD implements IDataSource, PotentialCallback {
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
    protected Vector[] latticePositions;
    protected Vector drShift;
    protected double EnShift;
    protected double dbeta;
    protected Vector[] rOrig;
    protected  int n;

    public MeterPIHMAcFD(PotentialCompute pcP1, double temperature, int nBeads, Box box, double dbeta) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.box = box;
        this.pcP1 = pcP1;
        this.nBeads = nBeads;
        this.beta = 1/temperature;
        this.dbeta = dbeta;

        dim = box.getSpace().D();
        numAtoms = box.getMoleculeList().size();

        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());
        rc = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

        drShift = box.getSpace().makeVector();
        EnShift = 0;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        if (numAtoms > 1) drShift = computeShift();
        double vir = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();
        pcP1.computeAll(true);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()].E(0);
            for (IAtom atom : molecule.getChildList()) {
                rc[molecule.getIndex()].PE(atom.getPosition());
            }
            rc[molecule.getIndex()].TE(1.0/nBeads);
            drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
            drc.PE(drShift);
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
            }
        }

        if (numAtoms == 1) {
            x[0] = dim/beta + pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift;
        } else {
            x[0] = -dim/2.0/beta + dim*numAtoms/beta + pcP1.getLastEnergy() + 1.0/2.0*vir - EnShift;
        }

        double dbUdb_p = dbUdb(beta+dbeta);
        double dbUdb_m = dbUdb(beta-dbeta);
        double d2bUdb2 = (dbUdb_p-dbUdb_m)/2/dbeta;

        if (numAtoms == 1) {
            x[1] = dim/beta/beta - d2bUdb2  + x[0]*x[0];
        } else {
            x[1] = -dim/2.0/beta/beta + dim*numAtoms/beta/beta - d2bUdb2  + x[0]*x[0];
        }
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
        //Scale coords
        Vector dri = box.getSpace().makeVector();
        Vector drc = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()].E(0);
            for (IAtom atom : molecule.getChildList()) {
                rc[molecule.getIndex()].PE(atom.getPosition());
            }
            rc[molecule.getIndex()].TE(1.0/nBeads);
            drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
            drc.PE(drShift);
            for (IAtom atom : molecule.getChildList()) {
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                Vector tmpV = box.getSpace().makeVector();
                tmpV.Ev1Mv2(dri, drc);
                tmpV.TE(Math.sqrt(beta_new/beta));
                tmpV.PEa1Tv1(Math.sqrt(beta/beta_new), drc);
                atom.getPosition().Ev1Pv2(latticePositions[molecule.getIndex()], tmpV);
            }
        }

        // recompute
        double vir = 0;
        pcP1.computeAll(true);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()].E(0);
            for (IAtom atom : molecule.getChildList()) {
                rc[molecule.getIndex()].PE(atom.getPosition());
            }
            rc[molecule.getIndex()].TE(1.0/nBeads);
            drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
            drc.PE(drShift);
            for (IAtom atom : molecule.getChildList()) {
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
            }
        }

        double dbUdb = pcP1.getLastEnergy() + 1.0/2.0*vir;
        // Bring atoms back to their original coord
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                atom.getPosition().E(rOrig[atom.getLeafIndex()]);
            }
        }

        return dbUdb;
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
