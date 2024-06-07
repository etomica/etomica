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

public class MeterPIHMAcFD implements IDataSource, PotentialCallback {
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
    protected double dbeta;
    protected Vector[] rOrig;

    public MeterPIHMAcFD(PotentialCompute pcP1, double temperature, int nBeads, Box box, double dbeta) {
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
        this.dbeta = dbeta;
        rOrig = box.getSpace().makeVectorArray(numAtoms*nBeads);
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        if (numAtoms > 1) {
            drShift = computeShift();
        }
        rHr = 0;
        double vir = 0;
        Vector drc = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();

        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            rc[molecule.getIndex()] = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
                dri.Ev1Mv2(rOrig[atom.getLeafIndex()], latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                box.getBoundary().nearestImage(dri);
                drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
                drc.PE(drShift);
                box.getBoundary().nearestImage(drc);

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
        double d2bUdb2 = (dbUdb_p - dbUdb_m)/(2*dbeta);

        if (numAtoms == 1) {
            x[1] = dim/beta/beta - d2bUdb2  + x[0]*x[0];
        } else {
            x[1] = -dim/2.0/beta/beta + dim*numAtoms/beta/beta - d2bUdb2  + x[0]*x[0];
        }

        return data;
    }

    public double dbUdb(double beta_new) {
        Vector dri = box.getSpace().makeVector();
        Vector drc = box.getSpace().makeVector();
        Vector tmpV = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                //Need fix shifts for LJ!
                dri.PE(drShift);
                box.getBoundary().nearestImage(dri);
                drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
                drc.PE(drShift);
                box.getBoundary().nearestImage(drc);
                tmpV.Ev1Mv2(dri, drc);
                box.getBoundary().nearestImage(tmpV);
                tmpV.TE(Math.sqrt(beta_new/beta));
                tmpV.PEa1Tv1(Math.sqrt(beta/beta_new), drc);
                atom.getPosition().E(tmpV);
            }
        }

        double vir = 0;
        pcP1.computeAll(true, null);
        Vector[] forces = pcP1.getForces();
        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                Vector ri = atom.getPosition();
                dri.Ev1Mv2(ri, latticePositions[molecule.getIndex()]);
                //Need fix shifts for LJ!
                dri.PE(drShift);
                box.getBoundary().nearestImage(dri);
                drc.Ev1Mv2(rc[molecule.getIndex()], latticePositions[molecule.getIndex()]);
                drc.PE(drShift);
                box.getBoundary().nearestImage(drc);
                vir -= forces[atom.getLeafIndex()].dot(dri);
                vir += 2.0*forces[atom.getLeafIndex()].dot(drc);
            }
        }

        double dbUdb = pcP1.getLastEnergy() + 1.0/2.0*vir;

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                atom.getPosition().E(rOrig[atom.getLeafIndex()]);
            }
        }

        return dbUdb;
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

    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) {
            return box.getSpace().makeVector();
        }
        int n = box.getMoleculeList().get(0).getChildList().size();
        Vector shift0 = box.getSpace().makeVector();
        Boundary boundary = box.getBoundary();
        Vector dr = box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        dr.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        boundary.nearestImage(dr);
        shift0.Ea1Tv1(-1, dr);
        // will shift ring0 back to lattice site; everything should be close and PBC should lock in
        // now determine additional shift needed to bring back to original COM
        Vector totalShift = box.getSpace().makeVector();
        for (int j = 0; j < box.getMoleculeList().size(); j++) {
            IMolecule m = box.getMoleculeList().get(j);
            for (int i = 0; i < n; i++) {
                Vector r = m.getChildList().get(i).getPosition();
                dr.Ev1Mv2(r, latticePositions[j]);
                dr.PE(shift0);
                boundary.nearestImage(dr);
                totalShift.PE(dr);
            }
        }
        totalShift.TE(-1.0/box.getLeafList().size());
        totalShift.PE(shift0);
        return totalShift;
    }

}
