/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAReal2 implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected double EnShift;
    protected final MCMoveHOReal2 move;
    protected int nShifts = 0;
    protected double rHr, dim;
    protected Vector[] rdot;
    protected double drdotHdrdot;
    protected Vector[] latticePositions;
    protected int numAtoms, nBeads;

    public MeterPIHMAReal2(PotentialMasterBonding pmBonding, PotentialCompute pcP1, int nBeads, double beta, MCMoveHOReal2 move) {
        this.move = move;
        this.beta = beta;
        int nData = 1;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.EnShift = 0;
        Box box = move.getBox();
        this.nBeads = nBeads;
        rdot = new Vector[this.nBeads];
        for (int i=0; i<this.nBeads; i++) {
            rdot[i] = box.getSpace().makeVector();
        }
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }
        numAtoms = box.getMoleculeList().size();
        dim = box.getSpace().D();
    }

    public void setNumShifts(int nShifts) {
        if (nShifts < 0) throw new RuntimeException("nShifts cannot be negative");
        this.nShifts = nShifts;
    }

    @Override
    public IData getData() {
        drdotHdrdot = 0 ;
        Box box = move.getBox();
        double[] x = data.getData();
//        System.out.println("******** HMA *************");

        pmBonding.computeAll(true);
        pcP1.computeAll(true);

        double En0 = dim*nBeads*numAtoms/2.0/beta + pcP1.getLastEnergy() - pmBonding.getLastEnergy();
        double Cvn0 = nBeads/2.0/beta/beta - 2.0/beta*pmBonding.getLastEnergy();

        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();

        IMoleculeList molecules = box.getMoleculeList();

        double[] gamma = move.getGamma();
        double[] dGamma = move.getdGamma();
        double[][] centerCoefficients = move.getCenterCoefficients();
        double[] f11 = centerCoefficients[0];
        double[] f1N = centerCoefficients[1];
        double[][] dCenterCoefficients = move.getDCenterCoefficients();
        double[] df11 = dCenterCoefficients[0];
        double[] df1N = dCenterCoefficients[1];
        double[][] d2CenterCoefficients = move.getD2CenterCoefficients();
        double[] d2f11 = d2CenterCoefficients[0];
        double[] d2f1N = d2CenterCoefficients[1];

        Vector tmp_r = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        Vector v = box.getSpace().makeVector();
        Vector a = box.getSpace().makeVector();
        Vector v0 = box.getSpace().makeVector();
        Vector a0 = box.getSpace().makeVector();
        Vector vPrev = box.getSpace().makeVector();
        Vector aPrev = box.getSpace().makeVector();
        double En = 0;
        double Cvn = 0;
        rHr = 0;

        Vector dr0ref = box.getSpace().makeVector();

        int ns = nShifts+1;
        for (int i = 0; i < molecules.size(); i++) {
            IAtomList beads = molecules.get(i).getChildList();
            if (ns > beads.size() || (beads.size() / ns) * ns != beads.size()) {
                throw new RuntimeException("# of beads must be a multiple of (# of shifts + 1)");
            }
            for (int indexShift=0; indexShift<beads.size(); indexShift += beads.size()/ns) {
                Vector dr0 = box.getSpace().makeVector();
                dr0.Ev1Mv2(beads.get(indexShift).getPosition(), latticePositions[i]);
                if(i==0){
                    dr0ref.E(dr0);
                }
                dr0.ME(dr0ref);
                box.getBoundary().nearestImage(dr0);

                Vector drPrev = box.getSpace().makeVector();

                for (int j = 0; j < beads.size(); j++) {
                    int aj = (j + indexShift) % beads.size();
                    IAtom atomj = beads.get(aj);
                    if (i != 0 || j != 0) {
                        En -= dim*gamma[j];
                        Cvn += dim*dGamma[j];
                    }
                    drj.Ev1Mv2(atomj.getPosition(), latticePositions[i]);
                    drj.Ev1Mv2(drj, dr0ref);
                    box.getBoundary().nearestImage(drj);
                    tmp_r.E(drj);
                    if (j > 0) {
                        tmp_r.PEa1Tv1(-f11[j], drPrev);
                        tmp_r.PEa1Tv1(-f1N[j], dr0);
                    }
                    v.Ea1Tv1(gamma[j], tmp_r);
                    a.Ea1Tv1(dGamma[j]+gamma[j]*gamma[j], tmp_r);
                    if (j > 0) {
                        v.PEa1Tv1(df11[j], drPrev);
                        v.PEa1Tv1(df1N[j], dr0);
                        v.PEa1Tv1(f11[j], vPrev);
                        v.PEa1Tv1(f1N[j], v0);

                        a.PEa1Tv1(d2f11[j], drPrev);
                        a.PEa1Tv1(d2f1N[j], dr0);
                        a.PEa1Tv1(2.0*df11[j], vPrev);
                        a.PEa1Tv1(2.0*df1N[j], v0);
                        a.PEa1Tv1(f11[j], aPrev);
                        a.PEa1Tv1(f1N[j], a0);

                        Vector tmpV = box.getSpace().makeVector();
                        tmpV.Ev1Mv2(v, vPrev);
                        Cvn -= beta/nBeads*move.mass*move.omegaN*move.omegaN*(tmpV.squared());
                        if (j == nBeads-1){
                            tmpV.Ev1Mv2(v0, v);
                            Cvn -= beta/nBeads*move.mass*move.omegaN*move.omegaN*(tmpV.squared());
                        }
                    }
                    int jj = atomj.getLeafIndex();
                    En -= beta*forcesU[jj].dot(v);
                    Cvn += beta*forcesU[jj].dot(a) + 2.0*forcesU[jj].dot(v) ;
                    if (nBeads > 1) {
                        En -= beta*forcesK[jj].dot(v);
                        Cvn += beta*forcesK[jj].dot(a) - 2.0*forcesK[jj].dot(v);
                    }

                    drPrev.E(drj);
                    if (j == 0) {
                        v0.E(v);
                        a0.E(a);
                    }
                    rdot[j].E(v);
                    vPrev.E(v);
                    aPrev.E(a);
                }
            }
        }

//        pcP1.computeAll(true, this);
        Cvn -= beta*drdotHdrdot;
        x[0] = En0 + En/ns - EnShift;
//        x[1] = Cvn0 + Cvn/ns;

        return data;
    }


    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector tmpV = move.getBox().getSpace().makeVector();
        tmpV.E(rdot[j]);
        Hij.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV); //drdot.H.drdot
    }

    public boolean wantsHessian() {
        return true;
    }


    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setShift(double EnShift){
        this.EnShift = EnShift;
    }

}