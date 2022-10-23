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
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputeField;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAReal2 implements IDataSource {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double beta;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected double EnShift;
    protected final MCMoveHOReal2 move;
    protected int nShifts = 0;

    public MeterPIHMAReal2(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double beta, MCMoveHOReal2 move) {
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
    }

    public void setNumShifts(int nShifts) {
        if (nShifts < 0) throw new RuntimeException("nShifts cannot be negative");
        this.nShifts = nShifts;
    }

    @Override
    public IData getData() {
        Box box = move.getBox();
        double[] x = data.getData();
//        System.out.println("******** HMA *************");

        int nBeads = box.getLeafList().size();
        double En0 = 0.5*nBeads/beta + pcP1.computeAll(true) - pmBonding.computeAll(true);
//        System.out.println("real: "+En);

        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();

        IMoleculeList molecules = box.getMoleculeList();

        double[] gamma = move.getGamma();
        double[][] centerCoefficients = move.getCenterCoefficients();
        double[] f11 = centerCoefficients[0];
        double[] f1N = centerCoefficients[1];
        double[][] dcenterCoefficients = move.getDCenterCoefficients();
        double[] df11 = dcenterCoefficients[0];
        double[] df1N = dcenterCoefficients[1];


        Vector v = box.getSpace().makeVector();
        Vector dr = box.getSpace().makeVector();
        Vector v0 = box.getSpace().makeVector();
        Vector vPrev = box.getSpace().makeVector();
        double En = 0;
        int ns = nShifts+1;
        for (int i = 0; i < molecules.size(); i++) {
            IAtomList beads = molecules.get(i).getChildList();
            if (ns > beads.size() || (beads.size() / ns) * ns != beads.size()) {
                throw new RuntimeException("# of beads must be a multiple of (# of shifts + 1)");
            }
            for (int indexShift=0; indexShift<beads.size(); indexShift += beads.size()/ns) {
                Vector rPrev = beads.get(indexShift).getPosition();
                Vector r0 = rPrev;

                for (int j = 0; j < beads.size(); j++) {
                    int aj = (j + indexShift) % beads.size();
                    IAtom atomj = beads.get(aj);
                    En -= gamma[j]; // Jacobian
                    Vector rj = atomj.getPosition();
                    dr.E(rj);
                    if (j > 0) {
                        dr.PEa1Tv1(-f11[j], rPrev);
                        dr.PEa1Tv1(-f1N[j], r0);
                    }
                    v.Ea1Tv1(gamma[j], dr);
                    if (j > 0) {
                        v.PEa1Tv1(df11[j], rPrev);
                        v.PEa1Tv1(df1N[j], r0);
                        v.PEa1Tv1(f11[j], vPrev);
                        v.PEa1Tv1(f1N[j], v0);
                    }
                    int jj = atomj.getLeafIndex();
                    En -= beta * forcesU[jj].dot(v);
                    if (nBeads > 1) {
                        En -= beta * forcesK[jj].dot(v);
                    }
                    rPrev = rj;
                    if (j == 0) v0.E(v);
                    vPrev.E(v);
                }
            }
        }
//        System.out.println("RealHMA: "+En);

        x[0] = En0 + En/ns - EnShift;

        return data;
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