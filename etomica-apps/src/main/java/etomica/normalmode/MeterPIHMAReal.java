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

public class MeterPIHMAReal implements IDataSource {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialComputeField pcP1;
    protected double beta;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected double EnShift;
    protected final MCMoveHOReal move;
    protected int nShifts = 0;

    public MeterPIHMAReal(PotentialMasterBonding pmBonding, PotentialComputeField pcP1, double beta, MCMoveHOReal move) {
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

        double[] sigma = move.getChainSigmas();
        double[] dSigma = move.getDChainSigmas();
        double[][] centerCoefficients = move.getCenterCoefficients();
        double[] R11 = centerCoefficients[0];
        double[] R1N = centerCoefficients[1];
        double[][] dcenterCoefficients = move.getDCenterCoefficients();
        double[] dR11 = dcenterCoefficients[0];
        double[] dR1N = dcenterCoefficients[1];

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
                    int N = beads.size() - j;
                    En -= dSigma[j] / sigma[j]; // Jacobian
    //                System.out.println("realJ: "+j+" "+dSigma[j]/sigma[j]);

                    Vector rj = atomj.getPosition();
                    dr.E(rj);
                    if (j > 0) {
                        dr.PEa1Tv1(-R11[N], rPrev);
                        dr.PEa1Tv1(-R1N[N], r0);
                    }
                    if (j == 0) {
                        v.Ea1Tv1(dSigma[0] / sigma[0], dr);
                    } else {
                        v.Ea1Tv1(dSigma[N] / sigma[N], dr);
                    }
                    if (j > 0) {
                        v.PEa1Tv1(dR11[N], rPrev);
                        v.PEa1Tv1(dR1N[N], r0);
                        v.PEa1Tv1(R11[N], vPrev);
                        v.PEa1Tv1(R1N[N], v0);
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
