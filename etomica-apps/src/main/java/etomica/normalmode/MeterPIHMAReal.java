/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

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

    @Override
    public IData getData() {
        Box box = move.getBox();
        double[] x = data.getData();
//        System.out.println("******** HMA *************");

        int nBeads = box.getLeafList().size();
        double En = 0.5*nBeads/beta + pcP1.computeAll(true) - pmBonding.computeAll(true);

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
        for (int i = 0; i<molecules.size(); i++) {
            IAtomList beads = molecules.get(i).getChildList();
            Vector rPrev = beads.get(0).getPosition();
            Vector r0 = rPrev;

            for (int j=0; j<beads.size(); j++) {
                int N = beads.size() - j;
                En += dSigma[j] / sigma[j]; // Jacobian

                dr.E(beads.get(j).getPosition());
                if (j>0) {
                    dr.PEa1Tv1(-R11[N], rPrev);
                    dr.PEa1Tv1(-R1N[N], r0);
                }
                v.Ea1Tv1(dSigma[j]/sigma[j], dr);
                if (j>0) {
                    v.PEa1Tv1(dR11[N], rPrev);
                    v.PEa1Tv1(dR1N[N], r0);
                    v.PEa1Tv1(R11[N], vPrev);
                    v.PEa1Tv1(R1N[N], v0);
                }
                En -= beta*(forcesU[i].dot(v));
                if (nBeads > 1) {
                    En -= beta*(forcesK[i].dot(v));
                }
                rPrev = beads.get(j).getPosition();
                if (j==0) v0.E(v);
                vPrev.E(v);
            }
        }

        x[0] = En - EnShift;

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
