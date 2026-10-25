/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSourcePotential;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidHMALJ implements IDataSourcePotential, PotentialCallback {
    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialCompute potentialMaster;
    protected double T;
    protected double latticeEnergy, latticePressure;
    protected final Box box;
    protected final CoordinateDefinition coordinteDefinition;
    protected final Space space;
    protected double uShift, pShift;
    protected double gV, gVV, mV, mVV;
    protected double fV , fVV, dP, dB;
    protected double gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
    protected double mx1, my1, my4, mx11, my11, mx12, mz12, mx44, my44;
    protected double mT, mVT, mx1T, my1T;
    protected Vector dri, drj, drij;
    protected boolean callComputeAll = true;
    protected Vector[] dFdeV;
    protected double Fr , Fdr , rHr , drHdr, rHdr;

    public MeterSolidHMALJ(Space space, PotentialCompute potentialCompute, CoordinateDefinition coordinateDefinition, double[] elasticParams, double T, boolean doD2) {
        int nData = 8;
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
        this.potentialMaster = potentialCompute;
        dim = space.D();
        this.space = space;
        this.T = T;
        box = coordinateDefinition.getBox();
        latticeEnergy = potentialCompute.computeAll(false);
        latticePressure = -potentialCompute.getLastVirial() / (box.getBoundary().volume() * dim);

        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{nData});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(nData);

        gV = elasticParams[0];
        gVV = elasticParams[1];

        int numAtoms = box.getMoleculeList().size();
        double V = box.getBoundary().volume();
        double rho = numAtoms / V;
        fV = (dP/T-rho)/(3*(numAtoms-1));
        fVV = fV*fV - (dB/T-rho)/(3*V*(numAtoms-1)) + 2.0*fV/(3.0*V);

        mV = gV - 1.0 / 3.0;
        mVV = gVV + gV * gV + 2.0 / 9.0;
        mT = 1.0 /2.0/T;
        mVT = mT * (mV + 1.0 / 3.0);
        mx1 = gx1 - 1.0;
        my1 = gy1;
        my4 = gy4 - 0.5;
        mx11 = gx11 + gx1 * gx1 + 1.0;
        my11 = gy11 + gy1 * gy1;
        mx12 = gx12 + gx1 * gy1; //gx2=gy1
        mz12 = gz12 + gy1 * gy1; //gz1=gy1 , gz2=gy1
        mx44 = gx44; //gx4=0
        my44 = gy44 + gy4 * gy4 + 1.0 / 4.0;
        mx1T = mT * gx1;
        my1T = mT * gy1;

        this.uShift = 0;
        this.pShift = 0;

        dri = space.makeVector();
        drj = space.makeVector();
        drij = space.makeVector();

        dFdeV = new Vector[0];

    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    public void setTemperature(double T) {
        this.T = T;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    @Override
    public IData getData() {
        Fr = 0 ;
        Fdr = 0;
        rHr = 0 ;
        drHdr = 0;
        rHdr = 0;

        double uSum;
        if (callComputeAll) {
            PotentialCallback callback = this;
            uSum = potentialMaster.computeAll(false, callback);
        } else {
            uSum = potentialMaster.getLastEnergy();
        }

        double[] x = data.getData();
        double V = box.getBoundary().volume();
        int N = box.getMoleculeList().size();
        double rho = N / V;

        x[0] = uSum / N - uShift;
        x[1] = uSum / N + 3.0 / 2.0 * (N - 1.0) / N * T + 1.0 / 2.0 / N * Fdr - uShift;

        //Pc and Pm
        x[2] = rho * T + Fr / (3.0 * V) - pShift;
        x[3] = 3.0 * (N - 1.0) / V * T * gV + Fr / (3.0 * V) + mV * Fdr / V + T / V - pShift;

        //Second derivatives
        //Cvm
        x[4] = 3.0 / 2.0 * (N - 1.0) / N  - 1.0 / 4.0 / T * (Fdr + drHdr) / N;

        // B
        x[5] = rho * T + (1.0 / 9.0 * rHr + 2.0 / 9.0 * Fr) / V;
        x[6] = - 3.0 * (N - 1.0) * T * gVV / V + (1.0 / 9.0 * rHr + 2.0 / 9.0 * Fr +  mV * mV * drHdr - mVV * Fdr + 2.0 / 3.0 * mV * rHdr) / V + T / V;

        //Gm
        x[7] = 3.0 * (N - 1.0) * gV / V + (gV * Fdr - mV * drHdr - 1.0 / 3.0 * rHdr) / 2.0 / V / T + 1.0 / V;

        return data;
    }

    @Override
    //pairCompute called by potentialcompute
    public void pairCompute(int i, int j, Vector rij, double[] u012) { //rij = rj - ri
//        if (i < j) System.out.println("Hessian: " + i + "  " + j);

        double du = u012[1]; // r u'
        double d2u = u012[2];// r2 u''
        double rij2 = rij.squared();
        IAtom iAtom = box.getLeafList().get(i);
        IAtom jAtom = box.getLeafList().get(j);
        dri.Ev1Mv2(iAtom.getPosition(), coordinteDefinition.getLatticePosition(iAtom));
        box.getBoundary().nearestImage(dri);
        drj.Ev1Mv2(jAtom.getPosition(), coordinteDefinition.getLatticePosition(jAtom));
        box.getBoundary().nearestImage(drj);
        drij.Ev1Mv2(drj, dri); // drij = rj - ri  (NOT ri - rj)!!!

        // Fr
        Fr += -du;
        // Fdr
        Fdr += -du / rij2 * rij.dot(drij);

        // Hij
        Tensor Hij = space.makeTensor();
        Hij.Ev1v2(rij, rij);
        Hij.TE((d2u - du) / rij2 / rij2);
        Hij.PE(0,0, du / rij2);
        Hij.PE(1,1, du / rij2);
        Hij.PE(2,2, du / rij2);

        Vector tmpV = space.makeVector();

        // rHr
        tmpV.E(rij);
        Hij.transform(tmpV);
        rHr += rij.dot(tmpV);

        // drHdr
        tmpV.E(drij);
        Hij.transform(tmpV); // tmpV = Hij.drj
        drHdr += drij.dot(tmpV);

        // rHdr
        tmpV.E(drij);
        Hij.transform(tmpV); // tmpV = Hij.drj
        rHdr += rij.dot(tmpV);
    }

    public void setShift(double uShift, double pShift) {
        this.uShift = uShift;
        this.pShift = pShift;
    }

    @Override
    public boolean needsForces() {
        return true;
    }

    @Override
    public PotentialCallback getPotentialCallback() {
        return this;
    }

    @Override
    public void doCallComputeAll(boolean callComputeAll) {
        this.callComputeAll = callComputeAll;
    }

    public boolean wantsHessian() {
        return true;
    }

    public boolean needsPairCallback() {
//        return doD2;
        return true;
    }

}