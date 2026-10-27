/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidHMAFD implements IDataSource {

    protected Vector[] latticePositions;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected double V, Vm, Vp;
    protected double rho;
    protected int N, nCells;
    protected double gV, gVV, mV;
    protected Vector dri;
    protected final double T;
    protected PotentialMasterList potentialMaster0, potentialMasterP, potentialMasterM;
    protected Vector[] rOrig;
    protected double dbeta, deV;
    protected double scaleVm, scaleVp, scaleBetaM, scaleBetaP;
    protected double beta, betaP, betaM;
    protected Box box0, boxP, boxM;
    protected double uShift, pShift;
    protected double U_0;

    public MeterSolidHMAFD(Space space, PotentialMasterList potentialMaster0, PotentialMasterList potentialMasterP, PotentialMasterList potentialMasterM, CoordinateDefinition coordinateDefinition, double T, double[] elasticParams, double dbeta, double deV) {
        int nData = 8;
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.potentialMaster0 = potentialMaster0;
        this.potentialMasterP = potentialMasterP;
        this.potentialMasterM = potentialMasterM;
        this.T = T;

        this.box0 = potentialMaster0.getNeighborManager().getBox();
        this.boxP = potentialMasterP.getNeighborManager().getBox();
        this.boxM = potentialMasterM.getNeighborManager().getBox();

        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        // dbeta
        beta = 1 / this.T;
        this.dbeta = dbeta;
        betaP = beta + dbeta;
        betaM = beta - dbeta;
        scaleBetaP = Math.sqrt(beta / betaP);
        scaleBetaM = Math.sqrt(beta / betaM);

        //dV
        this.deV = deV;
        V = box0.getBoundary().volume();
        N = box0.getLeafList().size();
        rho = N / V;
        Vm = V * (1 - deV);
        Vp = V * (1 + deV);
        scaleVp = Math.pow(Vp / V, 1.0 / 3.0);
        scaleVm = Math.pow(Vm / V, 1.0 / 3.0);

        gV   = elasticParams[0];
        gVV  = elasticParams[1];
        mV = gV - 1.0 / 3.0;

        dri = space.makeVector();
        rOrig = box0.getSpace().makeVectorArray(box0.getLeafList().size());

        latticePositions = box0.getSpace().makeVectorArray(N);
        for (int i = 0; i < N; i++) {
            latticePositions[i].E(coordinateDefinition.getLatticePosition(box0.getLeafList().get(i)));
        }

        this.uShift = 0;
        this.pShift = 0;
    }


    public IData getData() {
        double[] x = data.getData();

        for (IAtom atom : box0.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            rOrig[iAtom].E(atom.getPosition());
        }

        U_0 = potentialMaster0.computeAll(false);

        double[] dbUdbm = dbUdbm();
        double[] dUdeVc = dUdeVc();
        double[] dUdeVm = dUdeVm();

        //U
        x[0] = U_0 / N - uShift;
        x[1] = 3.0 * (N - 1.0) / 2.0 / N / beta  + dbUdbm[0] / N - uShift;

        //P
        x[2] = rho * T - dUdeVc[0] / V - pShift;
        x[3] = 3.0 * (N - 1.0) / V * gV * T - dUdeVm[0] / V + T / V  - pShift;

        //Cv
        x[4] = 3.0 / 2.0 * (N - 1.0) / N - dbUdbm[1] / N / T / T;

        //B
        x[5] = dUdeVc[1] / V + rho * T;
        x[6] = dUdeVm[1] / V - 3.0 * (N - 1.0) / V * T * gVV  + T / V;

        // G
        x[7] = dUdeVm[2] / V / T;

        return data;
    }


    public double[] dbUdbm() {
        // +dbeta
        for (IAtom atom : box0.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaP);
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }
        double U_p = potentialMaster0.computeAll(false);

        // -dbeta
        for (IAtom atom : box0.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaM / scaleBetaP);
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }
        double U_m = potentialMaster0.computeAll(false);

        for (IAtom atom : box0.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().E(rOrig[iAtom]);
        }

        double dbUdbm = (betaP * U_p - betaM * U_m) / 2.0 / dbeta;
        double d2bUdb2 = (betaP * U_p - 2.0 * beta * U_0 + betaM * U_m) / dbeta / dbeta;

        return new double[]{dbUdbm, d2bUdb2};
    }


    public double[] dUdeVc() {
        // +dV
        for (IAtom atom : boxP.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVp, box0.getLeafList().get(iAtom).getPosition());
        }
        double U_p = potentialMasterP.computeAll(false);

        for (IAtom atom : boxM.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVm, box0.getLeafList().get(iAtom).getPosition());
        }
        double U_m = potentialMasterM.computeAll(false);

        double dUdeV = (U_p - U_m) / 2.0 / deV;
        double d2UdeV2 = (U_p - 2.0 * U_0 +  U_m) / deV / deV;

        return new double[] {dUdeV, d2UdeV2};
    }


    public double[] dUdeVm() {
        // +deV
        for (IAtom atom : boxP.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVp, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(Math.exp(gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }
        double U_p = potentialMasterP.computeAll(false);

        // -deV
        for (IAtom atom : boxM.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVm, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(Math.exp(-gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }

        double U_m = potentialMasterM.computeAll(false);

        double dUdeV = (U_p - U_m) / 2.0 / deV;
        double d2UdeV2 = (U_p - 2.0 * U_0 + U_m) / deV / deV;

        //Cross derivatives for gamma_V

        //+deV +dbeta
        for (IAtom atom : boxP.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVp, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaP * Math.exp(gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }
        double U_pp = potentialMasterP.computeAll(false);

        //+deV -dbeta
        for (IAtom atom : boxP.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVp, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaM * Math.exp(gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }
        double U_pm = potentialMasterP.computeAll(false);

        //-deV +dbeta
        for (IAtom atom : boxM.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVm, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaP * Math.exp(-gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }
        double U_mp = potentialMasterM.computeAll(false);

        //-deV -dbeta
        for (IAtom atom : boxM.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().Ea1Tv1(scaleVm, latticePositions[iAtom]);
            dri.Ev1Mv2(box0.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(scaleBetaM * Math.exp(-gV * deV + 1.0 / 2.0 * gVV * deV * deV));
            atom.getPosition().PE(dri);
        }
        double U_mm = potentialMasterM.computeAll(false);

        double d2UdeVdbeta = (betaP * U_pp - betaM * U_pm - betaP * U_mp + betaM * U_mm) / 4.0 / deV / dbeta;

        return new double[] {dUdeV, d2UdeV2, d2UdeVdbeta};
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public void setShift(double uShift, double pShift) {
        this.uShift = uShift;
        this.pShift = pShift;
    }


}