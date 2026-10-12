/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
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
    protected final DataSourceScalar meterPE;
    protected double V, Vm, Vp;
    protected double rho;
    protected int N, nCells;
    protected double gV, gVV, mV;
    protected double dim;
    protected Vector dri, drj, drij;
    protected final double temperature;
    protected PotentialMasterList potentialMaster, potentialMasterp, potentialMasterm;
    protected Vector[] rOrig;
    protected int n;
    protected double dbeta, deV;
    protected double scaleVm, scaleVp, scaleTm, scaleTp;
    protected double beta, betap, betam;
    protected Box box, boxp, boxm;

    public MeterSolidHMAFD(Space space, DataSourceScalar meterPE, PotentialMasterList potentialMaster, PotentialMasterList potentialMasterp, PotentialMasterList potentialMasterm, Box boxp, Box boxm, CoordinateDefinition coordinateDefinition, double temperature, double[] elasticParams, double dbeta, double deV) {
        int nData = 7;
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        this.potentialMasterp = potentialMasterp;
        this.potentialMasterm = potentialMasterm;
        this.temperature = temperature;
        this.boxp = boxp;
        this.boxm = boxm;
        dim = space.D();
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        // dbeta
        this.dbeta = dbeta;
        beta = 1 / this.temperature;
        betap = beta + dbeta;
        betam = beta - dbeta;
        scaleTp = Math.sqrt(beta / betap);
        scaleTm = Math.sqrt(beta / betam);

        //dV
        this.deV = deV;
        box = coordinateDefinition.getBox();
        V = box.getBoundary().volume();
        N = box.getLeafList().size();
        rho = N / V;
        Vm = V * (1 - deV);
        Vp = V * (1 + deV);
        scaleVp = Math.pow(Vp / V, 1.0 / dim);
        scaleVm = Math.pow(Vm / V, 1.0 / dim);

        gV   = elasticParams[0];
        gVV  = elasticParams[1];
        mV = gV - 1.0 / dim;

        dri = space.makeVector();
        drj = space.makeVector();
        drij = space.makeVector();
        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());

        latticePositions = box.getSpace().makeVectorArray(N);
        for (int i = 0; i < N; i++) {
            latticePositions[i].E(coordinateDefinition.getLatticePosition(box.getLeafList().get(i)));
        }

        n = 0;
    }


    public IData getData() {
        n++;
        double[] x = data.getData();

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            rOrig[iAtom].E(atom.getPosition());
        }

        double Fdr = 0;
        Box box = coordinateDefinition.getBox();
        potentialMaster.computeAll(true);

        double uSum = potentialMaster.getLastEnergy();
        double virial = potentialMaster.getLastVirial();
        Vector[] forces = potentialMaster.getForces();

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            Vector F = forces[iAtom];
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            box.getBoundary().nearestImage(dri);
            Fdr += F.dot(dri);
        }

        // Analytical derivatives
        x[0] = uSum / N; // U_conv
        x[1] = uSum / N + dim * (N - 1.0) / 2.0 / N * temperature  + 0.5 * Fdr / N; // U_hma

        x[2] = rho * temperature - virial / (3 * V); //P_conv
        x[3] = dim * (N - 1.0) / V * temperature * gV - virial / (3 * V) + mV * Fdr / V; // P_conv

        //Finite-difference derivatives
        x[4] = dim * (N - 1.0) / 2.0 / N * temperature + dbUdb() / N; //U_hma
        x[5] = rho * temperature - dUdeVc() / V; // P_conv
        x[6] = dim * (N - 1.0) / V * gV * temperature - dUdeVm() / V; // P_hma


//        System.out.println(n + "    " + x[2] + "    " + x[3] + "    " +  x[5] + "    " + x[6]);

        return data;
    }



    /**
     * dU / dbeta
     * @return the energy derivative w.r.t. beta
     */
    public double dbUdb() {
        // +dbeta
        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            dri.TE(scaleTp);
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }

        double U_p = potentialMaster.computeAll(false);

        // -dbeta
        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            dri.Ev1Mv2(atom.getPosition(), latticePositions[iAtom]);
            dri.TE(scaleTm / scaleTp); //undo first
            atom.getPosition().Ev1Pv2(latticePositions[atom.getLeafIndex()], dri);
        }

        double U_m = potentialMaster.computeAll(false);

        double dbUdb = (betap * U_p - betam * U_m) / 2.0 / dbeta;

        for (IAtom atom : box.getLeafList()) {
            int iAtom = atom.getLeafIndex();
            atom.getPosition().E(rOrig[iAtom]);
        }

        return dbUdb;
    }



    /**
     * dU / deV
     * @return the energy derivative w.r.t. V
     */
    public double dUdeVc() {
        // +dV
        for (IAtom atomp : boxp.getLeafList()) {
            int iAtom = atomp.getLeafIndex();
            atomp.getPosition().Ea1Tv1(scaleVp, box.getLeafList().get(iAtom).getPosition());
        }
        double U_p = potentialMasterp.computeAll(false);

        for (IAtom atomm : boxm.getLeafList()) {
            int iAtom = atomm.getLeafIndex();
            atomm.getPosition().Ea1Tv1(scaleVm, box.getLeafList().get(iAtom).getPosition());
        }
        double U_m = potentialMasterm.computeAll(false);

        double dUdeV = (U_p - U_m) / 2.0 / deV;

        return dUdeV;
    }

    // -deV
    public double dUdeVm() {
        for (IAtom atomp : boxp.getLeafList()) {
            int iAtom = atomp.getLeafIndex();
            atomp.getPosition().Ea1Tv1(scaleVp, latticePositions[iAtom]);
            dri.Ev1Mv2(box.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(Math.exp(gV * deV));
            atomp.getPosition().PE(dri);
        }
        double U_p = potentialMasterp.computeAll(false);

        for (IAtom atomm : boxm.getLeafList()) {
            int iAtom = atomm.getLeafIndex();
            atomm.getPosition().Ea1Tv1(scaleVm, latticePositions[iAtom]);
            dri.Ev1Mv2(box.getLeafList().get(iAtom).getPosition(), latticePositions[iAtom]);
            dri.TE(Math.exp(-gV * deV));
            atomm.getPosition().PE(dri);
        }
        double U_m = potentialMasterm.computeAll(false);

        double dUdeV = (U_p - U_m) / 2.0 / deV;

        return dUdeV;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

}