/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.BoxInflate;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
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

public class MeterSolidHMA implements IDataSourcePotential, PotentialCallback {
    protected final int dim;
    protected final DataTag tag;
    protected DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final PotentialCompute potentialMaster;
    protected double temperature;
    protected double latticeEnergy, latticePressure;
    protected final Box box;
    protected double pRes, bpHarm;
    protected final CoordinateDefinition coordinteDefinition;
    protected double virial2, drHdr, rHdr;
    protected final Space space;
    protected double uShift, pShift;

    protected double virialx, virialy, virialz, virialxy, virialxz, virialyz;
    protected double x_Hxx_x, y_Hyy_y, z_Hzz_z;
    protected double x_Hxy_y, x_Hxz_z, y_Hyz_z;
    protected double x_Hyx_y, x_Hzx_z, y_Hzy_z;

    protected double dr1_H_dr1, dr2_H_dr2, dr3_H_dr3, dr1_H_dr2, dr1_H_dr3, dr2_H_dr3;
    protected double dr4_H_dr4, dr5_H_dr5, dr6_H_dr6, rx1_H_dr1, ry2_H_dr2, rz3_H_dr3;
    protected double rx1_H_dr2, rx1_H_dr3, ry2_H_dr3, ry2_H_dr1, rz3_H_dr1, rz3_H_dr2;
    protected double ry3z2_H_dr4, rx3z1_H_dr5, rx2y1_H_dr6;
    protected double drHdr1, drHdr2, drHdr3, rx1_H_dr, ry2_H_dr, rz3_H_dr;

    protected double gV, gVV, mV, mVV;


    protected double fV , fVV, dP, dB;

    protected double gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
    protected double mx1, my1, my4, mx11, my11, mx12, mz12, mx44, my44;
    protected double mT, mVT, mx1T, my1T;
    protected final boolean doD2;
    protected Vector dri, drj, drij;

    protected boolean callComputeAll = true;

    public MeterSolidHMA(Space space, PotentialCompute potentialCompute, CoordinateDefinition coordinateDefinition, double[] elasticParams, double temperature, boolean doD2) {
        this.doD2 = doD2;
        this.coordinteDefinition = coordinateDefinition;
        tag = new DataTag();
        this.potentialMaster = potentialCompute;
        dim = space.D();
        this.space = space;
        this.temperature = temperature;
        box = coordinateDefinition.getBox();
        latticeEnergy = potentialCompute.computeAll(false);
        latticePressure = -potentialCompute.getLastVirial() / (box.getBoundary().volume() * dim);
//        int n = doD2 ? 41 : 16;
        int n = doD2 ? 19 : 16;

        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{n});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(n);

        gV = elasticParams[0];
        gVV = elasticParams[1];
        gx1 = elasticParams[2];
        gy1 = elasticParams[3];
        gy4 = elasticParams[4];
        gx11 = elasticParams[5];
        gy11 = elasticParams[6];
        gx44 = elasticParams[7];
        gy44 = elasticParams[8];
        gx12 = elasticParams[9];
        gz12 = elasticParams[10];

        int numAtoms = box.getMoleculeList().size();
        double V = box.getBoundary().volume();
        double rho = numAtoms / V;

        mV = gV - 1.0 / 3.0;
        mVV = gVV + gV * gV + 2.0 / 9.0;

        // 256, 100K, rho=0.1
        dP = 39.73990629584506; //324.771547734307500-285.03164143846243
        dB = -4.16; //1473.43-1477.59;

        fV = (dP/temperature-rho)/(3*(numAtoms-1));
        fVV = fV*fV - 1.0/V*(dB/temperature-rho)/(3*(numAtoms-1));
        mVV = mV*mV - V*(dB/temperature-rho)/(3*(numAtoms-1));
//        gV = mV + 1/3.0;
        gVV = mVV - gV*gV -2.0/3.0;

        mT = 1.0 / (2.0 * temperature);
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
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

    public void setPRes(double pRes) {
        this.pRes = pRes;
        bpHarm = pRes / temperature;
    }

    /**
     * Computes total pressure in box by summing virial over all pairs, and adding
     * ideal-gas contribution.
     */
    @Override
    public IData getData() {
        double Fdr = 0, Frdot1 = 0, Frdot2 = 0, Frdot3 = 0, Frdot4 = 0, Frdot5 = 0, Frdot6 = 0;
        double Frddot11 = 0, Frddot22 = 0, Frddot33 = 0, Frddot12 = 0, Frddot13 = 0, Frddot23 = 0, Frddot44 = 0, Frddot55 = 0, Frddot66 = 0;
        double Frddot1T = 0, Frddot2T = 0, Frddot3T = 0;
        virialx = virialy = virialz = virialxy = virialxz = virialyz = 0;
        virial2 = drHdr = rHdr = 0;
        x_Hxx_x = y_Hyy_y = z_Hzz_z = x_Hxy_y = x_Hxz_z = y_Hyz_z = x_Hyx_y = x_Hzx_z = y_Hzy_z = 0;
        dr1_H_dr1 = dr2_H_dr2 = dr3_H_dr3 = dr1_H_dr2 = dr1_H_dr3 = dr2_H_dr3 = 0;
        dr4_H_dr4 = dr5_H_dr5 = dr6_H_dr6 = rx1_H_dr1 = ry2_H_dr2 = rz3_H_dr3 = 0;
        drHdr1 = drHdr2 = drHdr3 = rx1_H_dr = ry2_H_dr = rz3_H_dr = 0;
        rx1_H_dr2 = rx1_H_dr3 = ry2_H_dr3 = ry2_H_dr1 = rz3_H_dr1 = rz3_H_dr2 = 0;
        ry3z2_H_dr4 = rx3z1_H_dr5 = rx2y1_H_dr6 = 0;

        double uSum;
        if (callComputeAll) {
            PotentialCallback callback = doD2 ? this : null;
            System.out.println(" zero: V = " + box.getBoundary().volume());

            uSum = potentialMaster.computeAll(true, callback);
        } else {
            uSum = potentialMaster.getLastEnergy();
        }

        double[] x = data.getData();
        double V = box.getBoundary().volume();
        int N = box.getMoleculeList().size();
        double rho = N / V;
        double virial = potentialMaster.getLastVirial();
        double virial2 = potentialMaster.getLastVirial2();
        Vector[] forces = potentialMaster.getForces();
        Vector[] dFdeV = potentialMaster.getdFdeV();
        double dFdeVdr = 0;

//        if(false) {
        if(wantsHessian() && dFdeV != null){

            if(doD2){
                int nn = 0;
//                System.out.println(" dFdeV: " + dFdeV[nn].getX(0));

                double V0 = V;
                double deV = 0.000001;
                double Vp = V0*(1+deV);
                double rhop = box.getLeafList().size()/Vp;
                double Vm = V0*(1-deV);
                double rhom = box.getLeafList().size()/Vm;

                BoxInflate inflater = new BoxInflate(box, space);
                inflater.setTargetDensity(rhom);
                inflater.actionPerformed();

                System.out.println(" minus: V = " + Vm);
                potentialMaster.computeAll(true, this);
                double forces_m = potentialMaster.getForces()[nn].getX(0);

                System.out.println(" plus: V = " + Vm);
                inflater.setTargetDensity(rhop);
                inflater.actionPerformed();
                potentialMaster.computeAll(true, this);
                double forces_p = potentialMaster.getForces()[nn].getX(0);
            }
            System.exit(0);


        }



        IAtomList atoms = box.getLeafList();
        for (IAtom a : atoms) {
            Vector F = forces[a.getLeafIndex()];
            dri.Ev1Mv2(a.getPosition(), coordinteDefinition.getLatticePosition(a));
            box.getBoundary().nearestImage(dri);

            if (doD2) dFdeVdr += dFdeV[a.getLeafIndex()].dot(dri);


            Vector[] rdot = mapVel(dri);
            Vector[] rddot = mapAcc(dri);

            // F.rdot
            Fdr += F.dot(dri);
            Frdot1 += F.dot(rdot[0]);
            Frdot2 += F.dot(rdot[1]);
            Frdot3 += F.dot(rdot[2]);
            Frdot4 += F.dot(rdot[3]);
            Frdot5 += F.dot(rdot[4]);
            Frdot6 += F.dot(rdot[5]);

            // F.rddot
            Frddot11 += F.dot(rddot[0]);
            Frddot22 += F.dot(rddot[1]);
            Frddot33 += F.dot(rddot[2]);
            Frddot44 += F.dot(rddot[3]);
            Frddot55 += F.dot(rddot[4]);
            Frddot66 += F.dot(rddot[5]);
            Frddot23 += F.dot(rddot[6]);
            Frddot13 += F.dot(rddot[7]);
            Frddot12 += F.dot(rddot[8]);

            //b_ii
            Frddot1T += F.dot(rddot[9]);
            Frddot2T += F.dot(rddot[10]);
            Frddot3T += F.dot(rddot[11]);
        }

        //Props
        //U
        x[0] = uSum / N - uShift;

        x[1] = uSum / N + 3.0 * (N - 1.0) / N * temperature / 2.0 + 1.0 / 2.0 / N * Fdr - uShift;
        //P
        x[2] = rho * temperature - virial / (3 * V) - pShift;
//        x[3] = 3 * (N - 1) / V * temperature * gV - virial / (3 * V) + mV * Fdr / V - pShift;
//        gV = dP;

        x[3] = dP - virial / (3 * V) + fV * Fdr - pShift;


        // V*(dP/temperature-rho)/(dim*(numAtoms-1)) + 1/3.0;

        double p1Shift = -pShift;
        //P1
        x[4] = -rho * temperature + virialx / V - p1Shift;
        x[5] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialx + Frdot1) / V - p1Shift;
        //P2
        x[6] = -rho * temperature + virialy / V - p1Shift;
        x[7] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialy + Frdot2) / V - p1Shift;
        //P3
        x[8] = -rho * temperature + virialz / V - p1Shift;
        x[9] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialz + Frdot3) / V - p1Shift;
        //P4
        x[10] = virialyz / V;
        x[11] = (virialyz - Frdot4) / V;
        //P5
        x[12] = virialxz / V;
        x[13] = (virialxz - Frdot5) / V;
        //P6
        x[14] = virialxy / V;
        x[15] = (virialxy - Frdot6) / V;


        //Second derivatives
        if (doD2) {
            //Cv - HMA
            x[16] = -1.0/4.0/temperature*(Fdr + drHdr);

            // B
            x[17] = virial2/V + rho * temperature;
//            x[18] = (virial2 + mV * mV * drHdr - mVV * Fdr + 2.0 / 3.0 * mV * rHdr - 3.0 * (N - 1) * temperature * gVV) / V + temperature / V;
//            x[18] = (virial2 + mV * mV * drHdr - mVV * Fdr - 2.0* mV * dFdeVdr - 3.0 * (N - 1) * temperature * gVV) / V + temperature / V;
            x[18] = dB + virial2/V - V*fVV * Fdr - 2.0*fV*dFdeVdr  + V*fV*fV*drHdr ;
//
//            //C11
//            x[19] = x_Hxx_x / V + 2 * rho * temperature; //C11
//            x[20] = (x_Hxx_x - Frddot11 + dr1_H_dr1 + 2 * rx1_H_dr1 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
//            //C22
//            x[21] = y_Hyy_y / V + 2 * rho * temperature; //C11
//            x[22] = (y_Hyy_y - Frddot22 + dr2_H_dr2 + 2 * ry2_H_dr2 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
//            //C33
//            x[23] = z_Hzz_z / V + 2 * rho * temperature; //C11
//            x[24] = (z_Hzz_z - Frddot33 + dr3_H_dr3 + 2 * rz3_H_dr3 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
//
//            //C12
//            x[25] = x_Hxy_y / V;//C12
//            x[26] = (x_Hxy_y - Frddot12 + dr1_H_dr2 + rx1_H_dr2 + ry2_H_dr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
//            //C13
//            x[27] = x_Hxz_z / V;//C12
//            x[28] = (x_Hxz_z - Frddot13 + dr1_H_dr3 + rx1_H_dr3 + rz3_H_dr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
//            //C23
//            x[29] = y_Hyz_z / V;//C12
//            x[30] = (y_Hyz_z - Frddot23 + dr2_H_dr2 + ry2_H_dr3 + rz3_H_dr2 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
//
//            //C44
//            x[31] = y_Hzy_z / V + rho * temperature;
//            x[32] = (y_Hzy_z - Frddot44 + dr4_H_dr4 + ry3z2_H_dr4 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//            //C55
//            x[33] = x_Hzx_z / V + rho * temperature;//C55
//            x[34] = (x_Hzx_z - Frddot55 + dr5_H_dr5 + rx3z1_H_dr5 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//            //C6
//            x[35] = x_Hyx_y / V + rho * temperature;//C66
//            x[36] = (x_Hyx_y - Frddot66 + dr6_H_dr6 + rx2y1_H_dr6 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//
//            //b_V
//            x[37] = 1.0 / 2.0 / V / temperature * (gV * Fdr - mV * drHdr - 1.0 / 3.0 * rHdr) + 3.0 * (N - 1) * gV / V + 1.0 / V; //bV_hma
//            //b_mn
//            x[38] = 1.0 / V * (-Frddot1T + mT * drHdr1 + mT * rx1_H_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
//            x[39] = 1.0 / V * (-Frddot2T + mT * drHdr2 + mT * ry2_H_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
//            x[40] = 1.0 / V * (-Frddot3T + mT * drHdr3 + mT * rz3_H_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
        }

        return data;
    }

    @Override
    //pairCompute called by potentialcompute
    public void pairCompute(int iAtom, int jAtom, Vector rij, double[] u012) { //rij = rj - ri
        double du = u012[1];
        double d2u = u012[2];
        double rij2 = rij.squared();

        IAtom ia = box.getLeafList().get(iAtom);
        IAtom ja = box.getLeafList().get(jAtom);

        dri.Ev1Mv2(ia.getPosition(), coordinteDefinition.getLatticePosition(ia));
        box.getBoundary().nearestImage(dri);
        drj.Ev1Mv2(ja.getPosition(), coordinteDefinition.getLatticePosition(ja));
        box.getBoundary().nearestImage(drj);
        drij.Ev1Mv2(drj, dri);

        virialx += du / rij2 * rij.getX(0) * rij.getX(0);
        virialy += du / rij2 * rij.getX(1) * rij.getX(1);
        virialz += du / rij2 * rij.getX(2) * rij.getX(2);
        virialyz += du / rij2 * rij.getX(1) * rij.getX(2);
        virialxz += du / rij2 * rij.getX(0) * rij.getX(2);
        virialxy += du / rij2 * rij.getX(0) * rij.getX(1);


        if(doD2){
            virial2 += r1Hr2(rij, rij, rij, du, d2u);

//            r1Hr2(Vector rij, Vector r1ij, Vector r2ij, double du, double d2u) { //r1.H.r2
//            double rij2 = rij.squared();
//            double virial2 = du / rij2 * r1ij.dot(r2ij) + (d2u - du) / rij2 / rij2 * r1ij.dot(rij) * (r2ij.dot(rij));
//            return virial2;

            drHdr += r1Hr2(rij, drij, drij, du, d2u);
//            rHdr += r1Hr2(rij, rij, drij, du, d2u);

//            x_Hxx_x += r1xH2abr2y(rij, 0, 0, rij.getX(0), rij.getX(0), du, d2u);
//            y_Hyy_y += r1xH2abr2y(rij, 1, 1, rij.getX(1), rij.getX(1), du, d2u);
//            z_Hzz_z += r1xH2abr2y(rij, 2, 2, rij.getX(2), rij.getX(2), du, d2u);
//
//            x_Hxy_y += r1xH2abr2y(rij, 0, 1, rij.getX(0), rij.getX(1), du, d2u);
//            x_Hxz_z += r1xH2abr2y(rij, 0, 2, rij.getX(0), rij.getX(2), du, d2u);
//            y_Hyz_z += r1xH2abr2y(rij, 1, 2, rij.getX(1), rij.getX(2), du, d2u);
//
//            x_Hyx_y += r1xH2abr2y(rij, 1, 0, rij.getX(0), rij.getX(1), du, d2u);
//            x_Hzx_z += r1xH2abr2y(rij, 2, 0, rij.getX(0), rij.getX(2), du, d2u);
//            y_Hzy_z += r1xH2abr2y(rij, 2, 1, rij.getX(1), rij.getX(2), du, d2u);
//
//            Vector[] rdot = mapVel(drij);
//
//            dr1_H_dr1 += r1Hr2(rij, rdot[0], rdot[0], du, d2u);
//            dr2_H_dr2 += r1Hr2(rij, rdot[1], rdot[1], du, d2u);
//            dr3_H_dr3 += r1Hr2(rij, rdot[2], rdot[2], du, d2u);
//
//            dr1_H_dr2 += r1Hr2(rij, rdot[0], rdot[1], du, d2u);
//            dr1_H_dr3 += r1Hr2(rij, rdot[0], rdot[2], du, d2u);
//            dr2_H_dr3 += r1Hr2(rij, rdot[1], rdot[2], du, d2u);
//
//            dr4_H_dr4 += r1Hr2(rij, rdot[3], rdot[3], du, d2u);
//            dr5_H_dr5 += r1Hr2(rij, rdot[4], rdot[4], du, d2u);
//            dr6_H_dr6 += r1Hr2(rij, rdot[5], rdot[5], du, d2u);
//
//            Vector rx1 = space.makeVector();
//            rx1.setX(0, rij.getX(0));
//            Vector rx2 = space.makeVector();
//            rx2.setX(1, rij.getX(0));
//            Vector rx3 = space.makeVector();
//            rx3.setX(2, rij.getX(0));
//
//            Vector ry1 = space.makeVector();
//            ry1.setX(0, rij.getX(1));
//            Vector ry2 = space.makeVector();
//            ry2.setX(1, rij.getX(1));
//            Vector ry3 = space.makeVector();
//            ry3.setX(2, rij.getX(1));
//
//            Vector rz1 = space.makeVector();
//            rz1.setX(0, rij.getX(2));
//            Vector rz2 = space.makeVector();
//            rz2.setX(1, rij.getX(2));
//            Vector rz3 = space.makeVector();
//            rz3.setX(2, rij.getX(2));
//
//            Vector ry3z2 = space.makeVector();
//            ry3z2.Ev1Pv2(ry3, rz2);
//            Vector rx3z1 = space.makeVector();
//            rx3z1.Ev1Pv2(rx3, rz1);
//            Vector rx2y1 = space.makeVector();
//            rx2y1.Ev1Pv2(rx2, ry1);
//
//            rx1_H_dr1 += r1Hr2(rij, rx1, rdot[0], du, d2u);
//            ry2_H_dr2 += r1Hr2(rij, ry2, rdot[1], du, d2u);
//            rz3_H_dr3 += r1Hr2(rij, rz3, rdot[2], du, d2u);
//
//            rx1_H_dr2 += r1Hr2(rij, rx1, rdot[1], du, d2u);
//            rx1_H_dr3 += r1Hr2(rij, rx1, rdot[2], du, d2u);
//            ry2_H_dr3 += r1Hr2(rij, ry2, rdot[2], du, d2u);
//
//            ry2_H_dr1 += r1Hr2(rij, ry2, rdot[0], du, d2u);
//            rz3_H_dr1 += r1Hr2(rij, rz3, rdot[0], du, d2u);
//            rz3_H_dr2 += r1Hr2(rij, rz3, rdot[1], du, d2u);
//
//            ry3z2_H_dr4 += r1Hr2(rij, ry3z2, rdot[3], du, d2u);
//            rx3z1_H_dr5 += r1Hr2(rij, rx3z1, rdot[4], du, d2u);
//            rx2y1_H_dr6 += r1Hr2(rij, rx2y1, rdot[5], du, d2u);
//
//            drHdr1 += r1Hr2(rij, drij, rdot[0], du, d2u);
//            drHdr2 += r1Hr2(rij, drij, rdot[1], du, d2u);
//            drHdr3 += r1Hr2(rij, drij, rdot[2], du, d2u);
//
//            rx1_H_dr += r1Hr2(rij, rx1, drij, du, d2u);
//            ry2_H_dr += r1Hr2(rij, ry2, drij, du, d2u);
//            rz3_H_dr += r1Hr2(rij, rz3, drij, du, d2u);
        }
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // Add whatever you need to do with the hessian: elastic!
        IAtom ai = box.getLeafList().get(i);
        IAtom aj = box.getLeafList().get(j);
        dri.Ev1Mv2(ai.getPosition(), coordinteDefinition.getLatticePosition(ai));
        drj.Ev1Mv2(aj.getPosition(), coordinteDefinition.getLatticePosition(aj));
        box.getBoundary().nearestImage(dri);
        box.getBoundary().nearestImage(drj);

        //dri.Hij.drj
        Vector tmpV = space.makeVector();
        tmpV.E(drj);
        Hij.transform(tmpV); //Hij.drj
        drHdr += 2.0*dri.dot(tmpV);

        //self term
        //dri.Hij.dri
        tmpV.E(dri);
        Hij.transform(tmpV);
        drHdr -= dri.dot(tmpV);
        //drj.Hij.drj
        tmpV.E(drj);
        Hij.transform(tmpV);
        drHdr -= drj.dot(tmpV);


        // virial2
//        virial2 += r1Hr2(rij, rij, rij, du, d2u);

    }




    protected Vector[] mapVel(Vector dr) {
        Vector[] rdot = new Vector[6];
        for (int i = 0; i < 6; i++) {
            rdot[i] = space.makeVector();
        }
        rdot[0].setX(0, mx1 * dr.getX(0)); //rdot1
        rdot[0].setX(1, my1 * dr.getX(1));
        rdot[0].setX(2, my1 * dr.getX(2));

        rdot[1].setX(0, my1 * dr.getX(0)); //rdot2
        rdot[1].setX(1, mx1 * dr.getX(1));
        rdot[1].setX(2, my1 * dr.getX(2));

        rdot[2].setX(0, my1 * dr.getX(0));//rdot3
        rdot[2].setX(1, my1 * dr.getX(1));
        rdot[2].setX(2, mx1 * dr.getX(2));

        rdot[3].setX(1, dr.getX(2));//rdot4
        rdot[3].setX(2, dr.getX(1));
        rdot[3].TE(my4);

        rdot[4].setX(0, dr.getX(2));//rdot5
        rdot[4].setX(2, dr.getX(0));
        rdot[4].TE(my4);

        rdot[5].setX(0, dr.getX(1));//rdot6
        rdot[5].setX(1, dr.getX(0));
        rdot[5].TE(my4);

        return rdot;
    }

    protected Vector[] mapAcc(Vector dr) {
        Vector[] rddot = new Vector[12];
        for (int i = 0; i < 12; i++) {
            rddot[i] = space.makeVector();
        }

        rddot[0].setX(0, mx11 * dr.getX(0));//rddot1
        rddot[0].setX(1, my11 * dr.getX(1));
        rddot[0].setX(2, my11 * dr.getX(2));

        rddot[1].setX(0, my11 * dr.getX(0));//rddot2
        rddot[1].setX(1, mx11 * dr.getX(1));
        rddot[1].setX(2, my11 * dr.getX(2));

        rddot[2].setX(0, my11 * dr.getX(0));//rddot3
        rddot[2].setX(1, my11 * dr.getX(1));
        rddot[2].setX(2, mx11 * dr.getX(2));

        rddot[3].setX(0, mx44 * dr.getX(0));//rddot4
        rddot[3].setX(1, my44 * dr.getX(1));
        rddot[3].setX(2, my44 * dr.getX(2));

        rddot[4].setX(0, my44 * dr.getX(0));//rddot5
        rddot[4].setX(1, mx44 * dr.getX(1));
        rddot[4].setX(2, my44 * dr.getX(2));

        rddot[5].setX(0, my44 * dr.getX(0));//rddot6
        rddot[5].setX(1, my44 * dr.getX(1));
        rddot[5].setX(2, mx44 * dr.getX(2));

        rddot[6].setX(0, mz12 * dr.getX(0));//dr23
        rddot[6].setX(1, mx12 * dr.getX(1));
        rddot[6].setX(2, mx12 * dr.getX(2));

        rddot[7].setX(0, mx12 * dr.getX(0));//rddot13
        rddot[7].setX(1, mz12 * dr.getX(1));
        rddot[7].setX(2, mx12 * dr.getX(2));

        rddot[8].setX(0, mx12 * dr.getX(0));//rddot12
        rddot[8].setX(1, mx12 * dr.getX(1));
        rddot[8].setX(2, mz12 * dr.getX(2));

        rddot[9].setX(0, mx1T * dr.getX(0));
        rddot[9].setX(1, my1T * dr.getX(1));
        rddot[9].setX(2, my1T * dr.getX(2));

        rddot[10].setX(0, my1T * dr.getX(0));
        rddot[10].setX(1, mx1T * dr.getX(1));
        rddot[10].setX(2, my1T * dr.getX(2));

        rddot[11].setX(0, my1T * dr.getX(0));
        rddot[11].setX(1, my1T * dr.getX(1));
        rddot[11].setX(2, mx1T * dr.getX(2));

        return rddot;
    }

    protected double r1Hr2(Vector rij, Vector r1ij, Vector r2ij, double du, double d2u) { //r1.H.r2
        double rij2 = rij.squared();
        double virial2 = du / rij2 * r1ij.dot(r2ij) + (d2u - du) / rij2 / rij2 * r1ij.dot(rij) * (r2ij.dot(rij));
        return virial2;
    }

    protected double r1xH2abr2y(Vector rij, int a, int b, double x1, double y2, double du, double d2u) {
        double rij2 = rij.squared();
        double raH2rb = (d2u - du) / rij2 / rij2 * x1 * y2 * rij.getX(a) * rij.getX(b);
        return raH2rb;
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
        return doD2;
    }

}