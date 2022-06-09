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
import etomica.data.IDataSourcePotential;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;
import org.apache.xml.utils.IntVector;

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
    protected double rPhir, drPhidr, rPhidr;
    protected final Space space;
    protected double uShift, pShift;

    protected double virialx, virialy, virialz, virialxy, virialxz, virialyz;
    protected double x_Phixx_x, y_Phiyy_y, z_Phizz_z;
    protected double x_Phixy_y, x_Phixz_z, y_Phiyz_z;
    protected double x_Phiyx_y, x_Phizx_z, y_Phizy_z;

    protected double dr1_Phi_dr1, dr2_Phi_dr2, dr3_Phi_dr3, dr1_Phi_dr2, dr1_Phi_dr3, dr2_Phi_dr3;
    protected double dr4_Phi_dr4, dr5_Phi_dr5, dr6_Phi_dr6, rx1_Phi_dr1, ry2_Phi_dr2, rz3_Phi_dr3;
    protected double rx1_Phi_dr2, rx1_Phi_dr3, ry2_Phi_dr3, ry2_Phi_dr1, rz3_Phi_dr1, rz3_Phi_dr2;
    protected double ry3z2_Phi_dr4, rx3z1_Phi_dr5, rx2y1_Phi_dr6;
    protected double drPhidr1, drPhidr2, drPhidr3, rx1_Phi_dr, ry2_Phi_dr, rz3_Phi_dr;

    protected double gV, gVV, mV, mVV;
    protected double gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
    protected double mx1, my1, my4, mx11, my11, mx12, mz12, mx44, my44;
    protected double mT, mVT, mx1T, my1T;
    protected final boolean doD2;
    protected Vector dri, drj, drij, Rij;

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
        int n = doD2 ? 41 : 16;

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

        mV = gV - 1.0 / 3.0;
        mVV = gVV + gV * gV + 2.0 / 9.0;
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
        Rij = space.makeVector();

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
        rPhir = drPhidr = rPhidr = 0;
        x_Phixx_x = y_Phiyy_y = z_Phizz_z = x_Phixy_y = x_Phixz_z = y_Phiyz_z = x_Phiyx_y = x_Phizx_z = y_Phizy_z = 0;
        dr1_Phi_dr1 = dr2_Phi_dr2 = dr3_Phi_dr3 = dr1_Phi_dr2 = dr1_Phi_dr3 = dr2_Phi_dr3 = 0;
        dr4_Phi_dr4 = dr5_Phi_dr5 = dr6_Phi_dr6 = rx1_Phi_dr1 = ry2_Phi_dr2 = rz3_Phi_dr3 = 0;
        drPhidr1 = drPhidr2 = drPhidr3 = rx1_Phi_dr = ry2_Phi_dr = rz3_Phi_dr = 0;
        rx1_Phi_dr2 = rx1_Phi_dr3 = ry2_Phi_dr3 = ry2_Phi_dr1 = rz3_Phi_dr1 = rz3_Phi_dr2 = 0;
        ry3z2_Phi_dr4 = rx3z1_Phi_dr5 = rx2y1_Phi_dr6 = 0;


        double uSum;
        if (callComputeAll) {
            PotentialCallback callback = doD2 ? this : null;
            uSum = potentialMaster.computeAll(true, callback);
        } else {
            uSum = potentialMaster.getLastEnergy();
        }

        double[] x = data.getData();
        double V = box.getBoundary().volume();
        int N = box.getMoleculeList().size();
        double rho = N / V;
        double virial = potentialMaster.getLastVirial();
        Vector[] forces = potentialMaster.getForces();
        IAtomList atoms = box.getLeafList();
        for (IAtom a : atoms) {
            Vector F = forces[a.getLeafIndex()];
            dri.Ev1Mv2(a.getPosition(), coordinteDefinition.getLatticePosition(a));
            box.getBoundary().nearestImage(dri);

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
        x[3] = 3 * (N - 1) / V * temperature * gV - virial / (3 * V) + mV * Fdr / V - pShift;
        double p1Shift = -pShift;
        //P1
        x[4] = -rho * temperature + virialx / V - p1Shift;
        ;
        x[5] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialx + Frdot1) / V - p1Shift;
        //P2
        x[6] = -rho * temperature + virialy / V - p1Shift;
        ;
        x[7] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialy + Frdot2) / V - p1Shift;
        //P3
        x[8] = -rho * temperature + virialz / V - p1Shift;
        ;
        x[9] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialz + Frdot3) / V - p1Shift;
        ;
        //P4
        x[10] = virialyz / V;
        x[11] = (virialyz - Frdot4) / V;
        //P5
        x[12] = virialxz / V;
        x[13] = (virialxz - Frdot5) / V;
        //P6
        x[14] = virialxy / V;
        x[15] = (virialxy - Frdot6) / V;

        if (doD2) {
            //Cv - HMA
            x[16] = -1.0 / 4.0 / temperature * (Fdr + drPhidr);

            // B
            x[17] = (rPhir - 2.0 * virial) / (9 * V) + rho * temperature;
            x[18] = ((rPhir - 2.0 * virial) / 9.0 + mV * mV * drPhidr - mVV * Fdr + 2.0 / 3.0 * mV * rPhidr - 3.0 * (N - 1) * temperature * gVV) / V + temperature / V;

            //C11
            x[19] = x_Phixx_x / V + 2 * rho * temperature; //C11
            x[20] = (x_Phixx_x - Frddot11 + dr1_Phi_dr1 + 2 * rx1_Phi_dr1 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
            //C22
            x[21] = y_Phiyy_y / V + 2 * rho * temperature; //C11
            x[22] = (y_Phiyy_y - Frddot22 + dr2_Phi_dr2 + 2 * ry2_Phi_dr2 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
            //C33
            x[23] = z_Phizz_z / V + 2 * rho * temperature; //C11
            x[24] = (z_Phizz_z - Frddot33 + dr3_Phi_dr3 + 2 * rz3_Phi_dr3 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;

            //C12
            x[25] = x_Phixy_y / V;//C12
            x[26] = (x_Phixy_y - Frddot12 + dr1_Phi_dr2 + rx1_Phi_dr2 + ry2_Phi_dr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
            //C13
            x[27] = x_Phixz_z / V;//C12
            x[28] = (x_Phixz_z - Frddot13 + dr1_Phi_dr3 + rx1_Phi_dr3 + rz3_Phi_dr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
            //C23
            x[29] = y_Phiyz_z / V;//C12
            x[30] = (y_Phiyz_z - Frddot23 + dr2_Phi_dr2 + ry2_Phi_dr3 + rz3_Phi_dr2 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;

            //C44
            x[31] = y_Phizy_z / V + rho * temperature;
            x[32] = (y_Phizy_z - Frddot44 + dr4_Phi_dr4 + ry3z2_Phi_dr4 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
            //C55
            x[33] = x_Phizx_z / V + rho * temperature;//C55
            x[34] = (x_Phizx_z - Frddot55 + dr5_Phi_dr5 + rx3z1_Phi_dr5 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
            //C6
            x[35] = x_Phiyx_y / V + rho * temperature;//C66
            x[36] = (x_Phiyx_y - Frddot66 + dr6_Phi_dr6 + rx2y1_Phi_dr6 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;

            //b_V
            x[37] = 1.0 / 2.0 / V / temperature * (gV * Fdr - mV * drPhidr - 1.0 / 3.0 * rPhidr) + 3.0 * (N - 1) * gV / V + 1.0 / V; //bV_hma
            //b_mn
            x[38] = 1.0 / V * (-Frddot1T + mT * drPhidr1 + mT * rx1_Phi_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
            x[39] = 1.0 / V * (-Frddot2T + mT * drPhidr2 + mT * ry2_Phi_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
            x[40] = 1.0 / V * (-Frddot3T + mT * drPhidr3 + mT * rz3_Phi_dr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
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
        Vector Ri = coordinteDefinition.getLatticePosition(ia);
        Vector Rj = coordinteDefinition.getLatticePosition(ja);
        Rij.Ev1Mv2(Ri, Rj);

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

        rPhir += r1Phir2(rij, rij, rij, du, d2u);
        drPhidr += r1Phir2(rij, drij, drij, du, d2u);
        rPhidr += r1Phir2(rij, rij, drij, du, d2u);

        x_Phixx_x += r1xPhi2abr2y(rij, 0, 0, rij.getX(0), rij.getX(0), du, d2u);
        y_Phiyy_y += r1xPhi2abr2y(rij, 1, 1, rij.getX(1), rij.getX(1), du, d2u);
        z_Phizz_z += r1xPhi2abr2y(rij, 2, 2, rij.getX(2), rij.getX(2), du, d2u);

        x_Phixy_y += r1xPhi2abr2y(rij, 0, 1, rij.getX(0), rij.getX(1), du, d2u);
        x_Phixz_z += r1xPhi2abr2y(rij, 0, 2, rij.getX(0), rij.getX(2), du, d2u);
        y_Phiyz_z += r1xPhi2abr2y(rij, 1, 2, rij.getX(1), rij.getX(2), du, d2u);

        x_Phiyx_y += r1xPhi2abr2y(rij, 1, 0, rij.getX(0), rij.getX(1), du, d2u);
        x_Phizx_z += r1xPhi2abr2y(rij, 2, 0, rij.getX(0), rij.getX(2), du, d2u);
        y_Phizy_z += r1xPhi2abr2y(rij, 2, 1, rij.getX(1), rij.getX(2), du, d2u);

        Vector[] rdot = mapVel(drij);

        dr1_Phi_dr1 += r1Phir2(rij, rdot[0], rdot[0], du, d2u);
        dr2_Phi_dr2 += r1Phir2(rij, rdot[1], rdot[1], du, d2u);
        dr3_Phi_dr3 += r1Phir2(rij, rdot[2], rdot[2], du, d2u);

        dr1_Phi_dr2 += r1Phir2(rij, rdot[0], rdot[1], du, d2u);
        dr1_Phi_dr3 += r1Phir2(rij, rdot[0], rdot[2], du, d2u);
        dr2_Phi_dr3 += r1Phir2(rij, rdot[1], rdot[2], du, d2u);

        dr4_Phi_dr4 += r1Phir2(rij, rdot[3], rdot[3], du, d2u);
        dr5_Phi_dr5 += r1Phir2(rij, rdot[4], rdot[4], du, d2u);
        dr6_Phi_dr6 += r1Phir2(rij, rdot[5], rdot[5], du, d2u);

        Vector rx1 = space.makeVector();
        rx1.setX(0, rij.getX(0));
        Vector rx2 = space.makeVector();
        rx2.setX(1, rij.getX(0));
        Vector rx3 = space.makeVector();
        rx3.setX(2, rij.getX(0));

        Vector ry1 = space.makeVector();
        ry1.setX(0, rij.getX(1));
        Vector ry2 = space.makeVector();
        ry2.setX(1, rij.getX(1));
        Vector ry3 = space.makeVector();
        ry3.setX(2, rij.getX(1));

        Vector rz1 = space.makeVector();
        rz1.setX(0, rij.getX(2));
        Vector rz2 = space.makeVector();
        rz2.setX(1, rij.getX(2));
        Vector rz3 = space.makeVector();
        rz3.setX(2, rij.getX(2));

        Vector ry3z2 = space.makeVector();
        ry3z2.Ev1Pv2(ry3, rz2);
        Vector rx3z1 = space.makeVector();
        rx3z1.Ev1Pv2(rx3, rz1);
        Vector rx2y1 = space.makeVector();
        rx2y1.Ev1Pv2(rx2, ry1);

        rx1_Phi_dr1 += r1Phir2(rij, rx1, rdot[0], du, d2u);
        ry2_Phi_dr2 += r1Phir2(rij, ry2, rdot[1], du, d2u);
        rz3_Phi_dr3 += r1Phir2(rij, rz3, rdot[2], du, d2u);

        rx1_Phi_dr2 += r1Phir2(rij, rx1, rdot[1], du, d2u);
        rx1_Phi_dr3 += r1Phir2(rij, rx1, rdot[2], du, d2u);
        ry2_Phi_dr3 += r1Phir2(rij, ry2, rdot[2], du, d2u);

        ry2_Phi_dr1 += r1Phir2(rij, ry2, rdot[0], du, d2u);
        rz3_Phi_dr1 += r1Phir2(rij, rz3, rdot[0], du, d2u);
        rz3_Phi_dr2 += r1Phir2(rij, rz3, rdot[1], du, d2u);

        ry3z2_Phi_dr4 += r1Phir2(rij, ry3z2, rdot[3], du, d2u);
        rx3z1_Phi_dr5 += r1Phir2(rij, rx3z1, rdot[4], du, d2u);
        rx2y1_Phi_dr6 += r1Phir2(rij, rx2y1, rdot[5], du, d2u);

        drPhidr1 += r1Phir2(rij, drij, rdot[0], du, d2u);
        drPhidr2 += r1Phir2(rij, drij, rdot[1], du, d2u);
        drPhidr3 += r1Phir2(rij, drij, rdot[2], du, d2u);

        rx1_Phi_dr += r1Phir2(rij, rx1, drij, du, d2u);
        ry2_Phi_dr += r1Phir2(rij, ry2, drij, du, d2u);
        rz3_Phi_dr += r1Phir2(rij, rz3, drij, du, d2u);
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

    protected double r1Phir2(Vector rij, Vector r1ij, Vector r2ij, double du, double d2u) { //r1.Phi.r2
        double rij2 = rij.squared();
        double rPhir = du / rij2 * r1ij.dot(r2ij) + (d2u - du) / rij2 / rij2 * r1ij.dot(rij) * (r2ij.dot(rij));
        return rPhir;
    }

    protected double r1xPhi2abr2y(Vector rij, int a, int b, double x1, double y2, double du, double d2u) {
        double rij2 = rij.squared();
        double raPhi2rb = (d2u - du) / rij2 / rij2 * x1 * y2 * rij.getX(a) * rij.getX(b);
        return raPhi2rb;
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


    public void pairComputeHessian(int i, int j, Tensor Hij) { // Add whatever you need to do with the hessian: elastic!
        IAtom ai = box.getLeafList().getAtoms().get(i);
        IAtom aj = box.getLeafList().getAtoms().get(j);
        dri.Ev1Mv2(ai.getPosition(), coordinteDefinition.getLatticePosition(ai));
        drj.Ev1Mv2(aj.getPosition(), coordinteDefinition.getLatticePosition(aj));
        box.getBoundary().nearestImage(dri);
        box.getBoundary().nearestImage(drj);

        //i != j
        Vector tmpV = space.makeVector();
        tmpV.E(drj);
        Hij.transform(tmpV); //Hij.drj
        drPhidr += 2*dri.dot(tmpV); //dri.Hij.drj
        //self term
        tmpV.E(dri);
        Hij.transform(tmpV);
        drPhidr -= dri.dot(tmpV); //dri.Hij.drj

        tmpV.E(drj);
        Hij.transpose(); //Hij^T = Hji
        Hij.transform(tmpV);
        drPhidr -= dri.dot(tmpV); //drj.Hji.drj
    }
}