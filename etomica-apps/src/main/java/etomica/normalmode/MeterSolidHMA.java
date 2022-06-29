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
import etomica.units.Pascal;
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
    protected double drHdr, rHdr;
    protected final Space space;
    protected double uShift, pShift;

    protected double virial2;
    protected double virialx, virialy, virialz, virialxy, virialxz, virialyz;

    protected double[] virialXYZ, virial2XYZ;
    protected double x_Hxx_x, y_Hyy_y, z_Hzz_z;
    protected double x_Hxy_y, x_Hxz_z, y_Hyz_z;
    protected double x_Hyx_y, x_Hzx_z, y_Hzy_z;

    protected double[] rdotHrdot =  new double[9];
    protected double rx1Hdr1, ry2Hdr2, rz3Hdr3;
    protected double rx1Hdr2, rx1Hdr3, ry2Hdr3, ry2Hdr1, rz3Hdr1, rz3Hdr2;
    protected double ry3z2Hdr4, rx3z1Hdr5, rx2y1Hdr6;
    protected double drHdr1, drHdr2, drHdr3, rx1Hdr, ry2Hdr, rz3Hdr;

    protected double gV, gVV, mV, mVV;


    protected double fV , fVV, dP, dB;

    protected double gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
    protected double mx1, my1, my4, mx11, my11, mx12, mz12, mx44, my44;
    protected double mT, mVT, mx1T, my1T;
    protected final boolean doD2;
    protected Vector dri, drj, drij;

    protected boolean callComputeAll = true;


    protected double sumH = 0;
    protected int nn = 65;

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
        int nData = doD2 ? 37 : 16;

        dataInfo = new DataInfoDoubleArray("Stuff", Null.DIMENSION, new int[]{nData});
        dataInfo.addTag(tag);
        data = new DataDoubleArray(nData);

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

        dP = elasticParams[11];
        dB = elasticParams[12];

        System.out.println(" gV: " + gV + " gVV: " + gVV);
        System.out.println(" gx1: " + gx1 + " gy1: " + gy1 + " gy4 " + gy4);
        System.out.println(" gx11: " + gx11 + " gy11: " + gy11 + " gx44 " + gx44);
        System.out.println(" gy44: " + gy44 + " gx12: " + gx12 + " gz12 " + gz12);
        System.out.println(" dP: " + dP + " dB: "+ dB);



        int numAtoms = box.getMoleculeList().size();
        double V = box.getBoundary().volume();
        double rho = numAtoms / V;

//        gV = 0.894175378798851;// DB
//        gVV=1.658811488865302; //DB

//        gV = 1.282981055972710;//FB
//        gVV = 0.464931126443066;//FB

        mV = gV - 1.0 / 3.0;
        mVV = gVV + gV * gV + 2.0 / 9.0;

        //fitting: N=500, rho=0.1

        fV = (dP/temperature-rho)/(3*(numAtoms-1));
        fVV = fV*fV - (dB/temperature-rho)/(3*V*(numAtoms-1)) + 2.0*fV/(3.0*V);

//        mVV = mV*mV - V*(dB/temperature-rho)/(3*(numAtoms-1));
//        gVV = mVV - gV*gV -2.0/3.0;

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


        virial2 = 0 ;


        virialx = virialy = virialz = virialxy = virialxz = virialyz = 0;
        drHdr = rHdr = 0;
        x_Hxx_x = y_Hyy_y = z_Hzz_z = x_Hxy_y = x_Hxz_z = y_Hyz_z = x_Hyx_y = x_Hzx_z = y_Hzy_z = 0;

        for (int i=0; i<9; i++){
            rdotHrdot[i] = 0;
        }

        rx1Hdr1 = ry2Hdr2 = rz3Hdr3 = 0;
        drHdr1 = drHdr2 = drHdr3 = rx1Hdr = ry2Hdr = rz3Hdr = 0;
        rx1Hdr2 = rx1Hdr3 = ry2Hdr3 = ry2Hdr1 = rz3Hdr1 = rz3Hdr2 = 0;
        ry3z2Hdr4 = rx3z1Hdr5 = rx2y1Hdr6 = 0;

        double uSum;
        if (callComputeAll) {
            PotentialCallback callback = doD2 ? this : null;
            sumH = 0;
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
        double[] virialXYZ = potentialMaster.getLastVirialXYZ();
        double[] virial2XYZ = potentialMaster.getLastVirial2XYZ();

        Vector[] forces = potentialMaster.getForces();
        Vector[] dFdeV = potentialMaster.getdFdeV();
        Vector[][] dFde = potentialMaster.getdFde();


        boolean debug = false;
        if(debug && wantsHessian() && dFdeV != null){

            if (doD2){
//                Vector dr = space.makeVector();
//                double dx = 0.0001;
//                dr.setX(0, dx);
//                double forces_00 = potentialMaster.getForces()[nn].getX(0);
//
////                System.out.println("F00: " + forces_00);
//
//                box.getLeafList().get(0).getPosition().PEa1Tv1(-1.0, dr);
//                potentialMaster.computeAll(true, this);
//                double forces_m = potentialMaster.getForces()[nn].getX(0);
////                System.out.println("Fm: " + forces_m);
//
//                box.getLeafList().get(0).getPosition().PEa1Tv1(1.0, dr);
//                potentialMaster.computeAll(true, this);
//                double forces_0 = potentialMaster.getForces()[nn].getX(0);
////                System.out.println("F0: " + forces_0);
//
//                box.getLeafList().get(0).getPosition().PEa1Tv1(1.0, dr);
//                potentialMaster.computeAll(true, this);
//                double forces_p = potentialMaster.getForces()[nn].getX(0);
////                System.out.println("Fp: " + forces_p);
//
//                System.out.println("H01: " + (-(forces_p-forces_m)/2/dx));
//                System.out.println();
//                System.exit(0);

//
                double dFdeVnnx = dFdeV[nn].getX(0);
//                System.out.println(" dFdeV: " + dFdeV[nn].getX(0));
//                System.out.println(box.getLeafList().get(nn).getPosition());

                double deV = 0.00001;
                double Vm = V*(1.0-deV);
                double Vmm = V*(1.0-2.0*deV);
                double Vp = V*(1.0+deV);
                double Vpp = V*(1.0+2.0*deV);
                double rhom = N/Vm;
                double rhomm = N/Vmm;
                double rhop = N/Vp;
                double rhopp = N/Vpp;
                BoxInflate inflater = new BoxInflate(box, space);


                inflater.setTargetDensity(rhomm);
                inflater.actionPerformed();
                potentialMaster.init();
//                System.out.println(" minus2: V = " + Vmm);
                potentialMaster.computeAll(true, this);
                double U_mm = potentialMaster.getLastEnergy();
                double forces_mm = potentialMaster.getForces()[nn].getX(0);

                inflater.setTargetDensity(rhom);
                inflater.actionPerformed();
                potentialMaster.init();
//                System.out.println(" minus1: V = " + Vm);
                potentialMaster.computeAll(true, this);
                double forces_m = potentialMaster.getForces()[nn].getX(0);
                double U_m = potentialMaster.getLastEnergy();
//                System.out.println(box.getLeafList().get(nn).getPosition());


                inflater.setTargetDensity(rho);
                inflater.actionPerformed();
                potentialMaster.init();
                potentialMaster.computeAll(true, this);
                double forces_0 = potentialMaster.getForces()[nn].getX(0);
                double U_0 = potentialMaster.getLastEnergy();


//                System.out.println(" plus1: V = " + Vp);
                inflater.setTargetDensity(rhop);
                inflater.actionPerformed();
                potentialMaster.init();
                potentialMaster.computeAll(true, this);
                double forces_p = potentialMaster.getForces()[nn].getX(0);
                double U_p = potentialMaster.getLastEnergy();
//                System.out.println(box.getLeafList().get(nn).getPosition());




                inflater.setTargetDensity(rhopp);
                inflater.actionPerformed();
                potentialMaster.init();
//                System.out.println(" minus2: V = " + Vmm);
                potentialMaster.computeAll(true, this);
                double U_pp = potentialMaster.getLastEnergy();
                double forces_pp = potentialMaster.getForces()[nn].getX(0);

                System.out.println("\n FD-Fx: " + " "+ (((forces_p-forces_m)/2/deV)-dFdeVnnx)/dFdeVnnx*100);


//                System.out.println("\nn Virial2-FD: " + (U_m+U_p-2*U_0)/deV/deV + "  " + virial2);
                double fd2 = (U_m+U_p-2*U_0)/deV/deV;
                double fd5 = (-U_pp+16*U_p-30*U_0+16*U_m-U_mm)/(12*deV*deV);
//                System.out.println("\n Virial2-FD2: " + fd2 + "  " + virial2);
//                System.out.println(" Virial2-FD5: " + fd5 + "  " + virial2);
//                System.exit(0);

            }

        } //if

        double dFdeVdr = 0, dF1rdot1 = 0, dF2rdot2 = 0, dF3rdot3 = 0, dF4rdot4 = 0, dF5rdot5 = 0, dF6rdot6 = 0;
        double dF1rdot2 = 0, dF1rdot3 = 0, dF2rdot3 = 0, dF2rdot1 = 0, dF3rdot1 = 0, dF3rdot2 = 0;
        IAtomList atoms = box.getLeafList();
        for (IAtom a : atoms) {
            Vector F = forces[a.getLeafIndex()];
            dri.Ev1Mv2(a.getPosition(), coordinteDefinition.getLatticePosition(a));
            box.getBoundary().nearestImage(dri);

            Vector[] rdot = mapVel(dri);
            Vector[] rddot = mapAcc(dri);


            if (doD2) {
                dFdeVdr += dFdeV[a.getLeafIndex()].dot(dri);
                //00 , 01 , 02 , 11 , 12 , 22
                //0    1    2    3    4    5
                dF1rdot1 += dFde[0][a.getLeafIndex()].dot(rdot[0]);
                dF2rdot2 += dFde[3][a.getLeafIndex()].dot(rdot[1]);
                dF3rdot3 += dFde[5][a.getLeafIndex()].dot(rdot[2]);

                dF1rdot2 += dFde[0][a.getLeafIndex()].dot(rdot[1]);
                dF2rdot1 += dFde[3][a.getLeafIndex()].dot(rdot[0]);

                dF1rdot3 += dFde[0][a.getLeafIndex()].dot(rdot[2]);
                dF3rdot1 += dFde[5][a.getLeafIndex()].dot(rdot[0]);

                dF2rdot3 += dFde[3][a.getLeafIndex()].dot(rdot[2]);
                dF3rdot2 += dFde[5][a.getLeafIndex()].dot(rdot[1]);

                dF4rdot4 += dFde[4][a.getLeafIndex()].dot(rdot[3]); //dF/de4 * rdot3
                dF5rdot5 += dFde[2][a.getLeafIndex()].dot(rdot[4]);
                dF6rdot6 += dFde[1][a.getLeafIndex()].dot(rdot[5]);
            }

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
        //PRE
//        x[3] = dP - virial / (3 * V) + fV * Fdr - pShift; //PRE
        //PRB
         x[3] = 3 * (N - 1) / V * temperature * gV - virial / (3 * V) + mV * Fdr / V - pShift; //PRB


        // V*(dP/temperature-rho)/(dim*(numAtoms-1)) + 1/3.0;

        double p1Shift = -pShift;
        //P1
//        x[4] = -rho * temperature + virialx / V - p1Shift;
//        x[5] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialx + Frdot1) / V - p1Shift;
        x[4] = -rho * temperature + virialXYZ[0] / V - p1Shift;
        x[5] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialXYZ[0] + Frdot1) / V - p1Shift;


        //P2
        x[6] = -rho * temperature + virialXYZ[1] / V - p1Shift;
        x[7] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialXYZ[1] + Frdot2) / V - p1Shift;
        //P3
        x[8] = -rho * temperature + virialXYZ[2] / V - p1Shift;
        x[9] = -((N - 1) * temperature * (gx1 + 2 * gy1) - virialXYZ[2] + Frdot3) / V - p1Shift;
        //P4
//        x[10] = virialyz / V;
//        x[11] = (virialyz - Frdot4) / V;
        x[10] = virialXYZ[3] / V;
        x[11] = (virialXYZ[3] - Frdot4) / V;
        //P5
        x[12] = virialXYZ[4] / V;
        x[13] = (virialXYZ[4] - Frdot5) / V;
        //P6
        x[14] = virialXYZ[5] / V;
        x[15] = (virialXYZ[5] - Frdot6) / V;


        //Second derivatives
        if (doD2) {
            //Cv - HMA
            x[16] = -1.0/4.0/temperature*(Fdr + drHdr);

            // B
            x[17] = virial2/V + rho*temperature;
            // LJ - old
            // x[18] = (virial2 + mV * mV * drHdr - mVV * Fdr + 2.0 / 3.0 * mV * rHdr - 3.0 * (N - 1) * temperature * gVV) / V + temperature / V;
            //PRE
//            x[18] = dB + virial2/V - V*fVV*Fdr - 2.0*fV*dFdeVdr  + V*fV*fV*drHdr ;
            //PRB
            x[18] = (virial2 + mV*mV*drHdr - mVV*Fdr - 2.0*mV*dFdeVdr - 3.0*(N-1)*temperature*gVV)/V + temperature/V;

            // dr1Hdr1 = 0 , dr2Hdr2 = 3 , dr3Hdr3 = 5
            // dr1Hdr2 = 1 , dr1Hdr3 = 2 , dr2Hdr3 = 4
            // dr4Hdr4 = 6 , dr5Hdr5 = 7 , dr6Hdr6 = 8


            //C11
//            x[19] = x_Hxx_x / V + 2 * rho * temperature; //C11
//            x[20] = (x_Hxx_x + dr1Hdr1 - Frddot11 + 2 * rx1Hdr1 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
            x[19] = virial2XYZ[0]/V + 2*rho*temperature; //C11
            x[20] = (virial2XYZ[0] + rdotHrdot[0] - Frddot11  - 2*dF1rdot1 - (N-1.0)*temperature*(gx11+2*gy11))/V + 2*temperature/V;
            //C22
//            x[21] = y_Hyy_y / V + 2 * rho * temperature; //C22
//            x[22] = (y_Hyy_y + dr2Hdr2 - Frddot22 + 2 * ry2Hdr2 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
            x[21] = virial2XYZ[1]/V + 2*rho*temperature; //C22
            x[22] = (virial2XYZ[1] + rdotHrdot[3] - Frddot22  - 2*dF2rdot2 - (N-1.0)*temperature*(gx11+2*gy11))/V + 2*temperature/V;
//            //C33
//            x[23] = z_Hzz_z / V + 2 * rho * temperature; //C33
//            x[24] = (z_Hzz_z + dr3Hdr3 - Frddot33 + 2 * rz3Hdr3 - (N - 1.0) * temperature * (gx11 + 2 * gy11)) / V + 2 * temperature / V;
            x[23] = virial2XYZ[2]/V + 2*rho*temperature; //C33
            x[24] = (virial2XYZ[2] + rdotHrdot[5] - Frddot33  - 2*dF3rdot3 - (N-1.0)*temperature*(gx11+2*gy11))/V + 2*temperature/V;

//
//            //C12
//            x[25] = x_Hxy_y / V;//C12
//            x[26] = (x_Hxy_y + dr1Hdr2 - Frddot12 + rx1Hdr2 + ry2Hdr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
//            //C13
//            x[27] = x_Hxz_z / V;//C13
//            x[28] = (x_Hxz_z - Frddot13 + dr1Hdr3 + rx1Hdr3 + rz3Hdr1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
//            //C23
//            x[29] = y_Hyz_z / V;//C23
//            x[30] = (y_Hyz_z - Frddot23 + dr2Hdr2 + ry2Hdr3 + rz3Hdr2 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;

            // dr1Hdr2 = 1 , dr1Hdr3 = 2 , dr2Hdr3 = 4
            //C12
            x[25] = virial2XYZ[3] / V;
            x[26] = (virial2XYZ[3] + rdotHrdot[1] - Frddot12 - dF1rdot2 - dF2rdot1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
            //C13
            x[27] = virial2XYZ[4] / V;
            x[28] = (virial2XYZ[4] + rdotHrdot[2] - Frddot13 - dF1rdot3 - dF3rdot1 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;
            //C23
            x[29] = virial2XYZ[5] / V;
            x[30] = (virial2XYZ[5] + rdotHrdot[4] - Frddot23 - dF2rdot3 - dF3rdot2 - (N - 1.0) * temperature * (gz12 + 2 * gx12)) / V;

//            //C44
//            x[31] = y_Hzy_z / V + rho * temperature;
//            x[32] = (y_Hzy_z + dr4Hdr4 - Frddot44 + ry3z2Hdr4 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//            //C55
//            x[33] = x_Hzx_z / V + rho * temperature;//C55
//            x[34] = (x_Hzx_z + dr5Hdr5 - Frddot55 + rx3z1Hdr5 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//            //C6
//            x[35] = x_Hyx_y / V + rho * temperature;//C66
//            x[36] = (x_Hyx_y + dr6Hdr6 - Frddot66 + rx2y1Hdr6 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
//          C44 (yz,yz)
            x[31] = virial2XYZ[6] / V + rho * temperature;
            x[32] = (virial2XYZ[6] + rdotHrdot[6] - Frddot44 - 2.0*dF4rdot4 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
            //C55
            x[33] = virial2XYZ[7] / V + rho * temperature;
            x[34] = (virial2XYZ[7] + rdotHrdot[7] - Frddot55 - 2.0*dF5rdot5 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;
            //C66
            x[35] = virial2XYZ[8] / V + rho * temperature;
            x[36] = (virial2XYZ[8] + rdotHrdot[8] - Frddot66 - 2.0*dF6rdot6 - (N - 1) * temperature * (gx44 + 2 * gy44) + temperature) / V;

//            //b_V
//            x[37] = 1.0 / 2.0 / V / temperature * (gV * Fdr - mV * drHdr - 1.0 / 3.0 * rHdr) + 3.0 * (N - 1) * gV / V + 1.0 / V; //bV_hma




            //  rHdr -> -3 dFdeV/V
//            x[19] = 1.0 / 2.0 / V / temperature * (gV * Fdr - mV * drHdr + dFdeVdr) + 3.0 * (N - 1) * gV / V + 1.0 / V; //bV_hma



//            //b_mn
//            x[38] = 1.0 / V * (-Frddot1T + mT * drHdr1 + mT * rx1Hdr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
//            x[39] = 1.0 / V * (-Frddot2T + mT * drHdr2 + mT * ry2Hdr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
//            x[40] = 1.0 / V * (-Frddot3T + mT * drHdr3 + mT * rz3Hdr) - (N - 1) / V * (gx1 + 2 * gy1) - 1 / V;
        }

        return data;
    }

    @Override
    //pairCompute called by potentialcompute
    public void pairCompute(int iAtom, int jAtom, Vector rij, double[] u012) { //rij = rj - ri
        double du = u012[1];
        double d2u = u012[2];
        double rij2 = rij.squared();

        boolean debug = false;
        if (debug && doD2 && iAtom==0 && jAtom==nn) {
            Tensor unity = space.makeTensor();
            unity.setComponent(0, 0, 1.0);
            unity.setComponent(1, 1, 1.0);
            unity.setComponent(2, 2, 1.0);

            Tensor t = space.makeTensor();
            t.Ev1v2(rij,rij);
            double fac = (du-d2u)/rij2/rij2 ;
            t.TE(fac);
            fac = -du/rij2;
            t.PEa1Tt1(fac, unity);
            sumH += t.component(0,0);
            System.out.println(t.component(0,0));
        }

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


        if (doD2){
//            r1Hr2(Vector rij, Vector r1ij, Vector r2ij, double du, double d2u) { //r1.H.r2
//            double rij2 = rij.squared();

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
            Vector[] rdot = mapVel(drij);

            // dr1Hdr1 = 0 , dr2Hdr2 = 3 , dr3Hdr3 = 5
            // dr1Hdr2 = 1 , dr1Hdr3 = 2 , dr2Hdr3 = 4
            // dr4Hdr4 = 6 , dr5Hdr5 = 7 , dr6Hdr6 = 8
            rdotHrdot[0] += r1Hr2(rij, rdot[0], rdot[0], du, d2u);//C11
            rdotHrdot[3] += r1Hr2(rij, rdot[1], rdot[1], du, d2u);//C22
            rdotHrdot[5] += r1Hr2(rij, rdot[2], rdot[2], du, d2u);//C33
//
            rdotHrdot[1] += r1Hr2(rij, rdot[0], rdot[1], du, d2u);//C12
            rdotHrdot[2] += r1Hr2(rij, rdot[0], rdot[2], du, d2u);//C13
            rdotHrdot[4] += r1Hr2(rij, rdot[1], rdot[2], du, d2u);//C23
//
            rdotHrdot[6] += r1Hr2(rij, rdot[3], rdot[3], du, d2u);//C44
            rdotHrdot[7] += r1Hr2(rij, rdot[4], rdot[4], du, d2u);//C55
            rdotHrdot[8] += r1Hr2(rij, rdot[5], rdot[5], du, d2u);//C66
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
//            rx1Hdr1 += r1Hr2(rij, rx1, rdot[0], du, d2u);
//            ry2Hdr2 += r1Hr2(rij, ry2, rdot[1], du, d2u);
//            rz3Hdr3 += r1Hr2(rij, rz3, rdot[2], du, d2u);
//
//            rx1Hdr2 += r1Hr2(rij, rx1, rdot[1], du, d2u);
//            rx1Hdr3 += r1Hr2(rij, rx1, rdot[2], du, d2u);
//            ry2Hdr3 += r1Hr2(rij, ry2, rdot[2], du, d2u);
//
//            ry2Hdr1 += r1Hr2(rij, ry2, rdot[0], du, d2u);
//            rz3Hdr1 += r1Hr2(rij, rz3, rdot[0], du, d2u);
//            rz3Hdr2 += r1Hr2(rij, rz3, rdot[1], du, d2u);
//
//            ry3z2Hdr4 += r1Hr2(rij, ry3z2, rdot[3], du, d2u);
//            rx3z1Hdr5 += r1Hr2(rij, rx3z1, rdot[4], du, d2u);
//            rx2y1Hdr6 += r1Hr2(rij, rx2y1, rdot[5], du, d2u);
//
//            drHdr1 += r1Hr2(rij, drij, rdot[0], du, d2u);
//            drHdr2 += r1Hr2(rij, drij, rdot[1], du, d2u);
//            drHdr3 += r1Hr2(rij, drij, rdot[2], du, d2u);
//
//            rx1Hdr += r1Hr2(rij, rx1, drij, du, d2u);
//            ry2Hdr += r1Hr2(rij, ry2, drij, du, d2u);
//            rz3Hdr += r1Hr2(rij, rz3, drij, du, d2u);
        }
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // Add whatever you need to do with the hessian: elastic!

        boolean debug = false;
        if(debug && doD2 && (i==0 && j==nn) || (i==nn && j==0)){
            sumH += Hij.component(0,0);
        }


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


        Vector[] ridot = mapVel(dri);
        Vector[] rjdot = mapVel(drj);

//        dr1Hdr1 , dr2Hdr2 , dr3Hdr3
//        dr1Hdr2 , dr1Hdr3 , dr2Hdr3
//        dr4Hdr4 , dr5Hdr5 , dr6Hdr6
        int indCount=0;
        //00 , 01 , 02 , 11 , 12 , 22 , 33 , 44 , 55
        //0    1    2    3    4    5    6    7    8
        // dr1Hdr1 = 0 , dr2Hdr2 = 3 , dr3Hdr3 = 5
        // dr1Hdr2 = 1 , dr1Hdr3 = 2 , dr2Hdr3 = 4
        // dr4Hdr4 = 6 , dr5Hdr5 = 7 , dr6Hdr6 = 8
        for (int l=0; l<6; l++) {
            for (int ll=l; ll<6; ll++) {
                if((l <= 2 && ll > 2) || (l > 2 && ll > l)) continue;
                tmpV.E(rjdot[ll]);
                Hij.transform(tmpV); //Hij.drj
                rdotHrdot[indCount] += 2.0*ridot[l].dot(tmpV);
                tmpV.E(ridot[ll]);
                Hij.transform(tmpV);
                rdotHrdot[indCount] -= ridot[l].dot(tmpV);
                tmpV.E(rjdot[ll]);
                Hij.transform(tmpV);
                rdotHrdot[indCount] -= rjdot[l].dot(tmpV);
                indCount++;
            }
        }
        //dr1Hdr1
//        tmpV.E(rjdot[0]);
//        Hij.transform(tmpV); //Hij.drj
//        dr1Hdr1 += 2.0*ridot[0].dot(tmpV);
//        tmpV.E(ridot[0]);
//        Hij.transform(tmpV);
//        dr1Hdr1 -= ridot[0].dot(tmpV);
//        tmpV.E(rjdot[0]);
//        Hij.transform(tmpV);
//        dr1Hdr1 -= rjdot[0].dot(tmpV);
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

        rdot[3].setX(1, dr.getX(2));//rdot4: yz
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
        double r1Hr2 = du / rij2 * r1ij.dot(r2ij) + (d2u - du) / rij2 / rij2 * r1ij.dot(rij) * (r2ij.dot(rij));
        return r1Hr2;
    }

    protected double r1xH2abr2y(Vector rij, int a, int b, double x1, double y2, double du, double d2u) {
        double rij2 = rij.squared();
        double r1xH2abr2y = (d2u - du) / rij2 / rij2 * x1 * y2 * rij.getX(a) * rij.getX(b);
        return r1xH2abr2y;
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