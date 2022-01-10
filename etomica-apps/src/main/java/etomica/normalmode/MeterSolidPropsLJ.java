/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidPropsLJ implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected final DataSourceScalar meterPE;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double volume, density;
    protected int nMol, nCells;
    protected double gV, gVV, gx1, gy1, gy4, gx11, gy11, gx12, gz12, gx44, gy44;
    protected double Ushift, Pshift;

    protected final double temperature;
    private final PotentialCalculation pcSolidProps;
    protected  PotentialCalculationForceSum pcForceSum;

    public MeterSolidPropsLJ(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double temperature, double[] elasticParams) {
        int nData = 41;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), coordinateDefinition.getBox());
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        this.temperature = temperature;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().size();
        density = nMol/volume;

        gV   = elasticParams[0];  gVV  = elasticParams[1];
        gx1  = elasticParams[2];  gy1  = elasticParams[3];   gy4 = elasticParams[4];
        gx11 = elasticParams[5];  gy11 = elasticParams[6];  gx44 = elasticParams[7];
        gy44 = elasticParams[8];  gx12 = elasticParams[9];  gz12 = elasticParams[10];

        this.Ushift = 0;
        this.Pshift = 0;
        pcSolidProps = new PotentialCalculationLJSP(space,coordinateDefinition.getBox(),coordinateDefinition,temperature, elasticParams);
    }


    public IData getData() {
        Box box = coordinateDefinition.getBox();
        double[] sum;
        ((PotentialCalculationLJSP)pcSolidProps).reset();
        potentialMaster.calculate(box, id, pcSolidProps);
        sum = ((PotentialCalculationLJSP)pcSolidProps).getSum();
        potentialMaster.calculate(box, id, pcForceSum);
        double[] x = data.getData();

        double Fr      = sum[0];
        double Fdr     = sum[1];
        double rPhir   = sum[2];
        double drPhidr = sum[3];
        double rPhidr  = sum[4];

        double Fdr1 = sum[5];
        double Fdr2 = sum[6];
        double Fdr3 = sum[7];
        double Fdr4 = sum[8];
        double Fdr5 = sum[9];
        double Fdr6 = sum[10];

        double Fdr11 = sum[11];
        double Fdr22 = sum[12];
        double Fdr33 = sum[13];
        double Fdr12 = sum[14];
        double Fdr13 = sum[15];
        double Fdr23 = sum[16];
        double Fdr44 = sum[17];
        double Fdr55 = sum[18];
        double Fdr66 = sum[19];

        double Fxrx = sum[20];
        double Fyry = sum[21];
        double Fzrz = sum[22];

        double Fxry = sum[23];
        double Fxrz = sum[24];
        double Fyrz = sum[25];

        double x_Phixx_x = sum[26];
        double y_Phiyy_y = sum[27];
        double z_Phizz_z = sum[28];

        double x_Phixy_y = sum[29];
        double x_Phixz_z = sum[30];
        double y_Phiyz_z = sum[31];

        double x_Phiyx_y = sum[32];
        double x_Phizx_z = sum[33];
        double y_Phizy_z = sum[34];

        double dr1_Phi_dr1 = sum[35];
        double dr2_Phi_dr2 = sum[36];
        double dr3_Phi_dr3 = sum[37];

        double dr1_Phi_dr2 = sum[38];
        double dr1_Phi_dr3 = sum[39];
        double dr2_Phi_dr3 = sum[40];

        double dr4_Phi_dr4 = sum[41];
        double dr5_Phi_dr5 = sum[42];
        double dr6_Phi_dr6 = sum[43];

        double rx1_Phi_dr1 = sum[44];
        double ry2_Phi_dr2 = sum[45];
        double rz3_Phi_dr3 = sum[46];

        double rx1_Phi_dr2 = sum[47];
        double rx1_Phi_dr3 = sum[48];
        double ry2_Phi_dr3 = sum[49];

        double ry2_Phi_dr1 = sum[50];
        double rz3_Phi_dr1 = sum[51];
        double rz3_Phi_dr2 = sum[52];

        double ry3z2_Phi_dr4 = sum[53];
        double rx3z1_Phi_dr5 = sum[54];
        double rx2y1_Phi_dr6 = sum[55];
        //b_mn
        double Fdr1G = sum[56];
        double Fdr2G = sum[57];
        double Fdr3G = sum[58];

        double drPhidr1 = sum[59];
        double drPhidr2 = sum[60];
        double drPhidr3 = sum[61];

        double rx1_Phi_dr = sum[62];
        double ry2_Phi_dr = sum[63];
        double rz3_Phi_dr = sum[64];

        double U = meterPE.getDataAsScalar();
        //U
        x[0] = U - Ushift;
        x[1] = U + 1.0/2.0*Fdr + 3.0/2.0*(nMol-1.0)*temperature - Ushift;
        //P: shiftP=P
        x[2] = 1.0/3.0/volume*Fr + density*temperature - this.Pshift;
        double fV = gV-1.0/3.0;//gV/volume-1.0/3.0/volume;
        x[3] = 1.0/3.0/volume*Fr + fV/volume*Fdr + 3*(nMol-1)*temperature/volume*gV + temperature/volume - this.Pshift;

        //P1, P2, P3: shiftP = -P
        x[4] = -1.0/volume*Fxrx - density*temperature + this.Pshift;//P1
        x[5] = -1.0/volume*(Fxrx+Fdr1 + (nMol-1)*temperature*(gx1+2.0*gy1) + temperature) + this.Pshift;
        x[6] = -1.0/volume*Fyry - density*temperature + this.Pshift;//P2
        x[7] = -1.0/volume*(Fyry+Fdr2 + (nMol-1)*temperature*(gx1+2.0*gy1) + temperature) + this.Pshift;
        x[8] = -1.0/volume*Fzrz - density*temperature + this.Pshift;//P3
        x[9] = -1.0/volume*(Fzrz+Fdr3 + (nMol-1)*temperature*(gx1+2.0*gy1) + temperature) + this.Pshift;
        //P4, P5, P5: no shift (<P4>=<P5>=<P6>=0)
        x[10] = -1.0/volume*Fyrz;//P4
        x[11] = -1.0/volume*(Fyrz+Fdr4);
        x[12] = -1.0/volume*Fxrz;//P5
        x[13] = -1.0/volume*(Fxrz+Fdr5);
        x[14] = -1.0/volume*Fxry;//P6
        x[15] = -1.0/volume*(Fxry+Fdr6);
        //Cv
        x[16] = -1.0/4.0/temperature*(Fdr+drPhidr) + 3.0/2.0*(nMol-1.0);
        //B
        double d2UdeV_vir = 1.0/9.0*rPhir + 2.0/9.0*Fr;
        x[17] = 1.0/volume*d2UdeV_vir + density*temperature;
        double hVV = gVV+gV*gV+2.0/9.0;
        x[18] = 1.0/volume*(d2UdeV_vir + fV*fV*drPhidr - hVV*Fdr + 2.0/3.0*fV*rPhidr) - 3.0*(nMol-1)*temperature*gVV/volume + temperature/volume;
        //C11, C22, C33
        x[19] = 1.0/volume*x_Phixx_x + 2*density*temperature; //C11
        x[20] = 1.0/volume*(x_Phixx_x - Fdr11+dr1_Phi_dr1 + 2*rx1_Phi_dr1 - (nMol-1.0)*temperature*(gx11 +2*gy11)) + 2*temperature/volume;
        x[21] = 1.0/volume*y_Phiyy_y  + 2*density*temperature;//C22
        x[22] = 1.0/volume*(y_Phiyy_y - Fdr22+dr2_Phi_dr2 + 2*ry2_Phi_dr2) - (nMol-1.0)/volume*temperature*(gx11 +2*gy11) + 2*temperature/volume;
        x[23] = 1.0/volume*z_Phizz_z + 2*density*temperature;//C33
        x[24] = 1.0/volume*(z_Phizz_z - Fdr33+dr3_Phi_dr3 + 2*rz3_Phi_dr3) - (nMol-1.0)/volume*temperature*(gx11 +2*gy11) + 2*temperature/volume;
        //C12, C13, C23
        x[25] = 1.0/volume*x_Phixy_y;//C12
        x[26] = 1.0/volume*(x_Phixy_y - Fdr12 + dr1_Phi_dr2 + rx1_Phi_dr2 + ry2_Phi_dr1 - (nMol-1.0)*temperature*(gz12+2*gx12));
        x[27] = 1.0/volume*x_Phixz_z;//C13
        x[28] = 1.0/volume*(x_Phixz_z - Fdr13 + dr1_Phi_dr3 + rx1_Phi_dr3 + rz3_Phi_dr1 - (nMol-1.0)*temperature*(gz12+2*gx12));
        x[29] = 1.0/volume*y_Phiyz_z;//C23
        x[30] = 1.0/volume*(y_Phiyz_z - Fdr23 + dr2_Phi_dr3 + ry2_Phi_dr3 + rz3_Phi_dr2 - (nMol-1.0)*temperature*(gz12+2*gx12));
        //C44, C55, C66
        x[31] = 1.0/volume*y_Phizy_z + density*temperature;//C44
        x[32] = 1.0/volume*(y_Phizy_z - Fdr44 + dr4_Phi_dr4 + ry3z2_Phi_dr4 - (nMol-1)*temperature*(gx44+2*gy44)) + temperature/volume;
        x[33] = 1.0/volume*x_Phizx_z + density*temperature;//C55
        x[34] = 1.0/volume*(x_Phizx_z - Fdr55 + dr5_Phi_dr5 + rx3z1_Phi_dr5 - (nMol-1)*temperature*(gx44+2*gy44)) + temperature/volume;
        x[35] = 1.0/volume*x_Phiyx_y + density*temperature;//C66
        x[36] = 1.0/volume*(x_Phiyx_y - Fdr66 + dr6_Phi_dr6 + rx2y1_Phi_dr6 - (nMol-1)*temperature*(gx44+2*gy44)) + temperature/volume;
        //b_V
        x[37] = 1.0/2.0/volume/temperature*(gV*Fdr-fV*drPhidr-1.0/3.0*rPhidr) + 3.0*(nMol-1)*gV/volume + 1.0/volume; //bV_hma
        //b_mn
        x[38] = 1.0/2.0/volume/temperature*(-Fdr1G+drPhidr1+rx1_Phi_dr) - (nMol-1)/volume*(gx1+2*gy1)-1/volume;
        x[39] = 1.0/2.0/volume/temperature*(-Fdr2G+drPhidr2+ry2_Phi_dr) - (nMol-1)/volume*(gx1+2*gy1)-1/volume;
        x[40] = 1.0/2.0/volume/temperature*(-Fdr3G+drPhidr3+rz3_Phi_dr) - (nMol-1)/volume*(gx1+2*gy1)-1/volume;

        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public void setShift(double Ushift , double Pshift){
        this.Ushift = Ushift;
        this.Pshift = Pshift;
    }

}