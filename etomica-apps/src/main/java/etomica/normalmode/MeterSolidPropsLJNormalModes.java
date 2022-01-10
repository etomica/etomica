/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.lang.* ;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.math.Complex;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterSolidPropsLJNormalModes implements IDataSource {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected final DataSourceScalar meterPE;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<Vector> forceManager;
    protected  PotentialCalculationForceSum pcForceSum;
    private final PotentialCalculation pcSolidProps;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double ULat, PLat;
    protected double volume, density;
    protected int nMol, nRowsA, nColumnsA, nCells, nBasis, kdof;
    protected double dP, f1;
    protected  CoordinateDefinition.BasisCell[] cells0;
    protected double[][] A;
    protected Vector[] posCells;
    protected  int[] kDeg;
    protected Vector[] kVectors;
    protected double[][] gruneisen, gruneisen11, gruneisen22, gruneisen33;
    protected Complex[][][] evecs, evecs11, evecs22, evecs33;
    protected Complex[][][] constants_modal, constants_modal11, constants_modal22, constants_modal33;
    protected double[][]  force_real;
    protected Complex[][] force_complex;
    protected double[][]  dr_real;
    protected Complex[][] dr_complex;
    protected final double temperature, fV;
    protected double Pqh, Pqh11, Pqh22, Pqh33;
    protected Complex[] Qk_complex, Qk11_complex, Qk22_complex, Qk33_complex;
    protected Complex[] fk_complex, fk11_complex, fk22_complex, fk33_complex;


    public MeterSolidPropsLJNormalModes(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, CoordinateDefinition.BasisCell[] cells0, double temperature, double dP, double ULat, double PLat,
           double[][] A, double[][] B, double [][] C, double [][] D , double[][] C_modal, double[][] C_modal11, double[][] C_modal22, double[][] C_modal33) {
        int nData = 11;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        this.cells0 = cells0;
        id = new IteratorDirective();
        pcSolidProps = new PotentialCalculationLJSP(space,coordinateDefinition.getBox(),coordinateDefinition,temperature, new double[9]);
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<>(a -> space.makeVector(), coordinateDefinition.getBox());
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        this.temperature = temperature;
        this.dP = dP;
        this.ULat = ULat;
        this.PLat = PLat;
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().size();
        density = nMol/volume;

        nCells = cells0.length;
        nBasis = cells0[0].molecules.size();
        kdof = 3*nBasis;
        nRowsA = A.length;
        nColumnsA = A[0].length;
        kDeg = new int[nRowsA];
        gruneisen   = new double[nRowsA][kdof];
        gruneisen11 = new double[nRowsA][kdof];
        gruneisen22 = new double[nRowsA][kdof];
        gruneisen33 = new double[nRowsA][kdof];
        kVectors = new Vector[nRowsA];
        evecs   = new Complex[nRowsA][kdof][kdof];
        evecs11 = new Complex[nRowsA][kdof][kdof];
        evecs22 = new Complex[nRowsA][kdof][kdof];
        evecs33 = new Complex[nRowsA][kdof][kdof];
        constants_modal = new Complex[nRowsA][kdof][kdof];
        constants_modal11 = new Complex[nRowsA][kdof][kdof];
        constants_modal22 = new Complex[nRowsA][kdof][kdof];
        constants_modal33 = new Complex[nRowsA][kdof][kdof];
        double tmp_r;
        double tmp_i;
        double sum_gruneisen = 0;
        double sum_gruneisen11 = 0;
        double sum_gruneisen22 = 0;
        double sum_gruneisen33 = 0;

        Complex sumv, sum1, sum2, sum3;

        sumv = new Complex(0,0);

        for(int k=0; k<nRowsA; k++){
        	kDeg[k] = (int) A[k][0];
        	kVectors[k] = space.makeVector();
        	kVectors[k].setX(0, A[k][1]);
        	kVectors[k].setX(1, A[k][2]);
        	kVectors[k].setX(2, A[k][3]);
            for (int i=0; i<kdof; i++){
                gruneisen[k][i]   = A[k][4+i*(2*kdof+1)];
                gruneisen11[k][i] = B[k][4+i*(2*kdof+1)];
                gruneisen22[k][i] = C[k][4+i*(2*kdof+1)];
                gruneisen33[k][i] = D[k][4+i*(2*kdof+1)];
                sum_gruneisen   += kDeg[k]*gruneisen[k][i];
                sum_gruneisen11 += kDeg[k]*gruneisen11[k][i];
                sum_gruneisen22 += kDeg[k]*gruneisen22[k][i];
                sum_gruneisen33 += kDeg[k]*gruneisen33[k][i];
                for (int j=0; j<kdof; j++){
                    tmp_r = A[k][5 + i*(2*kdof+1) + 2*j];
                    tmp_i = A[k][5 + i*(2*kdof+1) + 2*j + 1];
                    evecs[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = B[k][5 + i*(2*kdof+1) + 2*j];
                    tmp_i = B[k][5 + i*(2*kdof+1) + 2*j + 1];
                    evecs11[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = C[k][5 + i*(2*kdof+1) + 2*j];
                    tmp_i = C[k][5 + i*(2*kdof+1) + 2*j + 1];
                    evecs22[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = D[k][5 + i*(2*kdof+1) + 2*j];
                    tmp_i = D[k][5 + i*(2*kdof+1) + 2*j + 1];
                    evecs33[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = C_modal[k][i*(2*kdof) + 2*j];
                    tmp_i = C_modal[k][i*(2*kdof) + 2*j + 1];
                    constants_modal[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = C_modal11[k][i*(2*kdof) + 2*j];
                    tmp_i = C_modal11[k][i*(2*kdof) + 2*j + 1];
                    constants_modal11[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = C_modal22[k][i*(2*kdof) + 2*j];
                    tmp_i = C_modal22[k][i*(2*kdof) + 2*j + 1];
                    constants_modal22[k][i][j] = new Complex(tmp_r, tmp_i);

                    tmp_r = C_modal33[k][i*(2*kdof) + 2*j];
                    tmp_i = C_modal33[k][i*(2*kdof) + 2*j + 1];
                    constants_modal33[k][i][j] = new Complex(tmp_r, tmp_i);
                }//j
            }//i


//            System.out.println("");
//            System.out.println(" k = " + "  " + (k+1));
//            for (int i=0; i<kdof; i++) {
//                for (int j = 0; j < kdof; j++) {
//                    sumv = new Complex(0,0);
//                    sum1 = new Complex(0,0);
//                    sum2 = new Complex(0,0);
//                    sum3 = new Complex(0,0);
//                    int M = 1;
//                    for (int m = 0; m < kdof; m++) {
////                        if(m==9||m==10||m==11){
////                            sumv  = sumv.plus(evecs[k][m][i].times(evecs[k][m][j].conjugate()).times(new Complex(gruneisen[k][m],0)));
////                            sum1  = sum1.plus(evecs11[k][m][i].times(evecs11[k][m][j].conjugate()).times(new Complex(gruneisen11[k][m],0)));
////                            sum2  = sum2.plus(evecs22[k][m][i].times(evecs22[k][m][j].conjugate()).times(new Complex(gruneisen22[k][m],0)));
////                            sum3  = sum3.plus(evecs33[k][m][i].times(evecs33[k][m][j].conjugate()).times(new Complex(gruneisen33[k][m],0)));
//                            sum1  = sum1.plus((evecs11[k][i][m].conjugate()).times(evecs11[k][j][m]));
////                        }
//                    }
////                    double xv = ((sumv.times(new Complex(3*volume,0))).minus(sum1.plus(sum2).plus(sum3))).modulus();
//
//                    if(sum1.modulus()<1e-10){
//                        System.out.print(" 0 ");
//                    }else if(Math.abs(sum1.modulus()-1.0)<1e-10) {
//                        System.out.print(" 1 ");
//                    }else{
//                        System.out.print(sum1.modulus());
//                        // System.out.println((i+1) + " "+ (j+1) + "  <<->>  " + sum1.modulus());
//                    }
//                }//j
//                System.out.println();
////                System.out.println(3*volume*gruneisen[k][i] - (gruneisen11[k][i]+gruneisen22[k][i]+gruneisen33[k][i]));
//            }//i

        }//k


//        System.exit(0);


        force_real    = new double[nCells][kdof];
        force_complex = new Complex[nRowsA][kdof];
        fk_complex    = new Complex[kdof];
        fk11_complex  = new Complex[kdof];
        fk22_complex  = new Complex[kdof];
        fk33_complex  = new Complex[kdof];

        dr_real       = new double[nCells][kdof];
        dr_complex    = new Complex[nRowsA][kdof];
        Qk_complex    = new Complex[kdof];
        Qk11_complex  = new Complex[kdof];
        Qk22_complex  = new Complex[kdof];
        Qk33_complex  = new Complex[kdof];

        Pqh   = temperature*(sum_gruneisen+1/volume);
        Pqh11 = temperature*(sum_gruneisen11/volume+1/volume);
        Pqh22 = temperature*(sum_gruneisen22/volume+1/volume);
        Pqh33 = temperature*(sum_gruneisen33/volume+1/volume);
        fV = (Pqh/temperature-density)/3.0/(nMol-1.0);
        System.out.println(" Pqh: " + Pqh + "  " + Pqh11 + "  " + Pqh22 + "   " +  Pqh33 );
        System.out.println(volume);
    }

    public IData getData() {
        Box box = coordinateDefinition.getBox();

        double[] sum;
        ((PotentialCalculationLJSP)pcSolidProps).reset();
        potentialMaster.calculate(box, id, pcSolidProps);
        sum = ((PotentialCalculationLJSP)pcSolidProps).getSum();

        pcForceSum.reset();
        potentialMaster.calculate(box, id, pcForceSum);
        double[] x   = data.getData();
        double fr    = sum[0];
        double fdr   = sum[1];
        double fxrx  = sum[2];
        double fyry  = sum[3];
        double fzrz  = sum[4];
        double fxdrx = sum[5];
        double fydry = sum[6];
        double fzdrz = sum[7];
        double dU    = meterPE.getDataAsScalar() - ULat;

        IAtom atomBasis;
        int iBasis;
        Vector  dri = space.makeVector();
        Vector ri_lat, ri;
        for(int l=0; l<nCells; l++){
            for(int i=0; i<nBasis; i++){
                iBasis  = cells0[l].molecules.get(i).getIndex();
                atomBasis = coordinateDefinition.getBox().getLeafList().getAtoms().get(iBasis);
                ri_lat = coordinateDefinition.getLatticePosition(atomBasis);
                ri = atomBasis.getPosition();
                dri.Ev1Mv2(ri, ri_lat);
                box.getBoundary().nearestImage(dri);
                for(int j=0; j<3; j++) {
                    force_real[l][3*i+j] = forceManager.getAgent(atomBasis).getX(j);
                    dr_real[l][3*i+j] = dri.getX(j);
                }
            }
        } //l

        // Fourier Transform
        force_complex = realToReciprocal(force_real);
        dr_complex    = realToReciprocal(dr_real);

        double sum_real = 0;
        double sum_real11 = 0;
        double sum_real22 = 0;
        double sum_real33 = 0;

        for(int k=0; k<nRowsA; k++) { // wavevector
            for (int i=0; i<kdof; i++) { // mode
                if (Math.abs(gruneisen[k][i]) < 1e-9)  continue; // skip COM
                fk_complex[i] = new Complex(0,0);
                fk11_complex[i] = new Complex(0,0);
                fk22_complex[i] = new Complex(0,0);
                fk33_complex[i] = new Complex(0,0);

                Qk_complex[i] = new Complex(0,0);
                Qk11_complex[i] = new Complex(0,0);
                Qk22_complex[i] = new Complex(0,0);
                Qk33_complex[i] = new Complex(0,0);

                for (int j=0; j<kdof; j++) {
                    //Forces
                    fk_complex[i]   = fk_complex[i].plus(  force_complex[k][j].times(evecs[k][i][j].conjugate()));
                    fk11_complex[i] = fk11_complex[i].plus(force_complex[k][j].times(evecs11[k][i][j].conjugate()));
                    fk22_complex[i] = fk22_complex[i].plus(force_complex[k][j].times(evecs22[k][i][j].conjugate()));
                    fk33_complex[i] = fk33_complex[i].plus(force_complex[k][j].times(evecs33[k][i][j].conjugate()));
                    //Coordinates
                    Qk_complex[i]   = Qk_complex[i].plus(  dr_complex[k][j].times(evecs[k][i][j].conjugate()));
                    Qk11_complex[i] = Qk11_complex[i].plus(dr_complex[k][j].times(evecs11[k][i][j].conjugate()));
                    Qk22_complex[i] = Qk22_complex[i].plus(dr_complex[k][j].times(evecs22[k][i][j].conjugate()));
                    Qk33_complex[i] = Qk33_complex[i].plus(dr_complex[k][j].times(evecs33[k][i][j].conjugate()));
                }

                    sum_real   += kDeg[k]*gruneisen[k][i]  *((fk_complex[i].conjugate()).times(Qk_complex[i])).real();
                    sum_real11 += kDeg[k]*gruneisen11[k][i]*((fk11_complex[i].conjugate()).times(Qk11_complex[i])).real();
                    sum_real22 += kDeg[k]*gruneisen22[k][i]*((fk22_complex[i].conjugate()).times(Qk22_complex[i])).real();
                    sum_real33 += kDeg[k]*gruneisen33[k][i]*((fk33_complex[i].conjugate()).times(Qk33_complex[i])).real();

            }//i

            // Eigenvector derivatives
            if(true){
                for (int i=0; i<kdof; i++) {
                    for (int j=0; j<kdof; j++) {
                        if (Math.abs(gruneisen[k][i]) == 0 || Math.abs(gruneisen[k][j]) == 0) continue; // skip COM
                        sum_real += kDeg[k]*((fk_complex[i].conjugate()).times(constants_modal[k][i][j]).times(Qk_complex[j])).real();
                    }
                }

                for (int i=0; i<kdof; i++) {
                    for (int j=0; j<kdof; j++) {
                        if (Math.abs(gruneisen11[k][i]) == 0 || Math.abs(gruneisen11[k][j]) == 0) continue;
                        sum_real11 += kDeg[k]*((fk11_complex[i].conjugate()).times(constants_modal11[k][i][j]).times(Qk11_complex[j])).real();
                    }
                }
                for (int i=0; i<kdof; i++) {
                    for (int j=0; j<kdof; j++) {
                        if (Math.abs(gruneisen22[k][i]) == 0 || Math.abs(gruneisen22[k][j]) == 0)  continue;
                        sum_real22 += kDeg[k]*((fk22_complex[i].conjugate()).times(constants_modal22[k][i][j]).times(Qk22_complex[j])).real();
                    }
                }
                for (int i=0; i<kdof; i++) {
                    for (int j=0; j<kdof; j++) {
                        if (Math.abs(gruneisen33[k][i]) == 0 || Math.abs(gruneisen33[k][j]) == 0)  continue;
                        sum_real33 += kDeg[k]*((fk33_complex[i].conjugate()).times(constants_modal33[k][i][j]).times(Qk33_complex[j])).real();
                    }
                }
            }

        }//k

        sum_real /= nCells;
        sum_real11 /= nCells;
        sum_real22 /= nCells;
        sum_real33 /= nCells;
// P
        double Pvir  = fr/3.0/volume;
        double Pvir11  = fxrx/volume;
        double Pvir22  = fyry/volume;
        double Pvir33  = fzrz/volume;

        double Pconv = Pvir + density*temperature;
        double Phma  = Pqh + Pvir + fV*fdr;
        double Pnm   = Pqh + Pvir - fdr/3.0/volume + sum_real ;

        double Pconv11 = Pvir11 + density*temperature;
        double Pconv22 = Pvir22 + density*temperature;
        double Pconv33 = Pvir33 + density*temperature;

        double Pnm11 = Pqh11 + Pvir11 - fxdrx/volume + sum_real11/volume;
        double Pnm22 = Pqh22 + Pvir22 - fydry/volume + sum_real22/volume;
        double Pnm33 = Pqh33 + Pvir33 - fzdrz/volume + sum_real33/volume;

        x[0] = Pconv;
        x[1] = Phma;
        x[2] = Pnm;
        x[3] = Pconv11;
        x[4] = Pconv22;
        x[5] = Pconv33;
        x[6] = Pnm11;
        x[7] = Pnm22;
        x[8] = Pnm33;
        x[9] = (Pconv11+Pconv22+Pconv33)/3.0;
        x[10]= (Pnm11+Pnm22+Pnm33)/3.0;
        return data;
    }


    private Complex[][] realToReciprocal(double[][] vec_real) {
        Complex[][] vec_complex = new Complex[nRowsA][kdof];

        for(int k=0; k<nRowsA; k++) {
            for (int i=0; i<kdof; i++) {
                double sum_real = 0;
                double sum_imag = 0;
                for (int l=0; l<nCells; l++) {
                    double kdotr = kVectors[k].dot(cells0[l].cellPosition);
                    sum_real += vec_real[l][i]*Math.cos(kdotr);
                    sum_imag -= vec_real[l][i]*Math.sin(kdotr);
                } //l
                vec_complex[k][i] = new Complex(sum_real, sum_imag);
            }
        }
        return vec_complex;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}