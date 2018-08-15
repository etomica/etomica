/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.potential;

import etomica.exception.MethodNotImplementedException;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

/**
 * pair potential for Helium Przybytek et al. (2017) Phys. Rev. Lett. 119, 123401.
 *
 * @author Navneeth Gokul
 */
public class P2HePCJS extends Potential2SoftSpherical {

    protected final int M_BO;
    protected final int I0_BO, I1_BO;
    protected final int N0_BO, N1_BO;
    protected final double[] a_BO;
    protected final double[][] P_BO;
    protected final double zeta_BO;
    protected final double[] C_BO;
    protected final int M_AD;
    protected final int I0_AD, I1_AD;
    protected final int N0_AD, N1_AD;
    protected final double[] a_AD;
    protected final double[][] P_AD;
    protected final double zeta_AD;
    protected final double[] C_AD;
    protected final int M_REL;
    protected final int I0_REL, I1_REL;
    protected final int N0_REL, N1_REL;
    protected final double[] a_REL;
    protected final double[][] P_REL;
    protected final double zeta_REL;
    protected final double[] C_REL;
    protected final int M_QED;
    protected final int I0_QED, I1_QED;
    protected final int N0_QED, N1_QED;
    protected final double[] a_QED;
    protected final double[][] P_QED;
    protected final double zeta_QED;
    protected final double[] C_QED;

    public P2HePCJS(Space space) {
        super(space);

        M_BO = 3;
        I0_BO = -1; I1_BO = 2;
        N0_BO =  6; N1_BO = 16;

        a_BO = new double[]{2.19076900924970,
                            3.91558362097102,
                            8.87314362466591};

        // P_BO (I0_BO:I1_BO, M_BO) = P_BO(-1:2 , 3)
        // shape=( I1_BO-I0_BO+1, M_BO ) = (2-(-1)+1 , 3) = (4,3)
        P_BO = new double[][]   {{2.69150966652494e1, -2.62457650950449e2, 2.39542554285199e2},
                                {1.74015968871153e1, 4.57300508401525e2, 6.73223355754347e2},
                                {-3.26566208189471e0, 6.66624947516144e2, 8.62651807552327e2},
                                {1.46700433938961e-1, -1.01541164656195e2, 7.86212667513655e2}};

        zeta_BO = 4.45565354034820;

        // C_BO(N0_BO:N1_BO) = C_BO(6,16)
        C_BO = new double[] {0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1.460977837725,
                            0.0,
                            14.11785737,
                            0.0,
                            183.691075,
                            -76.72571,
                            3.372e3,
                            -3808.3254,
                            8.534e4,
                            -1.707e5,
                            2.86e6};

        M_AD = 3;
        I0_AD = 0; I1_AD = 2;
        N0_AD =  6; N1_AD = 10;

        a_AD = new double[]{2.75067513975782,
                            5.36085266623527,
                            7.13731208196370};

        // P_AD ( I0_AD:I1_AD, M_AD ) = P_AD(0:2 , 3)
        // shape=( I1_AD-I0_AD+1, M_AD ) = (2-0 +1 , 3) = (3,3)
        P_AD = new double[][]   {{1.60034916931265e-2, -5.02322894521138e-1, -7.80222735840106e-2},
                                {-5.17915211148225e-3, -5.28178298205425e-2, 1.83457620987126e0},
                                {3.20021447656502e-3, 2.93999790052146e-1, 9.46016993696249e-1}};

        zeta_AD = 4.80971750762740;

        // C_AD(N0_AD:N1_AD) = C_AD(6,10)
        C_AD = new double[] {0,
                            0,
                            0,
                            0,
                            0,
                            0,
                            1.1445e-3,
                            0.0,
                            6.519e-3,
                            0.0,
                            6.68e-2};

        M_REL = 3;
        I0_REL = 0; I1_REL = 2;
        N0_REL =  4; N1_REL = 8;

        a_REL = new double[]{2.55442787544336,
                            4.26822206589066,
                            6.31801545306225};

        // P_REL( I0_REL:I1_REL, M_REL) = P_REL(0:2 , 3)
        // shape=( I1_REL-I0_REL+1,M_REL ) = (2-0 +1 , 3) = (3,3)
        P_REL = new double[][]  {{1.16259940182396e-2, -5.73711010584383e-1, 6.22980694403664e-1},
                                {-7.00805799572921e-3, 3.97167157865319e-1, 4.87089036456869e-1},
                                {7.83673396684808e-4, -1.10854901515699e-1, 6.81223617553767e-1}};

        zeta_REL = 5.00212547573053;

        // C_REL(N0_REL:N1_REL) = C_REL(4,8)
        C_REL = new double[]{0,
                            0,
                            0,
                            0,
                            -3.5322e-5,
                            0.0,
                            -3.434e-4,
                            0.0,
                            -3.55798027534830e-3};

        M_QED = 2;
        I0_QED = 0; I1_QED = 2;
        N0_QED = 3; N1_QED = 6;

        a_QED = new double[]{2.71385934118196,
                            4.85539104988208};

        // P_QED(I0_QED:I1_QED,M_QED) = P_QED(0:2 , 2)
        // shape=( I1_REL-I0_REL+1,M_REL ) = (2-0 +1 , 2) = (3,2)
        P_QED = new double[][]  {{1.34952154947189e-3, -1.96866479382540e-4},
                                {-1.02362538302740e-3, -9.61841561868819e-4},
                                {3.56862380020579e-4, -1.08014679908054e-3}};




        zeta_QED = 5.43735514285525;

        // C_QED(N0_QED:N1_QED) = C_QED(3,6)
        C_QED = new double[]{0,
                            0,
                            0,
                            5.772353e-7,
                            0.0,
                            1.377841e-6,
                            7.61187886972970e-5};

    }

    /**
     * @return energy.
     * @param r distance in Bohr
     * @param doRet boolean to compute energy including retardation.
     */
    private double V(double r, boolean doRet) {

        double elms1 = formShort(r,M_BO,I0_BO,I1_BO,a_BO,P_BO);
        double elms2 = formShort(r,M_AD,I0_AD,I1_AD,a_AD,P_AD);
        double elms3 = formShort(r,M_REL,I0_REL,I1_REL,a_REL,P_REL);
        double elms4 = formShort(r,M_QED,I0_QED,I1_QED,a_QED,P_QED);
        //System.out.println(elms1+" "+elms2+" "+elms3+" "+elms4);
        if(doRet){
            elms1 += formLong(r, 1, N0_BO, N1_BO, zeta_BO, C_BO);
            elms2 += formLong(r, 0, N0_AD, N1_AD, zeta_AD, C_AD);
            elms3 += formLong(r, 2, N0_REL, N1_REL, zeta_REL, C_REL);
            elms4 += formLong(r, 2, N0_QED, N1_QED, zeta_QED, C_QED);
            //System.out.println(elms1+" "+elms2+" "+elms3+" "+elms4);
        }
        else {
            elms1 += formLong(r, 0, N0_BO, N1_BO, zeta_BO, C_BO);
            elms2 += formLong(r, 0, N0_AD, N1_AD, zeta_AD, C_AD);
            elms3 += formLong(r, 0, N0_REL, N1_REL, zeta_REL, C_REL);
            elms4 += formLong(r, 0, N0_QED, N1_QED, zeta_QED, C_QED);
            //System.out.println(elms1+" "+elms2+" "+elms3+" "+elms4);
        }
        double V = elms1 + elms2 + elms3 + elms4;
        return V;
    }

    private double formShort(double r, int M, int I0, int I1,double[] a, double[][] P){

        //short-range formula
        //sum_{k=1}^{M} exp(-a(k) R) sum_{i=I0}^{I1} P(i,k) R^i

        double formShort = 0.0;
        for (int k = 0; k < M; k++){
            double term = P[I1-I0][k]; // subtracted I1 from i for appropriate index
            //if(I0==-1)System.out.println("k = " + k +" I1 = "+ I1 +" P[I1][k] = " + term);
            for (int i = I1-1; i >= I0; i--){
                term = term*r + P[i-I0][k]; // subtracted I0 from i for appropriate index
                //if(I0==-1)System.out.println("k = " + k + " i = "+ i +" P[i][k] = " + P[i-I0][k] + " term = " + term);
            }
            formShort += Math.exp(-a[k]*r)*term;
            //if(I0==-1)System.out.println("k = " + k +" fs = " + formShort);
        }
        if(I0 != 0) formShort = formShort*Math.pow(r,I0);

        return formShort;
    }

    private double formLong(double r, int ret_type, int N0, int N1,double zeta, double[] C) {
        //long-range formula
        //
        //       - sum_{n=N0}^{N1} C(n) dampTT(n,zeta R)/R^n
        //
        //       damping of the leading term
        // ret_type = 0  ->   1     + shiftTT = dampTT (no retardation)
        // ret_type = 1  ->   retar + shiftTT          (BO [C6] potential)
        // ret_type = 2  ->   0     + shiftTT          (REL [C4] and QED [C3] potentials)

        double iR = 1/r;
        //if(N1==16)System.out.println("iR = "+iR);
        double zetaR = zeta*r;
        //if(N1==16)System.out.println("zetaR = "+zetaR);
        double formLong;

        if( N0 < N1 ) {
            formLong = - C[N1]*dampTT(N1,zetaR);
            //if(N1==16)System.out.println("N1 = " + N1 +" C[N1] = " + C[N1] +" dampTT(N1,zetaR) = "+ dampTT(N1,zetaR) +" fl = " + formLong);
            for (int n = N1-1; n >= N0+1; n -- ) {
                if( C[n] != 0) {
                    formLong = formLong*iR - C[n]*dampTT(n,zetaR);
                    //if(N1==16)System.out.println("n = " + n + " C[n] = " + C[n] +" dampTT(n,zetaR) = "+ dampTT(n,zetaR) +" fl = " + formLong);
                }
                else {
                    formLong = formLong*iR;
                    //if(N1==16)System.out.println("else fl = " + formLong);
                }
            }
        }
        else {
            formLong = 0.0;
        }
        switch (ret_type){
            case 0:
                formLong = formLong*iR - C[N0]*dampTT(N0,zetaR);
                break;
            case 1:
                formLong = formLong*iR - C[N0]*(retar(r)+shftTT(N0,zetaR));
                //if(N1==16)System.out.println("N0 = " + N0 +" C[N0] = " + C[N0] +" retar(r) = "+ retar(r) +" shftTT(N0,zetaR) = "+ shftTT(N0,zetaR) +" case fl = " + formLong);
                break;
            case 2:
                formLong = formLong*iR - C[N0]*shftTT(N0,zetaR);
                break;
            default:
                throw new RuntimeException("!!!!!WHOPPS!!!!");
        }

        formLong = formLong*Math.pow(iR,N0);
        return formLong;
    }

    private double dampTT(int n, double x){
        //Tang-Toennies damping function
        //
        //         dampTT(n,x) = 1 - exp(-x) sum_{i=0}^{n} x^i/i!
        //
        //                     = exp(-x) sum_{i=n+1}^{infty} x^i/i!
        //

        double dampTT_MINlimit = Double.MIN_VALUE;

        double[] dampTT_MIDlimit = {0.10536052,
                                    0.53181161, 1.1020653, 1.7447696, 2.4325910,
                                    3.1518980 , 3.8947668, 4.6561182, 5.4324681,
                                    6.2213046 , 7.0207466, 7.8293420, 8.6459425,
                                    9.4696212 ,10.299617 ,11.135297 ,11.976127};

        double[] dampTT_MAXlimit =  {35.350506,
                                    39.040395, 42.189070, 45.049320, 47.719353,
                                    50.250606, 52.674332, 55.011322, 57.276297,
                                    59.480153, 61.631241, 63.736139, 65.800138,
                                    67.827579, 69.822072, 71.786664, 73.723952};

        double dampTT;
        double term,suma;
        if ( x < dampTT_MINlimit ) {
            dampTT = 0.0;
        }
        else if ( x < dampTT_MIDlimit[n]) {
            term = 1.0;
            for ( int i = 1; i <= n; i++) {
                term *= x/i;
                //if(n==16)System.out.println("1: i = " + i +" term = "+term);
            }
            suma = 0.0;
            for ( int i = n+1; i <= 2*n+17; i++) {
                term *= x/i;
                suma += term;
                //if(n==16)System.out.println("2: i = " + i +" suma = "+ suma + " term = "+term);
            }
            dampTT = Math.exp(-x)*suma;
            //if(n==16)System.out.println("3: suma = "+ suma + " dampTT = "+dampTT);
        }
        else if ( x < dampTT_MAXlimit[n]) {
            term = 1.0;
            suma = 1.0;
            for ( int i = 1; i <= n; i++) {
                term *= x/i;
                suma += term;
                //if(n==16)System.out.println("4: i = " + i +" suma = "+ suma + " term = "+term);
            }
            dampTT = 1.0 - Math.exp(-x)*suma;
            //if(n==16)System.out.println("5: suma = "+ suma + " dampTT = "+dampTT);
        }
        else {
            dampTT = 1.0;
        }
        return dampTT;
    }

    private double shftTT(int n, double x){
        //shifted Tang-Toennies damping function
        //
        //         shftTT(n,x) = dampTT(n,x) - 1
        //
        //             = - exp(-x) sum_{i=0}^{n} x^i/i!
        //

        double shftTT;
        double shftTT_MAXlimit = 709.782712893384; //WHAT IS THIS???
        double term,suma;

        if ( x < shftTT_MAXlimit ){
            term = 1.0;
            suma = 1.0;
            for ( int i = 1; i <= n; i++) {
                term *= x/i;
                suma += term;
            }
            shftTT = -Math.exp(-x)*suma;
        }
        else {
            shftTT = 0.0;
        }
        return shftTT;
    }

    private double retar(double r){
        //A(5)
        double[] A = new double[]   {6.169870091028090e-2,
                                    8.281954197590046e-4,
                                    2.936387935854446e-6,
                                    4.020288241131003e-9,
                                    2.948899823358023e-12,};
        //B(6)
        double[] B = new double[]   {6.169870091028090e-2,
                                    8.523723778811090e-4,
                                    4.032972820432622e-6,
                                    9.969788378653687e-9,
                                    1.224004819647673e-11,
                                    8.978131367598067e-15};

        double num,den;

        num = A[4];
        num = num*r + A[3];
        num = num*r + A[2];
        num = num*r + A[1];
        num = num*r + A[0];
        num = num*r + 1.0;

        den = B[5];
        den = den*r + B[4];
        den = den*r + B[3];
        den = den*r + B[2];
        den = den*r + B[1];
        den = den*r + B[0];
        den = den*r + 1.0;

        double retar = num/den;
        return retar;
    }

    /**
     * @return energy u.
     * @param r2 square of the distance in Angstrom^2.
     */
    public double u(double r2) {
        if (r2 < 1e-8) return Double.MAX_VALUE;

        double r = Math.sqrt(r2); //distance in Angstrom
        double rBohr = BohrRadius.UNIT.fromSim(r); //distance in Bohr

        double uHartree =  V(rBohr,true); //energy in Hartree. Always account for retardation

        return Hartree.UNIT.toSim(uHartree);
    }

    /**
     * The derivative r*du/dr.
     */
    public double du(double r2) { throw new MethodNotImplementedException(); }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) { throw new MethodNotImplementedException(); }

    /**
     *  Integral used for corrections to potential truncation.
     */
    public double uInt(double rC) { throw new MethodNotImplementedException(); }

    public static void main(String[] args){
        Space space = Space3D.getInstance();
        final P2HePCJS pnew = new P2HePCJS(space);
        final P2HePCKLJS pold = new P2HePCKLJS(space);
        //double V1 = p.V(5.6,false);
        double tempK = 300; // Kelvin
        double beta = 1/Kelvin.UNIT.toSim(tempK);
        for (double r = 0; r <=5; r+=0.1){
            double unew = pnew.u(r*r);
            double uold = pold.u(r*r);
            double uratio = unew/uold;
            double f = Math.exp(-unew*beta) - 1;
            System.out.println(r+" "+" "+f);
        }
        //double V2 = p.V(0.0,true);
        //double dampTT = p.dampTT(16,11.941151488133176);
        //System.out.println(V2);
    }

}
