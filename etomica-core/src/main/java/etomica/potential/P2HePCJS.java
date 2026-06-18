/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.potential;

import etomica.exception.MethodNotImplementedException;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

/**
 * pair potential for Helium Przybytek et al. (2017) Phys. Rev. Lett. 119, 123401.
 *
 * @author Navneeth Gokul
 */
public class P2HePCJS implements IPotential2 {

    public static IPotential2 makeTruncated(double sigma, TruncationFactory tf) {
        return tf.make(new P2HePCJS(sigma));
    }

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
    protected final double[][] aS_BO, aS_AD, aS_REL, aS_QED;
    protected final double sigma;

    public P2HePCJS() {
        this(0);
    }

    public P2HePCJS(double sigma) {
        super();
        this.sigma = sigma;

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

        aS_BO = new double[][]{{0.80740, 1.2195e-6},
                {1.2603, 2.6076e-5},
                {0.13365, -1.6989e-7},
                {0.057123, -5.9399e-8}};
        aS_AD = new double[][]{{0.30940, 7.5317e-12},
                {1.8485, 5.8720e-6},
                {0.27382, 4.8110e-8},
                {0.047007, 1.0593e-10}};
        aS_REL = new double[][]{{0.20535, 3.0628e-12},
                {1.6829, 2.3492e-7},
                {0.42205, -1.2725e-6},
                {0.39672, 1.2107e-6},
                {0.064498, 3.8763e-10}};
        aS_QED = new double[][]{{0.80740, 3.3214e-7},
                {0.50144, -1.0022e-7},
                {0.16251, -6.5718e-8},
                {0.062390, -2.1856e-8}};
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

    private double intPow(double x, int a) {
        if (a < 0) return 1.0 / intPow(x, -a);
        if (a < 4) {
            switch (a) {
                case 0:
                    return 1;
                case 1:
                    return x;
                case 2:
                    return x * x;
                case 3:
                    return x * x * x;
            }
        }
        double x2 = x * x;
        switch (a) {
            case 4:
                return x2 * x2;
            case 5:
                return x2 * x2 * x;
            case 6:
                return x2 * x2 * x2;
        }
        throw new RuntimeException("oops, don't know how to handle a=" + a);
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
        if (I0 != 0) formShort = formShort * intPow(r, I0);

        return formShort;
    }

    private double dformShort(double r, int M, int I0, int I1, double[] a, double[][] P) {
        //short-range formula
        //sum_{k=1}^{M} exp(-a(k) R) sum_{i=I0}^{I1} P(i,k) R^i
        if (I1 - 1 <= I0) return 0;
        double formShort = 0.0, dformShort = 0.0;
        for (int k = 0; k < M; k++) {
            double term = P[I1 - I0][k]; // subtracted I1 from i for appropriate index
            double dterm = (I1 - I0) * P[I1 - I0][k] / r;
            //if(I0==-1)System.out.println("k = " + k +" I1 = "+ I1 +" P[I1][k] = " + term);
            for (int i = I1 - 1; i >= I0; i--) {
                term = term * r + P[i - I0][k]; // subtracted I0 from i for appropriate index
                if (i >= I0 + 1) dterm += dterm * r + (i - I0) * P[i - I0][k];
                //if(I0==-1)System.out.println("k = " + k + " i = "+ i +" P[i][k] = " + P[i-I0][k] + " term = " + term);
            }
            formShort += Math.exp(-a[k] * r) * term;
            dformShort += Math.exp(-a[k] * r) * dterm - a[k] * Math.exp(-a[k] * r) * term;
            //if(I0==-1)System.out.println("k = " + k +" fs = " + formShort);
        }
        if (I0 != 0) dformShort = dformShort * intPow(r, I0) + formShort * I0 * intPow(r, I0 - 1);

        return dformShort;
    }

    private double d2formShort(double r, int M, int I0, int I1, double[] a, double[][] P) {
        //short-range formula
        //sum_{k=1}^{M} exp(-a(k) R) sum_{i=I0}^{I1} P(i,k) R^i
        if (I1 - 1 <= I0 - 1) return 0;
        double formShort = 0.0, dformShort = 0.0, d2formShort = 0.0;
        for (int k = 0; k < M; k++) {
            double term = P[I1 - I0][k]; // subtracted I1 from i for appropriate index
            double dterm = (I1 - I0) * P[I1 - I0][k] / r;
            double d2term = (I1 - I0) * (I1 - I0 - 1) * P[I1 - I0][k] / (r * r);

            //if(I0==-1)System.out.println("k = " + k +" I1 = "+ I1 +" P[I1][k] = " + term);
            for (int i = I1 - 1; i >= I0; i--) {
                term = term * r + P[i - I0][k]; // subtracted I0 from i for appropriate index
                if (i > I0) dterm += dterm * r + (i - I0) * P[i - I0][k];
                if (i > I0 + 1) d2term += d2term * r + (i - I0) * (i - I0 - 1) * P[i - I0][k];
                //if(I0==-1)System.out.println("k = " + k + " i = "+ i +" P[i][k] = " + P[i-I0][k] + " term = " + term);
            }
            formShort += Math.exp(-a[k] * r) * term;
            dformShort += Math.exp(-a[k] * r) * dterm - a[k] * Math.exp(-a[k] * r) * term;
            d2formShort += Math.exp(-a[k] * r) * d2term - 2 * a[k] * Math.exp(-a[k] * r) * dterm + a[k] * a[k] * Math.exp(-a[k] * r) * term;
            //if(I0==-1)System.out.println("k = " + k +" fs = " + formShort);
        }
        if (I0 == 0) return d2formShort;
        if (I0 == 1) {
            d2formShort = d2formShort * intPow(r, I0) + dformShort * I0 * intPow(r, I0 - 1);
        } else {
            d2formShort = d2formShort * intPow(r, I0) + dformShort * I0 * intPow(r, I0 - 1) + formShort * I0 * (I0 - 1) * intPow(r, I0 - 2);
        }

        return d2formShort;
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
        double formLong = 0;

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

        formLong = formLong * intPow(iR, N0);
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

    public double sigma_BO(double R) {
        return formSigma(R, aS_BO);
    }

    public double sigma_AD(double R) {
        return formSigma(R, aS_AD);
    }

    public double sigma_REL(double R) {
        return formSigma(R, aS_REL);
    }

    public double sigma_QED(double R) {
        return formSigma(R, aS_QED);
    }

    public double sigma(double R) {
        // sum the squared errors, because that makes sense.
        double sum = 0, x;
        x = sigma_BO(R);
        sum += x * x;
        x = sigma_AD(R);
        sum += x * x;
        x = sigma_REL(R);
        sum += x * x;
        x = sigma_QED(R);
        sum += x * x;
        return Math.sqrt(sum);
    }

    public double formSigma(double R, double[][] myAS) {
        double R2 = R * R;
        double rv = Math.exp(-myAS[0][0] * R) * myAS[0][1];
        for (int i = 1; i < myAS.length; i++) {
            rv += Math.exp(-myAS[i][0] * R2) * myAS[i][1];
        }
        return rv;
    }

    /**
     * @return energy u.
     * @param r2 square of the distance in Angstrom^2.
     */
    public double u(double r2) {
        if (r2 < 1e-6) return Double.MAX_VALUE;

        double r = Math.sqrt(r2); //distance in Angstrom
        double rBohr = BohrRadius.UNIT.fromSim(r); //distance in Bohr

        double uHartree =  V(rBohr,true); //energy in Hartree. Always account for retardation
        double uSigma = 0;
        if (sigma != 0) {
            uHartree += sigma * sigma(rBohr);
        }

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

    public static void main(String[] args){
        final P2HePCJS pnew = new P2HePCJS();
        final P2HePCJS pnewp = new P2HePCJS(1);
        final P2HePCKLJS pold = new P2HePCKLJS();
        final P2HePCKLJS poldp = new P2HePCKLJS(1);
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

        System.out.println("\nerror");
        for (double r = 0.5; r < 10.0001; r += 0.1) {
            double unew = pnew.u(r * r);
            double unewp = pnewp.u(r * r);
            double uold = pold.u(r * r);
            double uoldp = poldp.u(r * r);
            if (uold == Double.POSITIVE_INFINITY || unew == Double.POSITIVE_INFINITY) continue;
            System.out.println(r + " " + unew + " " + (unewp - unew) + " " + (uoldp - uold) + " " + (uold - unew));
        }
        //double V2 = p.V(0.0,true);
        //double dampTT = p.dampTT(16,11.941151488133176);
        //System.out.println(V2);
    }

}
