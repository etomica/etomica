/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math;

/**
 * This class calculates the incomplete gamma function and complementary
 * incomplete gamma function.  This is a Java translation of
 *   Algorithm 542: Incomplete Gamma Functions [S14] by W. Gautschi from
 *     ACM Transactions on Mathematical Software (TOMS)
 *     Volume 5 ,  Issue 4  (December 1979)
 *     Pages: 482 - 489
 *  see http://www.netlib.org/toms/542
 *
 * The original algorithm was "taylored to the accuracy requirements of the CDC
 * 6500 computer."  Little to no attempt has been made to re-optimize for Java,
 * other than to use Double.{MAX,MIN}_VALUE in place of defined constants.
 *
 * @author Andrew Schultz
 */
public class Gamma {

    protected static final float PREC = 28.8989F;
    protected static final int TOPEXP = (int)(Math.log(Double.MAX_VALUE)/Math.log(10));
    protected static final int BOTEXP = (int)(Math.log(Double.MIN_VALUE)/Math.log(10));
    protected static final double AL10 = Math.log(10);
    protected static final double ALPREC = AL10*PREC;
    protected static final double TOP = AL10*TOPEXP;
    protected static final double BOT = AL10*BOTEXP;

    protected static final double[] C = new double[]{
        .57721566490153286060651209008,-.65587807152025388107701951515,
        -4.200263503409523552900393488E-2,.1665386113822914895017007951,
        -4.21977345555443367482083013E-2,-9.6219715278769735621149217E-3,
         7.2189432466630995423950103E-3,-1.165167591859065112113971E-3,
        -2.15241674114950972815730E-4,1.2805028238811618615320E-4,
        -2.013485478078823865569E-5,-1.2504934821426706573E-6,
         1.1330272319816958824E-6,-2.056338416977607103E-7,
         6.1160951044814158E-9,5.0020076444692229E-9,-1.181274570487020E-9,
         1.04342671169110E-10,7.782263439905E-12,-3.696805618642E-12,
         5.1003702875E-13,-2.058326054E-14,-5.34812254E-15,1.2267786E-15,
        -1.181259E-16,1.187E-18,1.412E-18,-2.30E-19,1.7E-20
    };
    
    protected double G, GSTAR;
    protected int IFLGST;

    /**
     * Return the normalized complementary incomplete gamma function,
     * (integral from t=x to t=inf of exp(-t)*t^(a-1)) / gamma(a)  if a>0
     * or
     * (integral from t=x to t=inf of exp(-t)*t^(a-1)) * exp(x)*x^(-a) if a<=0
     */
    public double calcGamma(double a, double x) {
        calcBoth(a, x);
        return G;
    }
    
    /**
     * Return the normalized incomplete gamma function,
     * (integral from t=0 to t=x of exp(-t)*t^(a-1)) / x^a  if a>0
     * 
     * If a=0 and x=0, gammaStar is indeterminate and an exception is thrown.
     * gammaStar approaches 1 as x approaches 0 while a=0.
     */
    public double calcGammaStar(double a, double x) {
        calcBoth(a, x);
        if (IFLGST == 3) {
            throw new RuntimeException("GSTAR is indeterminate");
        }
        return GSTAR;
    }

    /**
     * Return the normalized complementary incomplete gamma function from the
     * last calculation.  This is useful if you called calcGammaStar, but now
     * want gamma (gamma and gammaStar are both calculated).
     */
    public double getLastGamma() {
        return G;
    }
    
    /**
     * Return the normalized incomplete gamma function from the last
     * calculation.  This is useful if you called calcGamma, but now want
     * gammaStar (gamma and gammaStar are both calculated).
     */
    public double getLastGammaStar() {
        if (IFLGST == 3) {
            throw new RuntimeException("GSTAR is indeterminate");
        }
        return GSTAR;
    }
    
    protected void calcBoth(double A, double X) {
        double ACC = 14;
        
        if (X<0) {
            throw new IllegalArgumentException("X must be non-negative");
        }

        // initialization
        G = 0;
        GSTAR = 0;
        double AA = Math.abs(A);
        
        IFLGST = 0;
        if (X == 0) {
            // EVALUATION OF GSTAR(A,0.) AND G(A,0.)
            if (A<0) {
                GSTAR = Double.POSITIVE_INFINITY;
                G = 1/AA;
            }
            else if (A==0) {
                IFLGST = 3;
                GSTAR = 1;
                G = Double.POSITIVE_INFINITY;
            }
            else {
                G = 1;
            }
            return;
        }
        

        int I = 0;
        double ALX = 0;
        if (X>0) ALX = Math.log(X);
        double ALPHA = X + .25;
        if (X<.25 && X>0.) ALPHA = Math.log(.5)/ALX;
        double EPS = .5*Math.pow(10,-ACC);
        double EPS1 = EPS/100;
        double SGA = 1;
        if (A<0) SGA = -SGA;
        double AE = A;
        double AP1 = A + 1;
        double AEP1 = AP1;
        int MA = (int)(.5-A);
        double FMA = MA;
        double AEPS = A + FMA;
        double SGAE = 1;
        if (AEPS<0) SGAE = -SGAE;
        double AAEPS = Math.abs(AEPS);
        double ALGP1 = 0;
        
        // EVALUATION OF THE LOGARITHM OF THE ABSOLUTE VALUE OF
        // GAMMA(A+1.) AND DETERMINATION OF THE SIGN OF GAMMA(A+1.)
        
        double SGGA = 1;
        if (MA>0) {
            if (AEPS!=0) {
                SGGA = SGAE;
                if (MA==2*(MA/2)) SGGA = -SGGA;
                ALGP1 = DLGA(AEPS+1) - Math.log(AAEPS);
                if (MA!=1) {
                    ALGP1 = ALGP1 + DLGA(1-AEPS) - DLGA(FMA-AEPS);
                }
            }
        }
        else {
            ALGP1 = DLGA(AP1);
        }
        double ALGEP1 = ALGP1;
        boolean skip_taylor = false;
        if (A<=ALPHA) {
            if (X>1.5) {
                // EVALUATION OF G(A,X) FOR X>1.5 AND A<=ALPHA(X) BY
                // MEANS OF THE LEGENDRE CONTINUED FRACTION
                
                GSTAR = 1;
                double XPA = X + 1 - A;
                double XMA = X - 1. - A;
                double P = 0;
                double Q = XPA*XMA;
                double R = 4.*XPA;
                double S = -A + 1;
                double TERM = 1;
                double SUM = 1;
                double RHO = 0;
                double K = 1;
                do {
                  K = K + 1;
                  if (K>600) {
                    throw new RuntimeException("Convergence failed");
                  }
                  P = P + S;
                  Q = Q + R;
                  R = R + 8;
                  S = S + 2;
                  double T = P*(1+RHO);
                  RHO = T/(Q-T);
                  TERM = RHO*TERM;
                  SUM = SUM + TERM;
                } while (Math.abs(TERM)>EPS*SUM);
                if (A<0) {
                  G = SUM/XPA;
                  double ALG = Math.log(G);
                  // EVALUATION OF GSTAR(A,X) IN TERMS OF G(A,X)
                  GSTAR = 1;
                  if (MA>=0 && AEPS==0) return;
                  double SGT = SGA*SGGA;
                  double T = Math.log(AA) - X + A*ALX + ALG - ALGP1;
                  if (T<-ALPREC) return;
                  if (T>=TOP) {
                    GSTAR = -SGT*Double.POSITIVE_INFINITY;
                    return;
                  }
                  GSTAR = 1 - SGT*Math.exp(T);
                  return;
                }
                if (A==0) {
                  G = SUM/XPA;
                  return;
                }
                double ALG = A*ALX - X + Math.log(A*SUM/XPA) - ALGP1;
                if (ALG<=BOT) {
                  return;
                }
                G = Math.exp(ALG);
                GSTAR = 1 - G;
                return;
            }
            boolean do_direct = true;
            double U = 0;
            if (A<-.5) {
                I = 1;
                AE = AEPS;
                AEP1 = AEPS + 1;
                if (X<0.25 && AE>ALPHA) do_direct = false;
            }
            if (do_direct) {
                //DIRECT EVALUATION OF G(A,X) AND GSTAR(A,X) FOR X<=1.5
                //AND -.5<=A<=ALPHA(X)
                GSTAR = 1;
                if (A>=.5) {
                    U = Math.exp(DLGA(A)) - (Math.pow(X,A))/A;
                }
                else {
                    double SUM = 0;
                    for (int K = C.length-1; K>-1; K--) {
                        SUM = AE*SUM + C[K];
                    }
                    double GA = -SUM/(1+AE*SUM);
                    double Y = AE*ALX;
                    if (Math.abs(Y)<1) {
                        SUM = 1;
                        double TERM = 1;
                        int K = 1;
                        do {
                            K = K + 1;
                            if (K>600) {
                                throw new RuntimeException("Convergence failed");
                            }
                            TERM = Y*TERM/K;
                            SUM = SUM + TERM;
                        } while (Math.abs(TERM)>EPS1*SUM);
                        U = GA - SUM*ALX;
                    }
                    else {
                        U = GA - (Math.exp(Y)-1)/AE;
                    }
                }
                double P = AE*X;
                double Q = AEP1;
                double R = AE + 3;
                double TERM = 1;
                double SUM = 1;
                int K = 1;
                do {
                    K = K + 1;
                    if (K>600) {
                        throw new RuntimeException("Convergence failed");
                    }
                    P = P + X;
                    Q = Q + R;
                    R = R + 2.;
                    TERM = -P*TERM/Q;
                    SUM = SUM + TERM;
                } while (Math.abs(TERM)>EPS1*SUM);
                double V = Math.pow(X,AEP1)*SUM/AEP1;
                G = U + V;
                if (I!=1) {
                    if (A<0) {
                        double T = Math.exp(X)*Math.pow(X,-A);
                        G = T*G;
                        GSTAR = 1 - A*G*Math.exp(-ALGP1)/T;
                    }
                    else if (A==0) {
                        G = Math.exp(X)*G;
                    }
                    else {
                        G = A*G*Math.exp(-ALGP1);
                        GSTAR = 1 - G;
                    }
                    return;
                }
                skip_taylor = true;
            }
            I = 2;
            ALGEP1 = DLGA(AEP1);
        }

        // EVALUATION OF GSTAR(A,X) FOR A>ALPHA(X) BY TAYLOR EXPANSION

        if (!skip_taylor) {
            G = 1;
            double TERM = 1;
            double SUM = 1;
            int K = 0;
            do {
              K = K + 1;
              if (K>600) {
                  throw new RuntimeException("Convergence failed");
              }
              TERM = X*TERM/(AE+K);
              SUM = SUM + TERM;
            } while (Math.abs(TERM)>EPS*SUM);
            double ALGS = AE*ALX - X + Math.log(SUM) - ALGEP1;
            if (ALGS<=BOT) {
              return;
            }
            GSTAR = Math.exp(ALGS);
            G = 1 - GSTAR;
            if (I!=2) return;
            G = G*Math.exp(ALGEP1)/AE;
        }
        
        G = G*Math.exp(X)*Math.pow(X,-AE);
        for (int K=1; K<MA+1; K++) {
            G = (1-X*G)/(K-AE);
        }
        double ALG = Math.log(G);
        
        GSTAR = 1;
        if (MA>=0 && AEPS==0) return;
        double SGT = SGA*SGGA;
        double T = Math.log(AA) - X + A*ALX + ALG - ALGP1;
        if (T<-ALPREC) return;
        if (T>=TOP) {
            GSTAR = -SGT*Double.POSITIVE_INFINITY;
        }
        GSTAR = 1 - SGT*Math.exp(T);
    }

    protected final static double[] DBNUM = new double[]{-3617,1,-691,1,-1,1,-1,1};
    protected final static double[] DBDEN = new double[]{12240,156,36036,1188,1680,1260,360,12};
    protected final static double DC = 0.5*Math.log(8*Math.atan(1));
    protected static double DLGA(double DX) {
        double DP = 1;
        double DY = DX;
        float Y = (float)DY;
        while (Y<=35.027) {
            DP = DY*DP;
            DY = DY + 1;
            Y = (float)DY;
        }
        double DT = 1/(DY*DY); 
        double DS = 4.3867/2.44188;
        for (int i=0; i<DBNUM.length; i++) {
            DS = DT*DS + DBNUM[i]/DBDEN[i];
        }
        return (DY-.5)*Math.log(DY) - DY + DC + DS/DY - Math.log(DP);
    }
    
    public static void main(String[] args) {
        // perform various calculations to test algorithm.

        Gamma gamma = new Gamma();
        
        System.out.println("erf");
        for (int I=1; I<32; I++) {
            double DX = (I-1)*0.05;
            double DXSQ = DX*DX;
            double DERF = gamma.calcGammaStar(0.5, DXSQ);
            System.out.println(DX+" "+DERF);
        }
        
        System.out.println();
        
        System.out.println("erfc");
        for (int I=1; I<51; I++) {
            double DXMSQ = I*5.E-3;
            double DXSQ = 1./DXMSQ;
            double DG = gamma.calcGamma(.5, DXSQ);
            double DERFC = Math.sqrt(DXSQ)*Math.exp(DXSQ)*DG;
            System.out.println(DXMSQ+" "+DERFC);
        }

        System.out.println();

        double[] DX1 = new double[]{0, 0.01, 0.05, 0.2, 0.5, 1.5, 5.1, 10, 14.7, 19.8};
        double[] DX2 = new double[]{0.01, 0.37, 1.44, 3.02, 6.57, 20};
        
        System.out.println("Ei");
        for (int I = 0; I<DX1.length; I++) {
            double DX = DX1[I];
            for (int J=1; J<22; J++) {
                int N = J-1;
                double DA = -N+1;
                double DG = gamma.calcGamma(DA, DX);
                double DESUBN;
                if (DX<=0) {
                    DESUBN = 0;
                    if (N>1) {
                        DESUBN = DG;
                    }
                }
                else {
                    if (N == 0) {
                        DESUBN = DG/DX;
                    }
                    else {
                        DESUBN = Math.exp(-DX)*DG;
                    }
                }
                System.out.println(N+" "+DA+" "+DX+" "+DESUBN);
            }
            System.out.println();
        }
        
        System.out.println("Ei, part 2");
        for (int I=0; I<DX2.length; I++) {
            double DX = DX2[I];
            for (int J=1; J<12; J++) {
                double DANU = (J-1)*0.1;
                double DA = -DANU + 1;
                double DG = gamma.calcGamma(DA, DX);
                double DESBNU;
                if (J==11) {
                    DESBNU = Math.exp(-DX)*DG;
                }
                else {
                    DESBNU = Math.pow(DX, -DA)*Math.exp(DLGA(DA+1))*DG/DA;
                }
                System.out.println(DANU+" "+DA+" "+DX+" "+DESBNU);
            }
            System.out.println();
        }
        
        System.out.println("chi^2");
        double[] CCHSQ = new double[]{0.1, 1, 2, 4, 6, 8, 15, 20, 60};
        double[] NUMAX = new double[]{6, 12, 16, 22, 27, 30, 30, 30, 30};
        for (int I=0; I<CCHSQ.length; I++) {
            double CHSQ = CCHSQ[I];
            double NUI = NUMAX[I];
            for (int NU=1; NU<NUI+1; NU++) {
                double A = 0.5*NU;
                double X = 0.5*CHSQ;
                double Q = gamma.calcGamma(A, X);
                double P = gamma.calcGammaStar(A, X);
                System.out.println(NU+" "+CHSQ+" "+X+" "+P+" "+Q);
            }
            System.out.println();
        }
    }
}
