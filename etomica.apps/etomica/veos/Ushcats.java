package etomica.veos;

import etomica.math.SpecialFunctions;

public class Ushcats {

    private BetaSource betaSource;
    int N;
    BigValue[] M;
    BigValue[] dMdV;
    private double rhoLast;
    
    public Ushcats(BetaSource bs, int N) {
        betaSource = bs;
        setN(N);
    }
    
    public final void setN(int n) {
        N = n;
        M = new BigValue[N+1];
        dMdV = new BigValue[N+1];
        rhoLast = -1.;
    }
    
    private void calcM(double rho) {
        if(rho == rhoLast) return;
        M[0] = new BigValue(1.0);
        double lnrhoN = Math.log(rho/N);
        BigValue term = new BigValue();
        BigValue t = new BigValue();
        for(int k=1; k<=N; k++) {
            M[k] = new BigValue(M[k-1]);
            for(int i=1; i<=k; i++) {
                term.E(i*betaSource.beta(i));
                term.TEln(SpecialFunctions.lnGamma(k)-SpecialFunctions.lnGamma(k-i+1) + i*lnrhoN);
                t.E(M[k-i]);
                t.TE(N);
                if(i<k) {
                    t.PEa1Tv1(-(k-i), M[k-1-i]);
                }
                term.TE(t);                
                M[k].PE(term);
            }
        }
        rhoLast = rho;
    }
    
    private void calcdMdV(double rho) {
        if(rho == rhoLast) return;
        calcM(rho);
        dMdV[0] = new BigValue(0.0);
        double V = N/rho;
        double lnrhoN = -Math.log(V);
        BigValue term = new BigValue();
        BigValue t = new BigValue();
        BigValue dt = new BigValue();
        for(int k=1; k<=N; k++) {
            dMdV[k] = new BigValue(dMdV[k-1]);
            for(int i=1; i<=k; i++) {
                term.E(i*betaSource.beta(i));
                term.TEln(SpecialFunctions.lnGamma(k)-SpecialFunctions.lnGamma(k-i+1) + i*lnrhoN);
                t.E(M[k-i]);
                dt.E(dMdV[k-i]);
                t.TE(N);
                dt.TE(N);
                if(i<k) {
                    t.PEa1Tv1(-(k-i), M[k-1-i]);
                    dt.PEa1Tv1(-(k-i), dMdV[k-1-i]);
                }
                t.TE(-i/V);
                t.PE(dt);
                term.TE(t);                
                dMdV[k].PE(term);
            }
        }
    }


    /**
     * Computes q, given as Q * N!/V^N, where Q is the configurational integral.
     */
    public BigValue qCalc(double rho, int n) {
        calcM(rho);
        double lnrhoN = Math.log(rho/N);
        BigValue qN = new BigValue(1.0);
        BigValue ksum = new BigValue();
        BigValue term = new BigValue();
        for(int i=1; i<=n; i++) {
            ksum.E(0.0);
            double lngammai = SpecialFunctions.lnGamma(i);
            for(int k=1; k<=i; k++) {
                term.E(k*betaSource.beta(k));
                term.TEln(k*lnrhoN+lngammai-SpecialFunctions.lnGamma(i-k+1));
                term.TE(M[i-k]);
                ksum.PE(term);
            }
            qN.PEa1Tv1((N-i),ksum);
        }
        return qN;
    }

    public BigValue dqdVCalc(double rho, int n) {
        calcdMdV(rho);
        double V = N/rho;
        double lnrhoN = -Math.log(V);
        BigValue dqN = new BigValue(0.0);
        BigValue ksum = new BigValue();
        BigValue term = new BigValue();
        BigValue t = new BigValue();
        for(int i=1; i<=n; i++) {
            ksum.E(0.0);
            double lngammai = SpecialFunctions.lnGamma(i);
            for(int k=1; k<=i; k++) {
                term.E(k*betaSource.beta(k));
                term.TEln(k*lnrhoN+lngammai-SpecialFunctions.lnGamma(i-k+1));
                t.E(M[i-k]);
                t.TE(-k/V);
                t.PE(dMdV[i-k]);
                term.TE(t);
                ksum.PE(term);
            }
            dqN.PEa1Tv1((N-i),ksum);
        }
        return dqN;
    }
    
    public double pCalc(double rho, double T) {
        BigValue dqN = dqdVCalc(rho, N);
        dqN.DE(qCalc(rho,N));
        return T*(rho + dqN.value());
    }

//    public double calcQold(double rho) {
//        //calcM(rho);
//        double qN = 1.0;
//        for(int i=1; i<=N; i++) {
//            double ksum = 0.0;
//            double lngammai = SpecialFunctions.lnGamma(i);
//            for(int k=1; k<=i; k++) {
//                ksum += k*betaSource.beta(k)*M[i-k] * 
//                        Math.exp(lngammai-SpecialFunctions.lnGamma(i-k+1)+k*Math.log(rho/N));
//            }
//            qN += (N-i)*ksum;
//        }
//        return qN;
//    }
        
    public static void main(String[] args) {
        System.out.println("Starting...");
        double T = 10;//.8 * 8.0 / 27.0;
        int N = 1000;
        BetaSource betaSource = new VanderWaals(T,20);
        Ushcats ush = new Ushcats(betaSource,N);
        double[] rho = new double[20];
        for(int i=0; i<rho.length; i++) {
            rho[i] = (i+1)*0.3/rho.length;
        }
//        for(int j=0; j<rho.length; j++) {
//            double lnq = ush.calcLnQ(rho[j],N).lnValue();
//            double q = Math.exp(lnq);
//            double qOld = ush.calcQold(rho[j]);
//            if(Math.abs(q-qOld)/qOld > 1e-10) System.out.println("uh oh:"+q+" "+qOld);
//            System.out.println("{"+rho[j]+",  "+q+",  "+qOld+", "+(-lnq/N+Math.log(rho[j]))+"},");
//        }
        QSeries qSeries = new QSeries(betaSource);
        if(false) {
            double[] qVals = qSeries.qCalc(N, rho); 
            for(int j=0; j<rho.length; j++) {
                double lnq = ush.qCalc(rho[j],N).lnValue();
                double q = Math.exp(lnq);
    //            double qOld = ush.calcQold(rho[j]);
    //            if(Math.abs(q-qOld)/qOld > 1e-10) System.out.println("uh oh:"+q+" "+qOld);
                System.out.println("{"+rho[j]+",  "+q+",  "+qVals[j]+"},");
            }
        }
        if(true) {
//            double[] pVals = qSeries.pCalc(N, rho, T); 
            for(int j=0; j<rho.length; j++) {
                double p = ush.pCalc(rho[j],T);
//                System.out.println("{"+rho[j]+",  "+p+",  "+pVals[j]+", "+(p-pVals[j])/p+"},");
                System.out.println("{"+rho[j]+",  "+p+"},");
            }
        }

    }
    
}
