/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.veos;

import java.math.BigDecimal;

import etomica.math.SpecialFunctions;

public class Ushcats {

    private BetaSource betaSource;
    int N;
    IBigValue[] M;
    IBigValue[] dMdV;
    private double rhoLast;
    
    public Ushcats(BetaSource bs, int N) {
        betaSource = bs;
        setN(N);
    }
    
    protected IBigValue makeBigValue(double value) {
//        return new BigValue(value);
        return new MyBigDecimal(value);
    }
    
    protected IBigValue makeBigValue(IBigValue value) {
//        return new BigValue(value);
        return new MyBigDecimal(value);
    }
    
    public final void setN(int n) {
        N = n;
        M = new IBigValue[N+1];
        dMdV = new IBigValue[N+1];
        rhoLast = -1.;
    }
    
    private void calcM(double rho) {
        if(rho == rhoLast) return;
        M[0] = makeBigValue(1.0);
        IBigValue rhoN = makeBigValue(rho/N);
        IBigValue term = makeBigValue(0.0);
        IBigValue t = makeBigValue(0.0);
        for(int k=1; k<=N; k++) {
            M[k] = makeBigValue(M[k-1]);
            IBigValue fact = makeBigValue(rhoN);// equals (k-1)!/(k-i)!(rho/N)^k for each iteration
            int imax = Math.min(k, betaSource.maxIndex());
            for(int i=1; i<=imax; i++) {
                term.E(betaSource.kbeta(i));
                term.TE(fact);
                t.E(M[k-i]);
                t.TE(N);
                if(i<k) {
                    t.PEa1Tv1(-(k-i), M[k-1-i]);
                    fact.TE(k-i);
                    fact.TE(rhoN);
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
        dMdV[0] = makeBigValue(0.0);
        IBigValue rhoN = makeBigValue(rho/N);
        IBigValue term = makeBigValue(0.0);
        IBigValue t = makeBigValue(0.0);
        IBigValue dt = makeBigValue(0.0);
        for(int k=1; k<=N; k++) {
            dMdV[k] = makeBigValue(dMdV[k-1]);
            IBigValue fact = makeBigValue(rhoN);// equals (k-1)!/(k-i)!(rho/N)^k for each iteration
            int imax = Math.min(k, betaSource.maxIndex());
            for(int i=1; i<=imax; i++) {
                term.E(betaSource.kbeta(i));
                term.TE(fact);
                t.E(M[k-i]);
                dt.E(dMdV[k-i]);
                t.TE(N);
                dt.TE(N);
                if(i<k) {
                    t.PEa1Tv1(-(k-i), M[k-1-i]);
                    dt.PEa1Tv1(-(k-i), dMdV[k-1-i]);
                    fact.TE(k-i);
                    fact.TE(rhoN);
                }
                t.TE(-i);
                t.TE(rhoN);
                t.PE(dt);
                term.TE(t);                
                dMdV[k].PE(term);
            }
        }
    }


    /**
     * Computes q, given as Q * N!/V^N, where Q is the configurational integral.
     */
    public IBigValue qCalc(double rho, int n) {
        calcM(rho);
        IBigValue rhoN = makeBigValue(rho/N);
        IBigValue qN = makeBigValue(1.0);
        IBigValue ksum = makeBigValue(0.0);
        IBigValue term = makeBigValue(0.0);
        for(int i=1; i<=n; i++) {
            ksum.E(0.0);
            IBigValue fact = makeBigValue(rhoN);// equals (i-1)!/(i-k)!(rho/N)^k for each iteration
            int kmax = Math.min(i, betaSource.maxIndex());
            for(int k=1; k<=kmax; k++) {
                term.E(betaSource.kbeta(k));
                term.TE(fact);
                term.TE(M[i-k]);
                ksum.PE(term);
                if(k<i) {
                    fact.TE(i-k);
                    fact.TE(rhoN);
                }
            }
            qN.PEa1Tv1((N-i),ksum);
        }
        return qN;
    }

    public IBigValue dqdVCalc(double rho, int n) {
        calcdMdV(rho);
        IBigValue rhoN = makeBigValue(rho/N);
        IBigValue dqN = makeBigValue(0.0);
        IBigValue ksum = makeBigValue(0.0);
        IBigValue term = makeBigValue(0.0);
        IBigValue t = makeBigValue(0.0);
        for(int i=1; i<=n; i++) {
            ksum.E(0.0);
            IBigValue fact = makeBigValue(rhoN);// equals (i-1)!/(i-k)!(rho/N)^k for each iteration
            int kmax = Math.min(i, betaSource.maxIndex());
            for(int k=1; k<=kmax; k++) {
                term.E(betaSource.kbeta(k));
                term.TE(fact);
                t.E(M[i-k]);
                t.TE(-k);
                t.TE(rhoN);
                t.PE(dMdV[i-k]);
                term.TE(t);
                ksum.PE(term);
                if(k<i) {
                    fact.TE(i-k);
                    fact.TE(rhoN);
                }
            }
            dqN.PEa1Tv1((N-i),ksum);
        }
        return dqN;
    }
    
    public double pCalc(double rho, double T) {
        IBigValue dqN = dqdVCalc(rho, N);
        //IBigValue qN = qCalc(rho,N);
        dqN.DE(qCalc(rho,N));
        return T*(rho + dqN.value());
    }

    private static void printList(String label, int n, double[] x, double[] y) {
        System.out.print(label + "["+Integer.toString(n)+"] = {");
        for(int i=0; i<x.length; i++) {
            if(Double.isInfinite(y[i]) || Double.isNaN(y[i])) continue;
            System.out.print("{"+x[i]+",  "+new BigDecimal(y[i]).toPlainString()+"}");
            if(i != x.length-1) System.out.println(",");
            else System.out.println("};");
        }
        System.out.println();
    }
        
    public static void main(String[] args) {
        System.out.println("Starting...");

        int N = 100;
        int nRho = 31;
        int maxIndex = 4;//Integer.MAX_VALUE;//largest nonzero-virial index
        boolean doQ = false;
        boolean doP = !doQ;
        boolean doSeries = (N < 0);
        boolean dovdW = false;
        int outputIndex = maxIndex;//used to writing identifying parameter value on output

        BetaSource betaSource;
        double T;
        if(dovdW) {
            T = .8 * 8.0 / 27.0;
            betaSource = new VanderWaals(T, maxIndex);
        } else {
            T = 1.0;
            betaSource = new HardRods(maxIndex);
        }
        Ushcats ush = new Ushcats(betaSource,N);
        double[] rho = new double[nRho];
        for(int i=0; i<rho.length; i++) {
            rho[i] = (i)*0.9/(rho.length-1);
        }

        double tSeries = 0.0;
        long t1, t2;
        double[] lnqSer = new double[rho.length];
        if(doQ) {
            if(doSeries) {
                t1 = System.currentTimeMillis();
                QSeries qSeries = new QSeries(betaSource);
                double[] qSer = qSeries.qCalc(N, rho);
                for(int j=1; j<rho.length; j++) {
                    lnqSer[j] = Math.log(qSer[j])/N;
                }
                t2 = System.currentTimeMillis();
                tSeries = (t2-t1)/1000.;
                printList("lnqSer",outputIndex,rho,lnqSer);
            }

            double[] lnqUsh = new double[rho.length];
            t1 = System.currentTimeMillis();
            for(int j=1; j<rho.length; j++) {
                if(rho[j]==0) continue;
                IBigValue qUsh = ush.qCalc(rho[j],N);
                if(!qUsh.isPositive()) {
                    System.out.println("Negative Q at rho = "+rho[j]);
                }
                lnqUsh[j] = qUsh.lnValue()/N;
            }
            t2 = System.currentTimeMillis();
            double tUsh = (t2-t1)/1000.;
            printList("lnqUsh",outputIndex,rho,lnqUsh);

            if(doSeries) {
                double[] difference = new double[rho.length];
                for(int j=0; j<rho.length; j++) {
                    difference[j] = Math.log10(Math.abs(lnqSer[j] - lnqUsh[j]));
                }
                printList("lnqDif",outputIndex,rho,difference);
            }

            System.out.println("Timings (series, Ushcats): "+tSeries+", "+tUsh);

        }
        
        if(doP) {
            double[] pUsh = new double[rho.length];
            t1 = System.currentTimeMillis();
            for(int j=0; j<rho.length; j++) {
                if(rho[j]==0) pUsh[j] = 0.0;
                else pUsh[j] = ush.pCalc(rho[j],T);
            }
            t2 = System.currentTimeMillis();
            double tUsh = (t2-t1)/1000.;
            printList("pUsh",outputIndex,rho,pUsh);

            System.out.println("Timing (Ushcats): "+tUsh);

        }

    }
    
}
