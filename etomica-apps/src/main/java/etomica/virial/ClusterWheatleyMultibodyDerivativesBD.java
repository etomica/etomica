/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IMoleculeList;
import etomica.atom.MoleculeArrayList;
import etomica.math.SpecialFunctions;

import java.math.BigDecimal;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams, but adds in non-additive contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibodyDerivativesBD extends ClusterWheatleySoftDerivativesBD implements ClusterAbstractMultivalue {

    protected final MayerFunctionNonAdditive fNonAdditive;
    protected final MayerFunctionNonAdditive[] fMulti;
    protected final int[] moleculeIndices;
    protected final double[] r2;
    protected final MoleculeArrayList molecules;
    protected boolean doMulti;
    protected double rCut2;
    protected ClusterWheatleyMultibodyDerivativesBD clusterMultiBDBD;
    protected final BigDecimal[] fQmulti;
    protected final int nDer;
    protected int precisionLimit;
    protected final boolean doTotal;
    protected final BigDecimal[] pairValues;


   
    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fNonAdditive Mayer function that returns non-additive value for
     *          any number of molecules.
     * @param fMulti array of non-additive Mayer functions.  fMulti[3] is the
     *          3-body Mayer function (exp(-beta*deltaU3)-1), fMulti[4] is the
     *          4-body Mayer function, etc.  fMulti null entries will be
     *          ignored and the array need to not be of size equal to nPoints.
     *          If only 3-body Mayer function is available, then fMulti can be
     *          of length 4 (0,1,2,3)
     */
    public ClusterWheatleyMultibodyDerivativesBD(int nPoints, MayerFunction f, MayerFunctionNonAdditive fNonAdditive, MayerFunctionNonAdditive[] fMulti, int precision, int nDer, boolean doTotal) {
        super(nPoints, f, precision, nDer);
        this.nDer = nDer;
        this.doTotal=doTotal;
        this.fNonAdditive = fNonAdditive;
        this.fMulti = fMulti;
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // pairwise shouldn't try to use BD                
        molecules = new MoleculeArrayList(nPoints);
        rCut2 = Double.POSITIVE_INFINITY;
        fQmulti = new BigDecimal[1<<n];
        fQmulti[0] = fQmulti[1] = fQmulti[2] = BDONE;
        pairValues = new BigDecimal[nDer+1];
    }

    public void setPrecisionLimit(int newLimit) {
        if (newLimit > 300) newLimit = 300;
        precisionLimit = newLimit;
        clusterMultiBDBD = null;
    }


    public ClusterAbstract makeCopy() {
        ClusterWheatleyMultibodyDerivativesBD c = new ClusterWheatleyMultibodyDerivativesBD(n, f, fNonAdditive, fMulti, mc.getPrecision(),nDer,doTotal);
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
        return c;
    }

    public void setTemperature(double newT) {
        super.setTemperature(newT);
        if (clusterMultiBDBD != null) {
            clusterMultiBDBD.setTemperature(newT);
        }
    }

    public void setRCut(double newRCut) {
        rCut2 = newRCut * newRCut;
    }

    public void calcValue(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        double rMax = 0;
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                if (cPairs.getr2(i,j) > rCut2) {
                    value[0] = 0;
                    return;
                }
                if (cPairs.getr2(i,j) > rMax) rMax = cPairs.getr2(i,j);
            }
        }
        
        // do (multi+pair) - pair here so that we avoid recomputing f bonds
        int nf = 1 << n;
        if (!doTotal) {
            doMulti = false;
            super.calcValue(box);            
            System.arraycopy(fB[nf-1], 0, pairValues, 0, nDer+1);
        }
        doMulti = true;
        super.calcValue(box);    
        double bfac = (1.0-n)/SpecialFunctions.factorial(n);
        if (!doTotal) {
            for (int m = 0; m <= nDer; m++) {
                fB[nf - 1][m] = fB[nf - 1][m].subtract(pairValues[m], mc);
                value[m] = bfac * fB[nf - 1][m].doubleValue();
            }
        }
        // we have our own BD cluster for multibody
        // if an individual integrand (pair/total) is small, it won't trigger a BD calculation
        // BD only gets triggered here        
        if (Math.log10(Math.abs(fB[nf-1][0].doubleValue())) < -(mc.getPrecision()-5)) {
            // value is too small for us to compute it precisely
            if (mc.getPrecision() >= precisionLimit) {
                for (int m=0; m<=nDer;m++){           
                    value[m] = 0;
                }
                return;
            }
            if (clusterMultiBDBD == null) {
                int p = mc.getPrecision() + 20;
                if (p > precisionLimit) p = precisionLimit;
                clusterMultiBDBD = new ClusterWheatleyMultibodyDerivativesBD(n, f,fNonAdditive, fMulti, p,nDer,doTotal);
                clusterMultiBDBD.setTemperature(1/beta);
                clusterMultiBDBD.setDoCaching(doCaching);
                clusterMultiBDBD.setPrecisionLimit(precisionLimit);
            }
            System.arraycopy(clusterMultiBDBD.getAllLastValues(box), 0, value, 0, nDer+1);
            return;
        }
    }

    protected void calcFullFQ(BoxCluster box) {
        super.calcFullFQ(box);
        boolean debugme = false; //box.getIndex()==1;
        if(debugme)System.out.println(fQ[3][0]+" "+fQ[5][0]+" "+fQ[6][0]+" "+fQ[7][0]);        
        if (!doMulti) return;
        for (int i=3; i<fMulti.length; i++) {
            if (fMulti[i]!=null) fMulti[i].setBox(box);
        }
        if (fNonAdditive != null) {
            fNonAdditive.setBox(box);
        }
        int nf = 1<<n;
        IMoleculeList boxMolecules = box.getMoleculeList();
        // FQ[i] now contains the exp(-bU2) where U2 is the pair-wise energy for set i.
        // we need to go around and add the non-additive energy for each set.

        for (int i=3; i<nf; i++) {
            fQmulti[i] = BDONE;
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) {
                if (fQ[i][0].doubleValue() == 0){
                    for (int m=1;m<=nDer;m++){
                        fQ[i][m]=BDZERO;
                    }
                    continue;
                }                
                BigDecimal c = BDlog(fQ[i][0]).divide(BDbeta,mc);
                for(int m=1; m<=nDer;m++){
                    fQ[i][m]=fQ[i][m-1].multiply(c, mc);      //calculates derivatives of fQ w.r.t. beta      
                }                
                continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            }
            
            if (fQ[i][0].doubleValue() == 0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=BDZERO;
                } continue;
            }            
            // we want to loop over subsets of i with at least 3 points
            // we could try to be clever and use the moleculeIndices constructed below to loop
            // more efficiently here.
            int iLowBit = (i & -i);//next lowest bit
            for (int isub=iLowBit; isub<i; isub+=iLowBit) {//sum over partitions of i
                while ((isub & ~i) != 0) {
                    // loop until isub is an actual subset of i
                    isub += iLowBit;
                }
                fQ[i][0] = fQ[i][0].multiply(fQmulti[isub],mc);
                if (fQ[i][0].doubleValue()==0){
                    break;
                }
                
            }
            if (fQ[i][0].doubleValue()==0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=BDZERO;
                } continue;
            }            
            int l = 0;
            molecules.clear();
            for (int a=0; a<n; a++) {
                if ((i & (1<<a)) != 0) {
                    moleculeIndices[l] = a;
                    molecules.add(boxMolecules.getMolecule(a));
                    l++;
                }
            }
            if ((fMulti.length <= l || fMulti[l] == null) && fNonAdditive == null) {
                fQmulti[i] = BDONE;
                continue;
            }
            int ll = 0;
            for (int a=0; a<l-1; a++) {
                for (int b=a+1; b<l; b++) {
                    r2[ll] = box.getCPairSet().getr2(moleculeIndices[a],moleculeIndices[b]);
                    ll++;
                }
            }
            if (fMulti.length > l && fMulti[l] != null) {
                fQmulti[i] = new BigDecimal(fMulti[l].f(molecules, l, moleculeIndices, r2, beta),mc).add(BDONE);
            }
            fQ[i][0] = fQ[i][0].multiply(fQmulti[i],mc);
            if (fNonAdditive != null) {
                // we don't want to include this in fQmulti because we would just include it again
                // for larger sets
                fQ[i][0] = fQ[i][0].multiply(new BigDecimal(fNonAdditive.f(molecules, l, moleculeIndices, r2, beta)).add(BDONE), mc);
            }
            if (fQ[i][0].doubleValue() == 0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=BDZERO;
                }
                continue;
            }
            BigDecimal c = BDlog((fQ[i][0])).divide(BDbeta,mc);
            for(int m=1; m<=nDer;m++){
                fQ[i][m]=fQ[i][m-1].multiply(c,mc);      //calculates derivatives of fQ w.r.t. beta      
            }
            if(debugme)System.out.println(fQ[7][0]);
        }
    }
    
    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }
    
    public static class ClusterRetrievePrimes implements ClusterAbstract{
        protected final ClusterWheatleyMultibodyDerivativesBD cluster;
        protected final int n;
        
        public ClusterRetrievePrimes(ClusterWheatleyMultibodyDerivativesBD cluster, int n){
            this.cluster=cluster;
            this.n = n;
        }
        
        public ClusterAbstract makeCopy() {
            
            return null;
        }

        
        public int pointCount() {
            return cluster.n;
        }

        
        public double value(BoxCluster box) {
            return cluster.value[n];
        }

        public void setTemperature(double temperature) {
            
            
        }
        
        
    }
    
}
