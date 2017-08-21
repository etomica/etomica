/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;

/**
 * This class uses Wheatley's recursion approach to calculating all biconnected
 * diagrams, but adds in non-additive contributions.
 * 
 * @author Andrew Schultz
 */
public class ClusterWheatleyMultibodyDerivatives extends ClusterWheatleySoftDerivatives implements ClusterAbstractMultivalue {

    protected final MayerFunctionNonAdditive fNonAdditive;
    protected final MayerFunctionNonAdditive[] fMulti;
    protected final int[] moleculeIndices;
    protected final double[] r2;
    protected final MoleculeArrayList molecules;
    protected boolean doMulti;
    protected double rCut2;
    protected ClusterWheatleyMultibodyDerivativesBD clusterMultiBD;
    protected double multiTol;
    protected final double[] fQmulti;
    protected final int nDer;
    protected final boolean doTotal;
    public static boolean pushme = false, pushmeval=false;
    protected final double[] pairValues;
    protected long totcount = 0;
    protected long MultiBDcount = 0;

    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fMulti array of non-additive Mayer functions.  fMulti[3] is the
     *          3-body Mayer function (exp(-beta*deltaU3)-1), fMulti[4] is the
     *          4-body Mayer function, etc.  fMulti null entries will be
     *          ignored and the array need to not be of size equal to nPoints.
     *          If only 3-body Mayer function is available, then fMulti can be
     *          of length 4 (0,1,2,3).
     */
    public ClusterWheatleyMultibodyDerivatives(int nPoints, MayerFunction f, MayerFunctionNonAdditive[] fMulti, double multiTol, int nDer, boolean doTotal) {
        this(nPoints, f, null, fMulti, multiTol, nDer, doTotal);
    }

    /**
     * @param nPoints number of points
     * @param f pair Mayer function
     * @param fNonAdditive Mayer function that returns non-additive value for
     *          any number of molecules.
     */
    public ClusterWheatleyMultibodyDerivatives(int nPoints, MayerFunction f, MayerFunctionNonAdditive fNonAdditive, double multiTol, int nDer, boolean doTotal) {
        this(nPoints, f, fNonAdditive, new MayerFunctionNonAdditive[0], multiTol, nDer,doTotal);
    }

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
     *          of length 4 (0,1,2,3).
     * @param tol if the magnitude of the computed cluster value is less than
     *          the value will be recomputed using BigDecimal.  Use tol=0 to
     *          prevent BigDecimal computations.
     */
    public ClusterWheatleyMultibodyDerivatives(int nPoints, MayerFunction f, MayerFunctionNonAdditive fNonAdditive, MayerFunctionNonAdditive[] fMulti, double tol, int nDer, boolean doTotal) {
        super(nPoints, f, 0, nDer);
        this.nDer = nDer;
        this.doTotal=doTotal;
        this.fNonAdditive = fNonAdditive;
        this.fMulti = fMulti;        
        moleculeIndices = new int[nPoints];
        r2 = new double[nPoints*(nPoints-1)/2];
        // pairwise shouldn't try to use BD
        clusterBD = null;
        if(tol!=0){setTolerance(tol);}
        molecules = new MoleculeArrayList(nPoints);
        rCut2 = Double.POSITIVE_INFINITY;
        fQmulti = new double[1<<n];
        fQmulti[0] = fQmulti[1] = fQmulti[2] = 1;
        pairValues = new double[nDer+1];
        
    }

    public void setTolerance(double newTol) {
        clusterMultiBD = new ClusterWheatleyMultibodyDerivativesBD(n, f,fNonAdditive, fMulti, -3*(int)Math.log10(newTol),nDer,doTotal);
        clusterMultiBD.setDoCaching(false);
        clusterMultiBD.setPrecisionLimit(300);
        multiTol = newTol;
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleyMultibodyDerivatives c = new ClusterWheatleyMultibodyDerivatives(n, f, fNonAdditive, fMulti, multiTol,nDer,doTotal);
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
        return c;
    }

    public void setTemperature(double newT) {
        super.setTemperature(newT);
        if (clusterMultiBD != null) {
            clusterMultiBD.setTemperature(newT);
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
        
        double maxR2 = 0.1;
        if (pushme) {
            // force the system to hang out between minMaxR2 and maxMaxR2
            for (int i=0; i<n-1; i++) {
                for (int j=i+1; j<n; j++) {
                    double r2 = box.getCPairSet().getr2(i,j);
                    if (r2 > maxR2) maxR2 = r2;
                }
            }
            double minMaxR2 = 80*80;
            if (maxR2 < minMaxR2) {
                value[0] = 1e-200;
                return;
            }
            double maxMaxR2 = 100*100;
            if (maxR2 > maxMaxR2) {
                value[0] = 0;
                return;
            }
        }
        
        // do (multi+pair) - pair here so that we avoid recomputing f bonds
        
        if(!doTotal){
            doMulti = false;
            super.calcValue(box);
            System.arraycopy(value, 0, pairValues, 0, nDer+1);
        }
        doMulti = true;
        super.calcValue(box);
        if(!doTotal){
            for (int m=0; m<=nDer;m++){
                value[m] -= pairValues[m];
            }
        }
        // we have our own BD cluster for multibody
        // if an individual integrand (pair/total) is small, it won't trigger a BD calculation
        // BD only gets triggered here
        double bfac = (1.0-n)/SpecialFunctions.factorial(n);
        totcount++;
//        System.out.println(box.getIndex()+" "+java.util.Arrays.toString(value));
        if (Math.abs(value[0]) > 0 && Math.abs(value[0]/bfac) < multiTol) {
            if (clusterMultiBD != null) {
                double[] foo = value.clone();
                System.arraycopy(clusterMultiBD.getAllLastValues(box), 0, value, 0, nDer+1);
                MultiBDcount+=1;
                if(pushme||pushmeval){
                    for( int i=0;i<=nDer;i++){
                        System.out.print(value[i]/bfac+" "+foo[i]/bfac+" ");
                    }
                    System.out.println();
                }
            }
            else {
                for(int m=0;m<=nDer;m++){
                    value[m] = 0;
                }
            }
        }else if(pushmeval){value[0] = 1e-200;}
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
            fQmulti[i] = 1;
            int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)) {
                if (fQ[i][0] == 0){
                    for (int m=1;m<=nDer;m++){
                        fQ[i][m]=0;
                    }
                    continue;
                }                
                double c = Math.log(fQ[i][0])/beta;
                for(int m=1; m<=nDer;m++){
                    fQ[i][m]=fQ[i][m-1]*c;      //calculates derivatives of fQ w.r.t. beta      
                }                
                continue; // 2-point set; these fQ's were filled when bonds were computed, so skip
            }
            
            if (fQ[i][0] == 0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=0;
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
                fQ[i][0] *= fQmulti[isub];
                if (fQ[i][0]==0){
                    break;
                }
                
            }
            if (fQ[i][0]==0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=0;
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
                fQmulti[i] = 1;
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
                fQmulti[i] = fMulti[l].f(molecules, l, moleculeIndices, r2, beta)+1;
            }
            fQ[i][0] *= fQmulti[i];
            if (fNonAdditive != null) {
                // we don't want to include this in fQmulti because we would just include it again
                // for larger sets
                fQ[i][0] *= fNonAdditive.f(molecules, l, moleculeIndices, r2, beta)+1;
            }
            if (fQ[i][0] == 0){
                for (int m=1;m<=nDer;m++){
                    fQ[i][m]=0;
                }
                continue;
            }
            double c = Math.log(fQ[i][0])/beta;
            for(int m=1; m<=nDer;m++){
                fQ[i][m]=fQ[i][m-1]*c;      //calculates derivatives of fQ w.r.t. beta      
            }
            if(debugme)System.out.println(fQ[7][0]);
        }
    }
    
    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }
    
    public long getMultiBDcount(){
        return MultiBDcount;
    }
    
    public long getMulticount(){
        return totcount;
    }
    
    public double getMultiBDfrac(){
        return ((double)MultiBDcount)/totcount;
    }
    
    public static class ClusterRetrievePrimes implements ClusterAbstract{
        protected final ClusterWheatleyMultibodyDerivatives cluster;
        protected final int n;
        
        public ClusterRetrievePrimes(ClusterWheatleyMultibodyDerivatives cluster, int n){
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
