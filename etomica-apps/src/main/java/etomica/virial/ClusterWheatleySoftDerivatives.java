/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.math.SpecialFunctions;
import etomica.util.random.IRandom;


/**
 * This class calculates the sum of all biconnected clusters using Wheatley's
 * recursive formulation.
 * 
 * @author David Kofke and Andrew Schultz 
 */
public class ClusterWheatleySoftDerivatives implements ClusterAbstract, ClusterAbstractMultivalue {

    protected final int n;
    protected final MayerFunction f;
    protected final int nDer;
    protected final double[][] fQ, fC;
    protected final double[][] fA, fB;
    protected long cPairID = -1, lastCPairID = -1;
    protected final double[] value, lastValue;
    protected double beta;
    public static boolean pushme = false, pushmeval=false;
    protected double tol;
    protected ClusterWheatleySoftDerivativesBD clusterBD;
    protected boolean debug = false;
    protected boolean doCaching = true;
    protected final int[][] binomial;
    protected long SoftBDcount = 0;
    protected double BDAccFrac = 1;
    protected IRandom random;
    protected int stepcount = 0;
    protected long totcount = 0;
    protected boolean count=false;
    protected double rCut2 = Double.POSITIVE_INFINITY;


    public ClusterWheatleySoftDerivatives(int nPoints, MayerFunction f, double tol, int nDer) {        
        this.n = nPoints;
        value = new double[nDer+1];
        lastValue = new double[nDer+1];
        this.f = f;
        int nf = 1<<n;  // 2^n
        fQ = new double[nf][nDer+1];
        fC = new double[nf][nDer+1];
        for(int i=0; i<n; i++) {
            fQ[1<<i][0] = 1.0;
    	}
        fA = new double[nf][nDer+1];
        fB = new double[nf][nDer+1];
        this.nDer = nDer;
        setTolerance(tol);
        this.binomial = new int[nDer+1][]; 
        for(int m=0;m<=nDer;m++){
            binomial[m] = new int[m+1];
            for(int l=0;l<=m;l++){
                binomial[m][l] = (int)(SpecialFunctions.factorial(m)/(SpecialFunctions.factorial(l)*SpecialFunctions.factorial(m-l)));
            }
        }
    }

    public void setTolerance(double newTol) {
        if (newTol > 0) {
            clusterBD = new ClusterWheatleySoftDerivativesBD(n, f, -3*(int)Math.log10(newTol),nDer);
            clusterBD.setDoCaching(false);
            clusterBD.setPrecisionLimit(300);
        } else {
            clusterBD = null;
        }
        tol = newTol;
    }

    /**
     * Directs this cluster to only compute p fraction of the time when the
     * value is too small (below tol).  When it is computed, the value will be
     * boosted by 1/p.
     *
     * @param p   the fraction of time BD values will be computed
     * @param rng the random number generated used to decide to do BD or not
     */
    public void setBDAccFrac(double p, IRandom rng) {
        BDAccFrac = p;
        random = rng;
    }


    public void setDoCaching(boolean newDoCaching) {
        doCaching = newDoCaching;
        if (clusterBD != null) {
            clusterBD.setDoCaching(doCaching);
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterWheatleySoftDerivatives c = new ClusterWheatleySoftDerivatives(n, f, tol, nDer);
        c.setTemperature(1/beta);
        c.setDoCaching(doCaching);
        if (BDAccFrac < 1) {
            c.setBDAccFrac(BDAccFrac, random);
        }
        return c;
    }

    public int pointCount() {
        return n;
    }

    public double value(BoxCluster box) {
        if (doCaching) {
            CoordinatePairSet cPairs = box.getCPairSet();
            long thisCPairID = cPairs.getID();
            if (thisCPairID == cPairID) {
                return value[0];
            }
            if (thisCPairID == lastCPairID) {
                // we went back to the previous cluster, presumably because the last
                // cluster was a trial that was rejected.  so drop the most recent value/ID
                cPairID = lastCPairID;
                for(int m=0;m<=nDer;m++){
                	value[m] = lastValue[m];
                }
                return value[0];
            }
    
            // a new cluster
            lastCPairID = cPairID;
            for(int m=0;m<=nDer;m++){
            	lastValue[m] = value[m];
            }
            cPairID = thisCPairID;
        }
        stepcount++;
        updateF(box);
      
        calcValue(box);
        if (Double.isNaN(value[0]) || Double.isInfinite(value[0])) {
            debug = true;
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops "+value[0]);
        }
        return value[0];
    }

    /**
     * This calculates all FQ values given that the entries for pairs have
     * already been populated.
     */
    protected void calcFullFQ(BoxCluster box) {
        int nf = 1<<n;
        // generate all partitions and compute product of e-bonds for all pairs in partition
        for (int i=3; i<nf; i++) {
        	
        	int j = i & -i;//lowest bit in i
            if (i==j) continue; // 1-point set
            int k = i&~j; //strip j bit from i and set result to k
            if (k == (k&-k)){
                // 2-point set; these fQ's were filled when bonds were computed, so skip
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
            fQ[i][0] = fQ[k][0]; //initialize with previously-computed product of all pairs in partition, other than j
            if (fQ[i][0] == 0){
            	for (int m=1;m<=nDer;m++){
            		fQ[i][m]=0;
            	}
            	continue;
            }
            //loop over pairs formed from j and each point in set i; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l=(j<<1); l<i; l=(l<<1)) {
                if ((l&i)==0) continue; //l is not in i
                fQ[i][0] *= fQ[l | j][0];
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
        }
    }

    public void setRCut(double newRCut) {
        rCut2 = newRCut * newRCut;
    }

    /**
     * Returns the cluster value for the given configuration.  You must call
     * doCheck(BoxCluster) before calling this method.
     */
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
            double minMaxR2 = 5*5;
            if (maxR2 < minMaxR2) {
                value[0] = 1e-200;
                return;
            }
            double maxMaxR2 = 7*7;
            if (maxR2 > maxMaxR2) {
                value[0] = 0;
                return;
            }
        }
        calcFullFQ(box);
        int nf = 1<<n;
        //Compute the fC's
        for(int i=1; i<nf; i++) {
        	for(int m=0;m<=nDer;m++){
        		fC[i][m] = fQ[i][m];
//        		System.out.println("fC["+i+"]["+m+"]"+fC[i][m]);
        	}
            int iLowBit = i & -i;
            int inc = iLowBit<<1;
            for(int j=iLowBit; j<i; j+=inc) {
                int jComp = i & ~j;
                while ((j|jComp) != i && j<i) {
                    int jHighBits = j^iLowBit;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;

//                System.out.println("jcomp = "+jComp);
                for(int m=0;m<=nDer;m++){                
                	for(int l=0;l<=m;l++){
//
                		fC[i][m] -= binomial[m][l]*fC[j][l] * fQ[jComp][m-l];//for fQ, flip the bits on j; use only those appearing in i
                		//computes fC and its derivatives w.r.t beta
                	}
                }
            }
        }

        // find fA1
        for (int i=2; i<nf; i+=2) {
            // all even sets don't contain 1
            //fA[i] = 0;
        	for(int m=0;m<=nDer;m++){
        		fB[i][m] = fC[i][m];
        	}
        }
        for(int m=0;m<=nDer;m++){
        	fA[1][m] = 0;
        	fB[1][m] = fC[1][m];
        }
        for (int i=3; i<nf; i+=2) {
            // every set will contain 1
        	for(int m=0;m<=nDer;m++){
        		fA[i][m] = 0;
        		fB[i][m] = fC[i][m];
        	}
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            if (jBits==i) continue;
            //jBits has 1 and next lowest bit in i
            int iii = ii ^ iLow2Bit;//i with 2 lowest bits off
            int jInc = (iii & -iii);//3rd lowest bit, also increment for j
            for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i containing jBits
                int jComp = (i & ~j); //subset of i complementing j
                while ((j|jComp) != i && j<i) {
                    int jHighBits = j^jBits;
                    int jlow = jHighBits & -jHighBits;
                    j += jlow;
                    jComp = (i & ~j);
                }
                if (j==i) break;
                for(int m=0;m<=nDer;m++){
                	for(int l=0;l<=m;l++){
                		fA[i][m] += binomial[m][l]*fB[j][l] * fC[jComp|1][m-l];
                	}
                }
            }
            for(int m=0;m<=nDer;m++){
            	fB[i][m] -= fA[i][m];//remove from B graphs that contain articulation point at 0
            }
        }

        for (int v=1; v<n; v++) {
            int vs1 = 1<<v;
            for (int i=vs1+1; i<nf; i++) {
            	for(int m=0;m<=nDer;m++){
            		fA[i][m] = 0;
            	}
//                fB[v][i] = fB[v-1][i];//no a.p. at v or below, starts with those having no a.p. at v-1 or below
                //rest of this is to generate A (diagrams having a.p. at v but not below), and subtract it from B
                if ((i & vs1) == 0) continue;//if i doesn't contain v, fA and fB are done
                int iLowBit = (i&-i);//lowest bit in i
                if (iLowBit == i) { //lowest bit is only bit; fA and fB are done
                    continue;
                }
                int jBits;
                int ii = i ^ iLowBit;
                int iLow2Bit = (ii & -ii);
                if (iLowBit != vs1 && iLow2Bit != vs1) {
                    //v is not in the lowest 2 bits
                    // jBits is the lowest bit and v
                    jBits = iLowBit | vs1;

                    // we can only increment by the 2nd lowest
                    int jInc = iLow2Bit;

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        if ((j & jBits) != jBits) {
                            //ensure jBits are in j
                            j |= vs1;
                            if (j==i) break;
                        }
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow; // this might knock out the v bit
                            j |= vs1;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        for(int m=0;m<=nDer;m++){                
                            for(int l=0;l<=m;l++){
                                fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
                            }
                        }
                    }
                }
                else {
                    //lowest 2 bits contain v
                    // jBits is the lowest 2 bits
                    // we can start at jBits and increment by the 3rd lowest bit
                    jBits = iLowBit | iLow2Bit;
                    if (jBits == i) continue; // no bits left for jComp
                    
                    int iii = ii ^ iLow2Bit;
                    int jInc = (iii & -iii);

                    //at this point jBits has (lowest bit + v) or (v + next lowest bit)
                    for (int j=jBits; j<i; j+=jInc) {//sum over partitions of i
                        // start=jBits and jInc ensure that every set includes jBits
                        int jComp = i & ~j;//subset of i complementing j
                        while ((j|jComp) != i && j<i) {
                            int jHighBits = j^jBits;
                            int jlow = jHighBits & -jHighBits;
                            j += jlow;
                            jComp = (i & ~j);
                        }
                        if (j==i) break;
                        for(int m=0;m<=nDer;m++){
                        	for(int l=0;l<=m;l++){
                        		fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
                        	}
                        }
                    }
                }
                for(int m=0;m<=nDer;m++){
                	fB[i][m] -= fA[i][m];//remove from B graphs that contain articulation point at v
                }
            }
        }
        double bfac = (1.0-n)/SpecialFunctions.factorial(n);     
        totcount++;
        double r = BDAccFrac < 1 ? random.nextDouble() : 1;
        if (Math.abs(fB[nf - 1][0]) < Math.abs(tol)) {
            // integrand is too small for recursion to compute accurately.  we ought to do
            // BD, but it's expensive.  only do BD BDAccFrac of the time.  If we do it, then
            // boost the returned value by 1/BDAccFrac to account for the missed configurations
            boolean doBD = clusterBD != null && (BDAccFrac == 1 || r < BDAccFrac);
            boolean returnDoubleVal = false;
            if (doBD) {
                for(int m = 0; m<=nDer; m++) {
                    value[m] = bfac * fB[nf - 1][m] / BDAccFrac;
                }
                double[] foo = value.clone();                            
                SoftBDcount+=1;                
                if(count&&box.getIndex()==1&&!returnDoubleVal)System.out.println(stepcount);                
                System.arraycopy(clusterBD.getAllLastValues(box), 0, value, 0, nDer+1);                
                if(count&&box.getIndex()==1 && false){
                    double err = Math.abs(Math.abs(value[0]/bfac)-Math.abs(foo[0]/bfac))*100/Math.abs(value[0]/bfac);
                    System.out.println(stepcount + " : fB double = " + Math.abs(foo[0]/bfac)+" ,fB BD = " + Math.abs(value[0]/bfac) + ", Error % = "+ err);
                    if(true){
                        double minr2 = 1000;
                        int mini=1000;
                        int minj=1000;
                        for (int i=0; i<n-1; i++) {
                            for (int j=i+1; j<n; j++) {
                                double r2 = box.getCPairSet().getr2(i,j);
                                if (r2<minr2){
                                    minr2=r2;
                                    mini=i;
                                    minj=j;
                                }
                            }
                        }System.out.println("R2 pair "+mini+","+minj+ " : " + minr2);
                    }
                }
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
            if (!returnDoubleVal) return;
        }

//        System.out.println("fQ"+" "+Arrays.toString(fQ[nf-1]));
//        for(int i=0;i<fC.length;i++){
//            System.out.println("fC"+" "+i+" "+Arrays.toString(fC[i]));}
//        System.out.println("fA"+" "+Arrays.toString(fA[nf-1]));
//        System.out.println("fB"+" "+Arrays.toString(fB[nf-1]));
//        if (fB[nf-1][0] != -1 && fB[nf-1][1] == 0) {
//            throw new RuntimeException("oops");
//        }
        for(int m=0;m<=nDer;m++){
            value[m] = bfac*fB[nf-1][m];
        }
        if (pushme && maxR2 > 2*2) {
//            value *= Math.pow(maxR2/4, 6);
        }
        if(pushmeval){ value[0] = 1e-200;}
    }

    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        f.setBox(box);
        // recalculate all f values for all pairs
        for(int i=0; i<n-1; i++) {
            for(int j=i+1; j<n; j++) {
                double ff = f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                if (false && Double.isNaN(ff)) {
                    f.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                    throw new RuntimeException("oops");
                }
//                if (Math.abs(ff) < 1e-14) ff = 0;
                fQ[(1<<i)|(1<<j)][0] = ff+1;
            }
        }
        
//        ANALYTICAL CHECK FOR B2 LJ
//        double sum = 0;
//        double sumpi = 0;
//        double dsum = 0;
//        double dr = 0.001;
//        for (int i=0; i<10000; i++) {
//            double r = dr*i;
//            double fval = f.f(null,  r*r, beta);
//            sum += fval*r*r*dr;
//            sumpi += Math.abs(fval)*r*r*dr;
//            double u = ((P2LennardJones)f.getPotential()).u(r*r);
//            if (u<Double.POSITIVE_INFINITY) {
//                dsum += -u*Math.exp(-u*beta)*r*r*dr;
//            }
//        }
//        System.out.println("B2 = "+(-0.5*4*Math.PI*sum) +" "+ (-0.5*4*Math.PI*dsum)+" "+(dsum/sumpi)+" "+2*Math.PI*sumpi);
    }

    public void setTemperature(double temperature) {
        beta = 1/temperature;
        if (clusterBD != null) {
            clusterBD.setTemperature(temperature);
        }
    }

    public int getNumValues() {
        return value.length;
    }

    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }
    
    public long getSoftBDcount(){
        return SoftBDcount;
    }
    
    public long getSoftcount(){
        return totcount;
    }
    
    public double getSoftBDfrac(){
        return ((double)SoftBDcount)/totcount;
    }
    
    public static class ClusterRetrievePrimes implements ClusterAbstract{
    	protected final ClusterWheatleySoftDerivatives cluster;
    	protected final int n;
    	
    	public ClusterRetrievePrimes(ClusterWheatleySoftDerivatives cluster, int n){
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
