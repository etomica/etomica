package etomica.virial;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.model.impl.NodeImpl;
import etomica.graph.property.IsBiconnected;
import etomica.math.SpecialFunctions;

public class ClusterWheatleyExtendSW implements ClusterAbstract {
	
	protected final int n, nf, npairs;
    protected final MayerFunction f1, e2;
    protected final double[][] fQ;
    protected final double[] fQQ;
    protected final double[][] fC, fA, fB;
    protected double beta;
    protected IsBiconnected isB;
    protected Graph g;
    
    public ClusterWheatleyExtendSW(int nPoints, MayerFunction f1, MayerFunction e2){
    	this.n = nPoints;
    	this.nf = 1<<n; //2^n
    	this.npairs = n*(n-1)/2;
    	this.f1 = f1;
    	this.e2 = e2;
    	fQ = new double[nf][npairs+1];
    	for(int i=0; i<n; i++){
     		fQ[1<<i][0] = 1.0;
     		
//     		System.out.print("fQ["+(1<<i)+"]["+0+"]="+fQ[1<<i][0]+" ");
     		
     		for(int if1=1; if1<=npairs; if1++){
     			fQ[1<<i][if1] = 0.0;
     			
//     			System.out.print("fQ["+(1<<i)+"]["+if1+"]="+fQ[1<<i][if1]+" ");
     		}
//     		System.out.println();
     	}
    	fQQ = new double[npairs+1];
    	fC = new double[nf][npairs+1];
    	fA =  new double[nf][npairs+1];
    	fB =  new double[nf][npairs+1];
    
    	Node[] nds = new Node[n];
    	for(int i=0; i<nds.length; i++){
    		nds[i] = NodeImpl.createFieldNode((byte)i, 'b');
    	}
    	isB = new IsBiconnected();
    	g = new GraphImpl(nds);
    	outDegreeCore = new int[n];
    	outDegreeWell = new int[n];
    }
    
	public ClusterAbstract makeCopy() {
		ClusterWheatleyExtendSW c = new ClusterWheatleyExtendSW(n, f1, e2);
		c.setTemperature(1/beta);
		return c;
	}
	
	public int pointCount() {
		return n;
	}
	
	public double value(BoxCluster box) {
		return 0.0;
	}
	
	public double[] valueArray(BoxCluster box) {
		if (!checkConfig(box)) {
	    	for(int j=0; j<fB[nf-1].length; j++){
	    		fB[nf-1][j] = 0.0;
	    	}
	    	return fB[nf-1];
		}
		calcValue(box);
	    return fB[nf-1];
	}
	
	protected void updateF(BoxCluster box){
		CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        // recalculate all f values for all pairs
        for(int i=0; i<n; i++){
        	for(int j=i+1; j<n; j++){
        		fQ[1<<i|1<<j][0] = e2.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
        		fQ[1<<i|1<<j][1] = f1.f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
        		
//        		System.out.print("fQ["+(1<<i|1<<j)+"]["+ 0 +"]="+fQ[1<<i|1<<j][0]+" ");
//        		System.out.print("fQ["+(1<<i|1<<j)+"]["+ 1 +"]="+fQ[1<<i|1<<j][1]+" ");
        		
        		for(int if1=2; if1<=npairs; if1++){
        			fQ[1<<i|1<<j][if1] = 0.0;
        			
//        			System.out.print("fQ["+(1<<i|1<<j)+"]["+ if1 +"]="+fQ[1<<i|1<<j][if1]+" ");
        		}        			   
//        		System.out.println();
        	}
        }       
	}
	
	/**
     * This calculates all FQ values given that the entries for pairs have
     * already been populated.
     */
	protected void calcFullFQ(BoxCluster box){
		for(int i=3; i<nf; i++){
			int j = i & -i;
			if (i==j) continue;
			int k = i & ~j;
			if (k == (k & -k)) continue;

	        for(int if1=0; if1<=npairs; if1++){//if1=index of f1 bond, when #of f1 bond1=0, 1, 2, ... NPAIRS
	        	fQ[i][if1] = fQ[k][if1];
            }

            for (int l=(j<<1); l<i; l=(l<<1)){//l=index of 1 point set
            	if ( (l&i) == 0 ) continue;
                for(int if1=0; if1<=npairs; if1++){
                	fQQ[if1] = fQ[i][if1];
                    fQ[i][if1] = fQ[i][if1] * fQ[l|j][0];//add a e2 bonds.
                    if(if1>0) fQ[i][if1] += fQQ[if1-1] * fQ[l|j][1];//add a f1 bonds
                }
            }
		}
//		System.out.println("fQ="); MyUtility.display2DArray(fQ);
	}
	
	/**
     * Returns the cluster value for the given configuration.  
     */
    public double[] calcValue(BoxCluster box) {
        calcFullFQ(box);
       //Compute the fC's
       for (int i=1; i<nf; i++){
    	   for (int if1=0; if1<=npairs; if1++){
    		   fC[i][if1] = fQ[i][if1];
           }
           int iLowBit = i & -i;
           int inc = iLowBit<<1;
           
           for (int j=iLowBit; j<i; j+=inc){
        	   int jComp = i & ~j;
               while ((j|jComp) != i && j<i){
            	   int jHighBits = j^iLowBit;
                   int jlow = jHighBits & -jHighBits;
                   j += jlow;
                   jComp = (i & ~j);
               }
               if (j==i) break;
               for (int if1=0; if1<=npairs; if1++){
            	   for (int k=0; k<=if1; k++){
            		   fC[i][if1] -= fC[j][k] * fQ[jComp][if1-k];
                   }  
               }
           }
       }
//       System.out.println("fC="); MyUtility.display2DArray(fC);

       //find fA1
       for (int i=2; i<nf; i+=2){//all even sets (2,4,6) don't contain 1
    	   for (int if1=0; if1<=npairs; if1++){
    		   fB[i][if1] = fC[i][if1];
//    		   System.out.println("fB[" + i + "][" + if1+ "] =" + fB[i][if1]);
           }
       }

       for (int if1=0; if1<=npairs; if1++){
    	   fA[1][if1] = 0;
           fB[1][if1] = fC[1][if1];
//           System.out.println("fB[" + 1 + "][" + if1+ "] =" + fB[1][if1]);
       }

       for (int i=3; i<nf; i+=2){//every set will contain 1.
    	   for (int if1=0; if1<=npairs; if1++){
    		   fA[i][if1] = 0;
               fB[i][if1] = fC[i][if1];
           }
           int ii = i - 1;//all bits in i but lowest
           int iLow2Bit = (ii & -ii);//next lowest bit
           int jBits = 1 | iLow2Bit;
           if (jBits == i) continue;

           int iii = ii ^ iLow2Bit; //i with 2 lowest bits off
           int jInc = (iii & -iii);//3rd lowest bit, alsso increment for j
           for (int j=jBits; j<i; j+=jInc){//sum over partitions of i containing j Bits
        	   int jComp = (i & ~j);//subset of i complementing j
               while ((j|jComp) != i && j<i){//if j is not a proper subset of i.
            	   int jHighBits = j ^ jBits;
                   int jlow = jHighBits & -jHighBits;
                   j += jlow;
                   jComp = (i & ~j);
               }
               if (j==i) break;
//               if(i==nf-1) System.out.println("fA["+i+"]["+0+"]="+fA[i][0]+",fB["+i+"]["+0+"]="+fB[i][0]); 
               for (int if1=0; if1<=npairs; if1++){
            	   for (int k=0; k<=if1; k++){
            		   fA[i][if1] += fB[j][k] * fC[jComp|1][if1-k];
//                     if((i==NF-1) && (if1==0)) printf("* %f, %f, %f\n", fA[i][if1], fB[j][k], fC[jComp|1][if1-k]);
//            		   if((i==nf-1) && (if1==0)) System.out.println("fA["+i+"]["+if1+"]="+fA[i][if1]+",[if1"+",fB["+j+"]["+k+"]="+fB[j][k]+",fC["+(jComp|1)+"]["+(if1-k)+"]="+fC[jComp|1][if1-k]);
            	   }
			    }
           }	
           for (int if1=0; if1<=npairs; if1++){
        	   fB[i][if1] -= fA[i][if1];//remove from B graphs that contain articulation point 0.
           }
//           if(i==nf-1) System.out.println("fB["+i+"]["+0+"]="+fB[i][0]); 
       }

       for (int v=1; v<n; v++){
    	   int vs1 = 1<<v;
           for (int i=vs1+1; i<nf; i++){
        	   for (int if1=0; if1<=npairs; if1++){
        		   fA[i][if1] = 0;
               }
               if ( (i & vs1) == 0 ) continue;
               int iLowBit = (i & -i);
               if ( iLowBit == i ) continue;

               int jBits;
               int ii = i ^ iLowBit;
               int iLow2Bit = (ii & -ii);
               if ( iLowBit!=vs1 && iLow2Bit!=vs1 ){
            	   jBits = iLowBit | vs1;//v is not in the lowest 2 bits
                   int jInc = iLow2Bit;    //we can only increment by the 2nd lowest
                   for (int j=jBits; j<i; j+=jInc){
                	   if ( (j&jBits) != jBits ){
                		   j |= vs1;
                           if (j==i) break;
                       }
                       int jComp = i & ~j;
                       while ((j|jComp) != i && j<i){
                    	   int jHighBits = j^jBits;
                           int jlow = jHighBits & -jHighBits;
                           j += jlow;
                           j |= vs1;
                           jComp = (i & ~j);
                       }
                       if (j==i) break;
                       for (int if1=0; if1<=npairs; if1++){
                    	   for (int k=0; k<=if1; k++){
                    		   fA[i][if1] += fB[j][k] * (fB[jComp|vs1][if1-k] + fA[jComp|vs1][if1-k]);
                           }
                       }
                   }
               }else{
            	   //lowest 2 bits contain v
                   jBits = iLowBit | iLow2Bit;
                   if (jBits == i) continue; // no bits left jComp

                   int iii = ii ^ iLow2Bit;
                   int jInc = ( iii & -iii);
                   //at this point jBits has (lowest bit + v)
                   for (int j=jBits; j<i; j+=jInc){//sum over partitions of i
                	   int jComp = i & ~j;
                       while ((j|jComp) != i && j<i){
                    	   int jHighBits = j^jBits;
                           int jlow = jHighBits & -jHighBits;
                           j += jlow;
                           jComp = (i & ~j);
                       }
                       if (j==i) break;
                       for (int if1=0; if1<=npairs; if1++){
                    	   for (int k=0; k<=if1; k++){
                    		   fA[i][if1] += fB[j][k] * (fB[jComp|vs1][if1-k] + fA[jComp|vs1][if1-k]);
                           }
                       }
                   }
               }
               for (int if1=0; if1<=npairs; if1++){
            	   fB[i][if1] -= fA[i][if1];//remove from B graphs
               }
//        	   if(i==nf-1) System.out.println("fB["+i+"]["+0+"]=" + fB[i][0]);
           }
       }
       // we're actually messing up fB here, but we don't need it anymore
       double[] rv = fB[nf-1];
       double fac = (1.0-n)/SpecialFunctions.factorial(n);
       for (int i=0; i<rv.length; i++) {
    	   rv[i] *= fac;
       }
       return rv;
//       System.out.println("fA="); MyUtility.display2DArray(fA);
//       System.out.println("fB="); MyUtility.display2DArray(fB);
    }
    
	public void setTemperature(double temperature) {
		beta = 1/temperature;
	}

	public int getCoreEdgeCount() {
		return edgeCountCore;
	}

	public int getWellEdgeCount() {
		return edgeCountWell;
	}

	protected int edgeCountCore, edgeCountWell;
	protected final int[] outDegreeCore, outDegreeWell;

	public int[] getOutDegreeCore() {
		return outDegreeCore;
	}

	public int[] getOutDegreeWell() {
		return outDegreeWell;
	}

	public boolean checkConfig(BoxCluster box) {
	    updateF(box);
	    edgeCountCore = edgeCountWell = 0;
	    for(int i=0; i<n; i++){
        	for(int j=i+1; j<n; j++){
        		double e2 = fQ[1<<i|1<<j][0];
        		double f1 = fQ[1<<i|1<<j][1];
        		if (e2==0 || f1>0){//Means when there is an edge.
        			g.putEdge((byte)i, (byte)j);
        			if (e2==1) {
        				edgeCountWell++;
        				outDegreeWell[i]++;
        				outDegreeWell[j]++;
        			}
        			else {
        				edgeCountCore++;
        				outDegreeCore[i]++;
        				outDegreeCore[j]++;
        			}
        		}else if(g.hasEdge((byte)i, (byte)j)){
        			g.deleteEdge((byte)i, (byte)j);//This is remove edges set by last step. It is a initilization. 
        		}
        	}
	    }
	    return isB.check(g);
	}
}
