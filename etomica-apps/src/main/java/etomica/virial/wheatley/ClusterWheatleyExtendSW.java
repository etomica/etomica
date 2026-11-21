/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.model.impl.NodeImpl;
import etomica.graph.property.IsBiconnected;
import etomica.math.SpecialFunctions;
import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.ClusterAbstract;

/**
 * Computes the integrand of the <i>n</i>-th virial coefficient of a square-well
 * fluid as a polynomial in the well variable, using Wheatley's recursion for the sum
 * over all biconnected Mayer diagrams. One simulation gives <i>B</i><sub>n</sub> at
 * every temperature.
 *
 * <h2>What is computed</h2>
 * For a square-well potential (hard core of diameter &sigma;, well of depth &epsilon;
 * out to &lambda;&sigma;), each Mayer e-bond can be written
 * <blockquote>
 *   <i>e</i>(<i>r</i>) = <i>e</i><sub>2</sub>(<i>r</i>) + <i>x</i> <i>f</i><sub>1</sub>(<i>r</i>),
 *   &nbsp;&nbsp; <i>x</i> = exp(&beta;&epsilon;) &minus; 1,
 * </blockquote>
 * where <i>e</i><sub>2</sub> is the hard-core e-bond (0 inside the core, 1 outside)
 * and <i>f</i><sub>1</sub> is 1 for pairs inside the well and 0 otherwise. The sum
 * over biconnected diagrams on <i>n</i> points is then a polynomial in <i>x</i> of
 * degree at most <i>n</i>(<i>n</i> &minus; 1)/2, the number of pairs. For one
 * configuration, {@link #valueArray(BoxCluster)} returns its coefficients,
 * <blockquote>
 *   {@code value[m]} = (1 &minus; <i>n</i>)/<i>n</i>! &times;
 *   [<i>x</i><sup>m</sup>] <i>f</i><sub>B</sub>(<i>x</i>),
 *   &nbsp;&nbsp; <i>m</i> = 0, 1, &hellip;, <i>n</i>(<i>n</i> &minus; 1)/2,
 * </blockquote>
 * where [<i>x</i><sup>m</sup>] denotes the coefficient of <i>x</i><sup>m</sup>.
 * {@code value[m]} collects the diagrams with exactly <i>m</i> well bonds, so
 * {@code value[0]} is the hard-sphere integrand with diameter &sigma;. Averaging each
 * coefficient by Mayer sampling and evaluating
 * &Sigma;<sub>m</sub> {@code value[m]} <i>x</i><sup>m</sup> gives
 * <i>B</i><sub>n</sub>(<i>T</i>) for any temperature. No truncation is involved: the
 * polynomial is exact. Drivers include
 * {@link etomica.virial.simulations.squarewell.VirialSWWE} (with
 * {@link etomica.virial.MeterVirialSWWE}) and
 * {@link etomica.virial.simulations.squarewell.VirialSQWBinMultiThreaded}.
 *
 * <h2>Algorithm</h2>
 * For each subset <i>S</i> of the points, the product of e-bonds over the pairs in
 * <i>S</i> is a polynomial in <i>x</i>; {@link #calcFullFQ(BoxCluster)} builds its
 * coefficients by multiplying in one pair at a time. These polynomials are propagated
 * through the recursion <i>f</i><sub>Q</sub> &rarr; <i>f</i><sub>C</sub> (connected)
 * &rarr; <i>f</i><sub>B</sub> (biconnected) of
 * R. J. Wheatley, Phys. Rev. Lett. <b>110</b>, 200601 (2013), in which every product
 * of two subset quantities becomes a product of polynomials (a convolution of their
 * coefficients). The cost is
 * O(3<sup><i>n</i></sup> <i>P</i><sup>2</sup>) per configuration, where
 * <i>P</i> = <i>n</i>(<i>n</i> &minus; 1)/2, and memory is
 * O(2<sup><i>n</i></sup> <i>P</i>). With {@code f1} and {@code e2} equal to 0 or 1,
 * every coefficient is an integer, so the double-precision results are exact as long
 * as the intermediate values stay below 2<sup>53</sup>.
 *
 * <h2>Screening and configuration statistics</h2>
 * {@link #valueArray(BoxCluster)} first calls {@link #checkConfig(BoxCluster)}. It
 * builds the graph whose edges are the pairs with a non-zero Mayer function at some
 * <i>x</i> (core overlaps and well pairs). If that graph is not biconnected, every
 * coefficient is zero and the recursion is skipped. {@code checkConfig} also records
 * statistics used by binning meters such as
 * {@link etomica.virial.MeterVirialEBinMultiThreaded}: the numbers of core and well
 * edges, each point's core and well degree, and lists of the cliques (sets of three or
 * more points that are not the full set) formed by well edges, by core edges, by edges
 * of either kind, and by non-edges. These cliques are recorded, not used for
 * screening.
 *
 * <h2>Inputs and conventions</h2>
 * <ul>
 * <li>{@code e2} must return the hard-core e-bond, 0 or 1 (for example
 *     {@link etomica.virial.MayerEHardSphere}).</li>
 * <li>{@code f1} must return 1 for a pair inside the well and 0 otherwise. Its value
 *     must not depend on temperature, because the temperature dependence is carried
 *     entirely by <i>x</i>.</li>
 * <li>The edge classification in {@code checkConfig} assumes these values: a pair is a
 *     core edge if {@code e2} is 0, a well edge if {@code e2} is 1 and {@code f1} is
 *     positive.</li>
 * </ul>
 *
 * <h2>Limitations</h2>
 * <ul>
 * <li>{@link #value(BoxCluster)} always returns 0. Use
 *     {@link #valueArray(BoxCluster)}; this cluster cannot be sampled directly as a
 *     target cluster.</li>
 * <li>There is no caching: every call to {@code valueArray} recomputes.</li>
 * <li>The array returned by {@code valueArray} is an internal array that is
 *     overwritten by the next call.</li>
 * <li>The Mayer functions' {@code setBox} methods are not called.</li>
 * </ul>
 *
 * @see ClusterWheatleyPT
 */
public class ClusterWheatleyExtendSW implements ClusterAbstract {

	/** Number of points (molecules), n. */
	protected final int n;
	/** Number of subsets of the points, 2<sup>n</sup>. */
	protected final int nf;
	/** Number of pairs, P = n(n &minus; 1)/2; the polynomials have degree at most P. */
	protected final int npairs;
	/** Well indicator: 1 for a pair inside the well, 0 otherwise. */
	protected final MayerFunction f1;
	/** Hard-core e-bond: 0 inside the core, 1 outside. */
	protected final MayerFunction e2;
	/**
	 * Polynomial coefficients of the product of e-bonds for each subset, indexed
	 * [subset bitmask][power of x]. Pair entries are set by {@link #updateF(BoxCluster)};
	 * entries for larger subsets are set by {@link #calcFullFQ(BoxCluster)}.
	 */
	protected final double[][] fQ;
	/** Work array used by {@link #calcFullFQ(BoxCluster)}. */
	protected final double[] fQQ;
	/** Polynomial coefficients of the connected-diagram sum f<sub>C</sub> for each subset. */
	protected final double[][] fC;
	/** Work array for the articulation-point recursion, indexed like {@link #fB}. */
	protected final double[][] fA;
	/**
	 * Polynomial coefficients of the biconnected-diagram sum f<sub>B</sub> for each
	 * subset. {@code fB[nf-1]} holds the scaled result after {@link #calcValue(BoxCluster)}.
	 */
	protected final double[][] fB;
	/** Work array holding f<sub>B</sub> + f<sub>A</sub> for one subset during the recursion. */
	protected final double[] fBA;
	/** Reciprocal temperature, passed to the Mayer functions (which should ignore it). */
	protected double beta;
	/** Biconnectivity test used by {@link #checkConfig(BoxCluster)}. */
	protected IsBiconnected isB;
	/** Graph of core and well edges for the current configuration. */
	protected Graph g;

    /**
     * Creates the cluster.
     *
     * @param nPoints number of points (molecules), n; each work array has
     *                2<sup>n</sup> rows of n(n &minus; 1)/2 + 1 entries
     * @param f1      well indicator: 1 for a pair inside the well, 0 otherwise
     * @param e2      hard-core e-bond: 0 inside the core, 1 outside
     */
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
    	fBA = new double[npairs+1];

    	Node[] nds = new Node[n];
    	for(int i=0; i<nds.length; i++){
    		nds[i] = NodeImpl.createFieldNode((byte)i, 'b');
    	}
    	isB = new IsBiconnected();
    	g = new GraphImpl(nds);
    	outDegreeCore = new int[n];
    	outDegreeWell = new int[n];
        cliqueSetF1 = new boolean[nf];
        cliqueSetE2 = new boolean[nf];
        cliqueSetEF = new boolean[nf];
        cliqueSetNone = new boolean[nf];
        cliqueListF1 = new int[nf];
        cliqueListE2 = new int[nf];
        cliqueListEF = new int[nf];
        cliqueListNone = new int[nf];
    }

	/** Returns a new cluster with the same points, Mayer functions and temperature. */
	public ClusterAbstract makeCopy() {
		ClusterWheatleyExtendSW c = new ClusterWheatleyExtendSW(n, f1, e2);
		c.setTemperature(1/beta);
		return c;
	}

	/** Returns the number of points, n. */
	public int pointCount() {
		return n;
	}

	/**
	 * Always returns 0; use {@link #valueArray(BoxCluster)} instead.
	 *
	 * @param box the box holding the configuration (not used)
	 * @return 0
	 */
	public double value(BoxCluster box) {
		return 0.0;
	}

	/**
	 * Returns the polynomial coefficients for the box's current configuration. Calls
	 * {@link #checkConfig(BoxCluster)} first, and returns all zeros if the
	 * configuration is screened out.
	 *
	 * @param box the box holding the configuration
	 * @return the coefficients, indexed by the power of x (the number of well bonds),
	 *         of length n(n &minus; 1)/2 + 1; this is an internal array that the next
	 *         call overwrites
	 */
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

	/**
	 * Evaluates both Mayer functions for every pair and sets the pair entries of
	 * {@link #fQ}: coefficient 0 is {@code e2}, coefficient 1 is {@code f1}, and the
	 * rest are 0.
	 *
	 * @param box the box holding the configuration
	 */
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
	 * Fills {@link #fQ} for every subset of three or more points, given that the pair
	 * entries have already been set by {@link #updateF(BoxCluster)}. Each subset's
	 * polynomial is the previous subset's polynomial multiplied by
	 * (e<sub>2</sub> + x f<sub>1</sub>) for each added pair.
	 *
	 * @param box the box holding the configuration (not used)
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
	 * Runs the recursion for the current configuration and returns the polynomial
	 * coefficients. {@link #updateF(BoxCluster)} must have been called first for this
	 * configuration (as {@link #checkConfig(BoxCluster)} does). Most callers should use
	 * {@link #valueArray(BoxCluster)}.
	 *
	 * <p>Each sum over subsets visits only the proper subsets j of a set i that contain
	 * a fixed set of bits, using the submask step {@code s = (s - rest) & rest}.
	 *
	 * @param box the box holding the configuration (not used directly)
	 * @return the coefficients; this is {@code fB[nf-1]}, scaled in place by
	 *         (1 &minus; n)/n!
	 */
    public double[] calcValue(BoxCluster box) {
        calcFullFQ(box);
        //Compute the fC's
        for (int i=1; i<nf; i++){
            double[] fCi = fC[i];
            System.arraycopy(fQ[i], 0, fCi, 0, npairs+1);
            int iLowBit = i & -i;
            int rest = i ^ iLowBit;
            // sum over proper subsets j of i that contain iLowBit; jComp = i & ~j
            for (int s=0; s!=rest; s=(s-rest)&rest){
                double[] fCj = fC[iLowBit|s];
                double[] fQjComp = fQ[rest^s];
                for (int if1=0; if1<=npairs; if1++){
                    double t = 0;
                    for (int k=0; k<=if1; k++){
                        t += fCj[k] * fQjComp[if1-k];
                    }
                    fCi[if1] -= t;
                }
            }
        }

        //find fA1
        for (int i=2; i<nf; i+=2){//all even sets (2,4,6) don't contain 1
            System.arraycopy(fC[i], 0, fB[i], 0, npairs+1);
        }

        for (int if1=0; if1<=npairs; if1++){
            fA[1][if1] = 0;
            fB[1][if1] = fC[1][if1];
        }

        for (int i=3; i<nf; i+=2){//every set will contain 1.
            double[] fAi = fA[i];
            double[] fBi = fB[i];
            System.arraycopy(fC[i], 0, fBi, 0, npairs+1);
            for (int if1=0; if1<=npairs; if1++){
                fAi[if1] = 0;
            }
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            int rest = i ^ jBits;
            // sum over proper subsets j of i that contain jBits; jComp = i & ~j
            for (int s=0; s!=rest; s=(s-rest)&rest){
                double[] fBj = fB[jBits|s];
                double[] fCjComp = fC[(rest^s)|1];
                for (int if1=0; if1<=npairs; if1++){
                    double t = 0;
                    for (int k=0; k<=if1; k++){
                        t += fBj[k] * fCjComp[if1-k];
                    }
                    fAi[if1] += t;
                }
            }
            for (int if1=0; if1<=npairs; if1++){
                fBi[if1] -= fAi[if1];//remove from B graphs that contain articulation point 0.
            }
        }

        for (int v=1; v<n; v++){
            int vs1 = 1<<v;
            for (int i=vs1+1; i<nf; i++){
                double[] fAi = fA[i];
                for (int if1=0; if1<=npairs; if1++){
                    fAi[if1] = 0;
                }
                if ( (i & vs1) == 0 ) continue;
                int iLowBit = (i & -i);
                if ( iLowBit == i ) continue;

                // jBits is v plus the lowest bit of i other than v
                int jBits;
                if (iLowBit == vs1) {
                    int ii = i ^ iLowBit;
                    jBits = vs1 | (ii & -ii);
                }
                else {
                    jBits = iLowBit | vs1;
                }
                int rest = i ^ jBits;
                // sum over proper subsets j of i that contain jBits; jComp = i & ~j
                for (int s=0; s!=rest; s=(s-rest)&rest){
                    double[] fBj = fB[jBits|s];
                    int jCompv = (rest^s) | vs1;
                    double[] fBjComp = fB[jCompv];
                    double[] fAjComp = fA[jCompv];
                    for (int if1=0; if1<=npairs; if1++){
                        fBA[if1] = fBjComp[if1] + fAjComp[if1];
                    }
                    for (int if1=0; if1<=npairs; if1++){
                        double t = 0;
                        for (int k=0; k<=if1; k++){
                            t += fBj[k] * fBA[if1-k];
                        }
                        fAi[if1] += t;
                    }
                }
                double[] fBi = fB[i];
                for (int if1=0; if1<=npairs; if1++){
                    fBi[if1] -= fAi[if1];//remove from B graphs
                }
            }
        }
        // we're actually messing up fB here, but we don't need it anymore
        double[] rv = fB[nf-1];
        double fac = (1.0-n)/SpecialFunctions.factorial(n);
        for (int i=0; i<rv.length; i++) {
            rv[i] *= fac;
        }
        return rv;
    }

	/**
	 * Sets the temperature passed to the Mayer functions. Because the temperature
	 * dependence is carried by x, the Mayer functions should not use it.
	 *
	 * @param temperature the temperature, in simulation energy units (kT)
	 */
	public void setTemperature(double temperature) {
		beta = 1/temperature;
	}

	/**
	 * {@return the number of core edges (overlapping pairs) found by the most recent
	 * {@link #checkConfig(BoxCluster)}}
	 */
	public int getCoreEdgeCount() {
		return edgeCountCore;
	}

	/**
	 * {@return the number of well edges (pairs inside the well) found by the most recent
	 * {@link #checkConfig(BoxCluster)}}
	 */
	public int getWellEdgeCount() {
		return edgeCountWell;
	}

	/** Numbers of core and well edges from the most recent {@link #checkConfig(BoxCluster)}. */
	protected int edgeCountCore, edgeCountWell;
	/** Number of core and well edges at each point, from the most recent {@link #checkConfig(BoxCluster)}. */
	protected final int[] outDegreeCore, outDegreeWell;

	/**
	 * {@return the number of core edges at each point, from the most recent
	 * {@link #checkConfig(BoxCluster)}} This is an internal array.
	 */
	public int[] getOutDegreeCore() {
		return outDegreeCore;
	}

	/**
	 * {@return the number of well edges at each point, from the most recent
	 * {@link #checkConfig(BoxCluster)}} This is an internal array.
	 */
	public int[] getOutDegreeWell() {
		return outDegreeWell;
	}

	/**
	 * {@return the number of cliques of well edges (sets of three or more points, not the
	 * full set, in which every pair is a well edge) found by the most recent
	 * {@link #checkConfig(BoxCluster)}}
	 */
	public int getF1CliqueCount() {
	    return cliqueCountF1;
	}

	/**
	 * {@return the number of cliques of core edges (sets of three or more points, not the
	 * full set, in which every pair overlaps) found by the most recent
	 * {@link #checkConfig(BoxCluster)}}
	 */
	public int getE2CliqueCount() {
        return cliqueCountE2;
    }

    /**
     * {@return the number of cliques of edges of either kind (sets of three or more
     * points, not the full set, in which every pair is a core or well edge) found by the
     * most recent {@link #checkConfig(BoxCluster)}}
     */
    public int getEFCliqueCount() {
        return cliqueCountEF;
    }

    /**
     * {@return the number of sets of three or more points (not the full set) in which no
     * pair is an edge, found by the most recent {@link #checkConfig(BoxCluster)}}
     */
    public int getNoneCliqueCount() {
        return cliqueCountNone;
    }

    /**
     * {@return the well-edge cliques as subset bitmasks} Only the first
     * {@link #getF1CliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getF1Cliques() {
        return cliqueListF1;
    }

    /**
     * {@return the core-edge cliques as subset bitmasks} Only the first
     * {@link #getE2CliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getE2Cliques() {
        return cliqueListE2;
    }

    /**
     * {@return the cliques of edges of either kind as subset bitmasks} Only the first
     * {@link #getEFCliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getEFCliques() {
        return cliqueListEF;
    }

    /**
     * {@return the edge-free sets as subset bitmasks} Only the first
     * {@link #getNoneCliqueCount()} entries are valid. This is an internal array.
     */
    public int[] getNoneCliques() {
        return cliqueListNone;
    }

	/** Clique counts from the most recent {@link #checkConfig(BoxCluster)}. */
	protected int cliqueCountF1, cliqueCountE2, cliqueCountEF, cliqueCountNone;
	/** Whether each subset is a clique of each kind, indexed by subset bitmask. */
    protected final boolean[] cliqueSetF1, cliqueSetE2, cliqueSetEF, cliqueSetNone;
    /** Lists of cliques of each kind, as subset bitmasks. */
    protected final int[] cliqueListF1, cliqueListE2, cliqueListEF, cliqueListNone;

	/**
	 * Evaluates the Mayer functions (via {@link #updateF(BoxCluster)}), builds the graph
	 * of core and well edges, and tests it for biconnectivity. If it is biconnected,
	 * also records the edge counts, degrees and clique lists returned by the getters.
	 *
	 * @param box the box holding the configuration
	 * @return false if every coefficient is zero for this configuration (the edge graph
	 *         is not biconnected); true otherwise
	 */
	public boolean checkConfig(BoxCluster box) {
	    updateF(box);
	    edgeCountCore = edgeCountWell = 0;
	    outDegreeWell[0]=outDegreeCore[0]=0;
	    for(int i=0; i<n; i++){
        	for(int j=i+1; j<n; j++){
        	    if (i==0) {
        	        outDegreeWell[j]=outDegreeCore[j]=0;
        	    }
        	    int k = (1<<i)|(1<<j);
        	    double e2 = fQ[k][0];
        		double f1 = fQ[k][1];
        		if (e2==0 || f1>0){//Means when there is an edge.
        			g.putEdge((byte)i, (byte)j);
        			if (e2==1) {
        			    //g.getEdge((byte)i, (byte)j).setColor('B');
                        cliqueSetF1[k] = true;
                        cliqueSetE2[k] = false;
        				edgeCountWell++;
        				outDegreeWell[i]++;
        				outDegreeWell[j]++;
        			}
        			else {
                        cliqueSetF1[k] = false;
                        cliqueSetE2[k] = true;
        				edgeCountCore++;
        				outDegreeCore[i]++;
        				outDegreeCore[j]++;
        			}
        		}
        		else {
                    cliqueSetF1[k] = cliqueSetE2[k] = false;
                    if(g.hasEdge((byte)i, (byte)j)){
                        g.deleteEdge((byte)i, (byte)j);//This is remove edges set by last step. It is a initilization.
                    }
        		}
        		cliqueSetEF[k] = (cliqueSetF1[k] || cliqueSetE2[k]);
        		cliqueSetNone[k] = !cliqueSetEF[k];
        	}
	    }
	    if (!isB.check(g)) return false;
        // the existence of a clique separator indicates that the value for
        // this configuration is zero.  Loop through all sets, considering
        // each as a clique separator.
        cliqueCountF1 = cliqueCountE2 = cliqueCountEF = cliqueCountNone = 0;
        for (int i=3; i<nf-1; i++) {
            int j = i & -i;//lowest bit in i
            if (i==j) {
                // 1-point set.
                continue;
            }
            int k = i&~j; //strip j bit from i and set result to k
            int jj = k & -k; // 2nd lowest bit
            if (k == jj) {
                // 2-point set.  cliqueSet[i] was set in the above loop
                // over pairs.
                continue;
            }

            int kk = i&~jj; // strip jj bit from i

            cliqueSetF1[i] = cliqueSetF1[k] && cliqueSetF1[kk] && cliqueSetF1[j|jj];
            cliqueSetE2[i] = cliqueSetE2[k] && cliqueSetE2[kk] && cliqueSetE2[j|jj];
            cliqueSetEF[i] = cliqueSetEF[k] && cliqueSetEF[kk] && cliqueSetEF[j|jj];
            cliqueSetNone[i] = cliqueSetNone[k] && cliqueSetNone[kk] && cliqueSetNone[j|jj];

            if (cliqueSetF1[i]) {
                // i is a clique
                cliqueListF1[cliqueCountF1] = i;
                cliqueCountF1++;
            }
            else if (cliqueSetE2[i]) {
                // i is a clique
                cliqueListE2[cliqueCountE2] = i;
                cliqueCountE2++;
            }
            if (cliqueSetEF[i]) {
                // i is a clique
                cliqueListEF[cliqueCountEF] = i;
                cliqueCountEF++;
            }
            else if (cliqueSetNone[i]) {
                // i is a clique
                cliqueListNone[cliqueCountNone] = i;
                cliqueCountNone++;
            }
        }

	    return true;
	}
}
