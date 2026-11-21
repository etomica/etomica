/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.wheatley;

import etomica.math.SpecialFunctions;
import etomica.virial.AtomPairSet;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.ClusterAbstract;
import etomica.virial.cluster.ClusterAbstractMultivalue;

/**
 * Computes the integrand of the <i>n</i>-th virial coefficient together with its
 * perturbation-theory coefficients, using Wheatley's recursion for the sum over
 * all biconnected Mayer diagrams.
 *
 * <h2>What is computed</h2>
 * The pair potential is split into a reference part and a perturbation,
 * <i>u</i>(<i>r</i>) = <i>u</i><sub>R</sub>(<i>r</i>) + &lambda; <i>u</i><sub>P</sub>(<i>r</i>),
 * so that each Mayer e-bond is
 * <i>e</i>(&lambda;) = <i>e</i><sub>R</sub> exp(&minus;&lambda;&beta;<i>u</i><sub>P</sub>).
 * For one configuration of the <i>n</i> points, this cluster returns the Taylor
 * coefficients in &lambda; of the biconnected-diagram sum,
 * <blockquote>
 *   {@code value[m]} = (1 &minus; <i>n</i>)/<i>n</i>! &times;
 *   [&lambda;<sup>m</sup>] <i>f</i><sub>B</sub>(&lambda;),
 *   &nbsp;&nbsp; <i>m</i> = 0, 1, &hellip;, {@code nOrder},
 * </blockquote>
 * where [&lambda;<sup>m</sup>] denotes the coefficient of &lambda;<sup>m</sup> and
 * <i>f</i><sub>B</sub> is the sum over biconnected diagrams on all <i>n</i> points.
 * {@code value[0]} is the integrand for <i>B</i><sub>n</sub> of the reference potential
 * alone, and {@code value[m]} is the integrand for the <i>m</i>-th order perturbation
 * correction. At &lambda; = 1, the sum of the coefficients approaches the integrand
 * for the full potential as {@code nOrder} increases, provided the series converges.
 * Each coefficient is exact for the given configuration (up to round-off); the
 * expansion is truncated only at {@code nOrder}. Integration over
 * configurations is done by Mayer-sampling Monte Carlo (for example
 * {@link etomica.virial.simulations.perturbative.VirialLJPT}, which splits the
 * Lennard-Jones potential into a WCA reference and its attractive remainder).
 *
 * <h2>Algorithm</h2>
 * For each subset <i>S</i> of the points, the product of e-bonds expands as
 * <i>f</i><sub>Q</sub>(<i>S</i>) = <i>e</i><sub>R</sub>(<i>S</i>)
 * exp(&minus;&lambda;&beta;<i>U</i><sub>P</sub>(<i>S</i>)), whose &lambda;<sup>m</sup>
 * coefficient is <i>e</i><sub>R</sub>(<i>S</i>) (&minus;&beta;<i>U</i><sub>P</sub>(<i>S</i>))<sup>m</sup>/<i>m</i>!,
 * where <i>e</i><sub>R</sub>(<i>S</i>) is the product of reference e-bonds and
 * <i>U</i><sub>P</sub>(<i>S</i>) the total perturbation energy within <i>S</i>.
 * These are propagated through the recursion
 * <i>f</i><sub>Q</sub> &rarr; <i>f</i><sub>C</sub> (connected) &rarr;
 * <i>f</i><sub>B</sub> (biconnected) of
 * R. J. Wheatley, Phys. Rev. Lett. <b>110</b>, 200601 (2013).
 * Every step of the recursion is a sum of products of two subset quantities, so each
 * product of power series becomes a Cauchy product (convolution) of their
 * coefficients, truncated at {@code nOrder}. Working with Taylor coefficients rather
 * than &lambda;-derivatives means no binomial factors are needed (compare
 * {@link ClusterWheatleySoftDerivatives}, which carries &beta;-derivatives). The cost
 * is O(3<sup><i>n</i></sup> {@code nOrder}<sup>2</sup>) per configuration, and memory
 * is O(2<sup><i>n</i></sup> {@code nOrder}).
 *
 * <h2>Inputs and conventions</h2>
 * <ul>
 * <li>{@code fR} must return the reference Mayer function
 *     <i>f</i><sub>R</sub> = <i>e</i><sub>R</sub> &minus; 1.</li>
 * <li>{@code fP} must return exp(&minus;&beta;<i>u</i><sub>P</sub>) &minus; 1, that is,
 *     a Mayer function built from the perturbation potential alone. The class
 *     recovers &beta;<i>u</i><sub>P</sub> from it as &minus;ln(<i>f</i><sub>P</sub> + 1).
 *     <i>u</i><sub>P</sub> must be finite wherever <i>e</i><sub>R</sub> is non-zero.</li>
 * <li>Pairs with <i>e</i><sub>R</sub> = 0 (inside a hard core) are handled exactly:
 *     every coefficient for a set containing that pair is zero.</li>
 * <li>Both Mayer functions are evaluated with the &beta; set by
 *     {@link #setTemperature(double)}.</li>
 * </ul>
 *
 * <h2>Limitations</h2>
 * <ul>
 * <li>There is no BigDecimal fallback. The {@code tol} constructor argument is
 *     stored but not used, so values whose magnitude is near round-off (for example
 *     where many diagrams cancel at large separations) are returned as computed in
 *     double precision. {@link ClusterWheatleySoftDerivatives} shows how a
 *     BigDecimal check can be added.</li>
 * <li>The Mayer functions' {@code setBox} methods are not called, so they should not
 *     depend on the box (spherical pair potentials are fine).</li>
 * <li>Values are cached by the configuration ID of the box's
 *     {@link CoordinatePairSet}, including the configuration before the most recent
 *     one, so that a rejected Monte Carlo trial does not trigger a recomputation.</li>
 * </ul>
 *
 * <p>Individual orders can be used as ordinary clusters through
 * {@link #makeClusterOrder(int)}, for example as extra target clusters in an overlap
 * sampling simulation.
 *
 * @see ClusterWheatleySoftDerivatives
 * @see ClusterWheatleyExtendSW
 */
public class ClusterWheatleyPT implements ClusterAbstractMultivalue {

    /** Configuration ID of the box when the current values were computed. */
    protected long cPairID = -1;
    /** Configuration ID for the values held in {@link #lastValue}. */
    protected long lastCPairID = -1;
    /** Most recently computed coefficients, {@code value[m]} for order m. */
    protected final double[] value;
    /** Coefficients for the configuration before the current one. */
    protected final double[] lastValue;
    /** Number of points (molecules). */
    protected final int n;
    /** Number of subsets of the points, 2<sup>n</sup>. */
    protected final int nf;
    /** Highest perturbation order computed. */
    protected final int nOrder;
    /** Reference Mayer function, e<sub>R</sub> &minus; 1. */
    protected final MayerFunction fR;
    /** Mayer function of the perturbation potential alone, exp(&minus;&beta;u<sub>P</sub>) &minus; 1. */
    protected final MayerFunction fP;
    /**
     * Taylor coefficients of the product of e-bonds for each subset, indexed
     * [subset bitmask][order]. Pair entries are set by {@link #updateF(BoxCluster)};
     * entries for larger subsets are set by {@link #calcFullFQ()}.
     */
    protected final double[][] fQ;
    /**
     * Total reduced perturbation energy &beta;U<sub>P</sub> within each subset, indexed
     * by subset bitmask. Positive infinity marks a pair inside the reference hard core.
     */
    protected final double[] fU;
    /** Taylor coefficients of the connected-diagram sum f<sub>C</sub> for each subset. */
    protected final double[][] fC;
    /** Work array for the articulation-point recursion, indexed like {@link #fB}. */
    protected final double[][] fA;
    /**
     * Taylor coefficients of the biconnected-diagram sum f<sub>B</sub> for each subset.
     * After {@link #calcValue(BoxCluster)}, {@code fB[nf-1]} holds the full-set result.
     */
    protected final double[][] fB;
    /** Work array holding f<sub>B</sub> + f<sub>A</sub> for one subset during the recursion. */
    protected final double[] fBA;
    /** Reciprocal temperature, 1/kT. */
    protected double beta;
    /** Stored but currently unused; see the class description. */
    protected double tol;
    /** If true (the default), values are cached by configuration ID. */
    protected boolean doCaching = true;

    /**
     * Creates the cluster. Call {@link #setTemperature(double)} before use.
     *
     * @param nPoints number of points (molecules), n; each work array has 2<sup>n</sup>
     *                rows, and the cost grows as 3<sup>n</sup>
     * @param nOrder  highest order in &lambda; to compute; {@link #getNumValues()}
     *                returns nOrder + 1
     * @param fR      reference Mayer function, e<sub>R</sub> &minus; 1
     * @param fP      Mayer function of the perturbation potential alone,
     *                exp(&minus;&beta;u<sub>P</sub>) &minus; 1
     * @param tol     currently unused; see the class description
     */
    public ClusterWheatleyPT(int nPoints, int nOrder, MayerFunction fR, MayerFunction fP, double tol) {
        this.n = nPoints;
        this.nf = 1 << n; //2^n
        this.tol = tol;
        this.nOrder = nOrder;
        this.fR = fR;
        this.fP = fP;
        value = new double[nOrder + 1];
        lastValue = new double[nOrder + 1];
        fQ = new double[nf][nOrder + 1];
        for (int i = 0; i < n; i++) {
            fQ[1 << i][0] = 1.0;
            for (int if1 = 1; if1 <= nOrder; if1++) {
                fQ[1 << i][if1] = 0.0;
            }
        }
        fU = new double[nf];
        fC = new double[nf][nOrder + 1];
        fA = new double[nf][nOrder + 1];
        fB = new double[nf][nOrder + 1];
        fBA = new double[nOrder + 1];
    }

    /**
     * Returns a new cluster with the same points, order, Mayer functions and
     * temperature. Cached values are not copied.
     */
    public ClusterAbstract makeCopy() {
        ClusterWheatleyPT c = new ClusterWheatleyPT(n, nOrder, fR, fP, tol);
        c.setTemperature(1 / beta);
        return c;
    }

    /** Returns the number of points, n. */
    public int pointCount() {
        return n;
    }

    /**
     * Returns the zeroth-order coefficient, {@code value[0]}: the integrand for the
     * reference potential. All orders are computed and stored at the same time; use
     * {@link #getAllLastValues(BoxCluster)} or {@link #makeClusterOrder(int)} to get
     * the others.
     *
     * @param box the box holding the configuration
     * @return the zeroth-order coefficient
     * @throws RuntimeException if the result is NaN or infinite
     */
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
                System.arraycopy(lastValue, 0, value, 0, value.length);
                return value[0];
            }

            // a new cluster
            lastCPairID = cPairID;
            System.arraycopy(value, 0, lastValue, 0, value.length);
            cPairID = thisCPairID;
        }
        updateF(box);

        calcValue(box);
        if (Double.isNaN(value[0]) || Double.isInfinite(value[0])) {
            updateF(box);
            calcValue(box);
            throw new RuntimeException("oops " + value[0]);
        }
        return value[0];
    }

    /**
     * Evaluates both Mayer functions for every pair and fills the pair entries of
     * {@link #fQ} with e<sub>R</sub> (&minus;&beta;u<sub>P</sub>)<sup>m</sup>/m!, and of
     * {@link #fU} with &beta;u<sub>P</sub>.
     *
     * @param box the box holding the configuration
     */
    protected void updateF(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();

        // recalculate all f values for all pairs
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int ij = 1 << i | 1 << j;
                fQ[ij][0] = fR.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta) + 1;
                if (fQ[ij][0] == 0) {
                    fU[ij] = Double.POSITIVE_INFINITY;
                    for (int k = 1; k <= nOrder; k++) {
                        fQ[ij][k] = 0;
                    }
                    continue;
                }
                double fp = fP.f(aPairs.getAPair(i, j), cPairs.getr2(i, j), beta);

                double buij = -Math.log(fp + 1);
                fU[ij] = buij;
                for (int m = 1; m <= nOrder; m++) {
                    fQ[ij][m] = fQ[ij][m - 1] * (-buij) / m;
                }
            }
        }
    }

    /**
     * Fills {@link #fQ} and {@link #fU} for every subset of three or more points,
     * given that the pair entries have already been set by {@link #updateF(BoxCluster)}.
     * For subset S, fQ[S][m] = e<sub>R</sub>(S) (&minus;&beta;U<sub>P</sub>(S))<sup>m</sup>/m!.
     */
    protected void calcFullFQ() {
        for (int i = 3; i < nf; i++) {
            int j = i & -i;
            if (i == j) continue; // 1-point set
            int k = i & ~j;
            if (k == (k & -k)) continue; // 2-point set

            fQ[i][0] = fQ[k][0];
            fU[i] = fU[k];

            //loop over pairs formed from j and each point in set i; multiply by bond for each pair
            //all such pairs will be with bits higher than j, as j is the lowest bit in i
            for (int l = (j << 1); l < i; l = (l << 1)) {
                if ((l & i) == 0) continue; //l is not in i
                fQ[i][0] = fQ[i][0] * fQ[l | j][0];
                fU[i] += fU[l | j];
            }
            if (fQ[i][0] == 0) {
                for (int m = 1; m <= nOrder; m++) {
                    fQ[i][m] = 0;
                }
                continue;
            }
            for (int m = 1; m <= nOrder; m++) {
                fQ[i][m] = fQ[i][m - 1] * (-fU[i] / m);
            }
        }
//		System.out.println("fQ="); MyUtility.display2DArray(fQ);
    }

    /**
     * Runs the recursion for the current configuration and stores the coefficients in
     * {@link #value}. {@link #updateF(BoxCluster)} must have been called first for this
     * configuration. Most callers should use {@link #value(BoxCluster)} or
     * {@link #getAllLastValues(BoxCluster)}, which handle caching and call
     * {@code updateF}.
     *
     * <p>Each sum over subsets visits only the proper subsets j of a set i that contain
     * a fixed set of bits, using the submask step {@code s = (s - rest) & rest}.
     *
     * @param box the box holding the configuration (not used directly)
     */
    public void calcValue(BoxCluster box) {
        calcFullFQ();
        //Compute the fC's
        for (int i = 1; i < nf; i++) {
            double[] fCi = fC[i];
            System.arraycopy(fQ[i], 0, fCi, 0, nOrder + 1);
            int iLowBit = i & -i;
            int rest = i ^ iLowBit;
            // sum over proper subsets j of i that contain iLowBit; jComp = i & ~j
            for (int s = 0; s != rest; s = (s - rest) & rest) {
                double[] fCj = fC[iLowBit | s];
                double[] fQjComp = fQ[rest ^ s];
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    double t = 0;
                    for (int k = 0; k <= if1; k++) {
                        t += fCj[k] * fQjComp[if1 - k];
                    }
                    fCi[if1] -= t;
                }
            }
        }

        //find fA1
        for (int i = 2; i < nf; i += 2) {//all even sets (2,4,6) don't contain 1
            System.arraycopy(fC[i], 0, fB[i], 0, nOrder + 1);
        }

        for (int if1 = 0; if1 <= nOrder; if1++) {
            fA[1][if1] = 0;
            fB[1][if1] = fC[1][if1];
        }

        for (int i = 3; i < nf; i += 2) {//every set will contain 1.
            double[] fAi = fA[i];
            double[] fBi = fB[i];
            System.arraycopy(fC[i], 0, fBi, 0, nOrder + 1);
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fAi[if1] = 0;
            }
            int ii = i - 1;//all bits in i but lowest
            int iLow2Bit = (ii & -ii);//next lowest bit
            int jBits = 1 | iLow2Bit;
            int rest = i ^ jBits;
            // sum over proper subsets j of i that contain jBits; jComp = i & ~j
            for (int s = 0; s != rest; s = (s - rest) & rest) {
                double[] fBj = fB[jBits | s];
                double[] fCjComp = fC[(rest ^ s) | 1];
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    double t = 0;
                    for (int k = 0; k <= if1; k++) {
                        t += fBj[k] * fCjComp[if1 - k];
                    }
                    fAi[if1] += t;
                }
            }
            for (int if1 = 0; if1 <= nOrder; if1++) {
                fBi[if1] -= fAi[if1];//remove from B graphs that contain articulation point 0.
            }
        }

        for (int v = 1; v < n; v++) {
            int vs1 = 1 << v;
            for (int i = vs1 + 1; i < nf; i++) {
                double[] fAi = fA[i];
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    fAi[if1] = 0;
                }
                if ((i & vs1) == 0) continue;
                int iLowBit = (i & -i);
                if (iLowBit == i) continue;

                // jBits is v plus the lowest bit of i other than v
                int jBits;
                if (iLowBit == vs1) {
                    int ii = i ^ iLowBit;
                    jBits = vs1 | (ii & -ii);
                } else {
                    jBits = iLowBit | vs1;
                }
                int rest = i ^ jBits;
                // sum over proper subsets j of i that contain jBits; jComp = i & ~j
                for (int s = 0; s != rest; s = (s - rest) & rest) {
                    double[] fBj = fB[jBits | s];
                    int jCompv = (rest ^ s) | vs1;
                    double[] fBjComp = fB[jCompv];
                    double[] fAjComp = fA[jCompv];
                    for (int if1 = 0; if1 <= nOrder; if1++) {
                        fBA[if1] = fBjComp[if1] + fAjComp[if1];
                    }
                    for (int if1 = 0; if1 <= nOrder; if1++) {
                        double t = 0;
                        for (int k = 0; k <= if1; k++) {
                            t += fBj[k] * fBA[if1 - k];
                        }
                        fAi[if1] += t;
                    }
                }
                double[] fBi = fB[i];
                for (int if1 = 0; if1 <= nOrder; if1++) {
                    fBi[if1] -= fAi[if1];//remove from B graphs
                }
            }
        }

        double fac = (1.0 - n) / SpecialFunctions.factorial(n);
        for (int i = 0; i < value.length; i++) {
            value[i] = fB[nf - 1][i] * fac;
        }
    }

    /**
     * Sets the temperature used to evaluate both Mayer functions.
     *
     * @param temperature the temperature, in simulation energy units (kT)
     */
    public void setTemperature(double temperature) {
        beta = 1 / temperature;
    }

    /**
     * Returns a cluster whose value is the order-o coefficient of this cluster.
     * The returned cluster does not compute anything itself: it reports
     * {@code value[o]} from this cluster's most recent evaluation, so this cluster must
     * be evaluated for the same configuration first (as it is when it is the sampled
     * target cluster).
     *
     * @param o the order, 0 &le; o &le; nOrder
     * @return a cluster reporting the order-o coefficient
     */
    public ClusterRetrieveOrders makeClusterOrder(int o) {
        return new ClusterRetrieveOrders(o);
    }

    /** Returns nOrder + 1, the number of coefficients computed. */
    @Override
    public int getNumValues() {
        return nOrder + 1;
    }

    /**
     * Computes (or retrieves from the cache) all coefficients for the box's current
     * configuration.
     *
     * @param box the box holding the configuration
     * @return the coefficients, indexed by order; this is the internal array, so copy
     *         it if it must survive the next evaluation
     */
    @Override
    public double[] getAllLastValues(BoxCluster box) {
        value(box);
        return value;
    }

    /**
     * A view of one perturbation order of the enclosing {@link ClusterWheatleyPT}.
     * Created by {@link ClusterWheatleyPT#makeClusterOrder(int)}; see there for the
     * requirement that the enclosing cluster be evaluated first.
     */
    public class ClusterRetrieveOrders implements ClusterAbstract {
        /** The order reported by this cluster. */
        protected final int iOrder;

        /**
         * @param iOrder the order to report, 0 &le; iOrder &le; nOrder
         */
        public ClusterRetrieveOrders(int iOrder) {
            this.iOrder = iOrder;
        }

        /** Not supported; returns null. */
        public ClusterAbstract makeCopy() {
            return null;
        }

        /** Returns the number of points of the enclosing cluster. */
        public int pointCount() {
            return n;
        }

        /**
         * Returns the order-iOrder coefficient from the enclosing cluster's most recent
         * evaluation. The box is not examined.
         */
        public double value(BoxCluster box) {
            return value[iOrder];
        }

        /** Does nothing; the temperature is set on the enclosing cluster. */
        public void setTemperature(double temperature) {
        }
    }

    /**
     * Reports the order-n derivative from a {@link ClusterWheatleySoftDerivatives}
     * cluster when its most recent value was computed with BigDecimal arithmetic, and
     * 0 otherwise.
     *
     * <p>Note: this class wraps a {@code ClusterWheatleySoftDerivatives}, not a
     * {@code ClusterWheatleyPT}, and is not used by {@code ClusterWheatleyPT}. It
     * duplicates {@link ClusterWheatleySoftDerivatives.ClusterRetrievePrimesBD}.
     */
    public static class ClusterRetrievePrimesBD implements ClusterAbstract {
        /** The cluster whose values are reported. */
        protected final ClusterWheatleySoftDerivatives cluster;
        /** The derivative order reported. */
        protected final int n;

        /**
         * @param cluster the cluster whose values are reported
         * @param n       the derivative order to report
         */
        public ClusterRetrievePrimesBD(ClusterWheatleySoftDerivatives cluster, int n) {
            this.cluster = cluster;
            this.n = n;
        }

        /** Not supported; returns null. */
        public ClusterAbstract makeCopy() {
            return null;
        }

        /** Returns the number of points of the wrapped cluster. */
        public int pointCount() {
            return cluster.n;
        }

        /**
         * Returns the wrapped cluster's order-n derivative if its most recent value was
         * computed with BigDecimal arithmetic, and 0 otherwise.
         */
        public double value(BoxCluster box) {
            return cluster.valueBD ? cluster.value[n] : 0;
        }

        /** Does nothing. */
        public void setTemperature(double temperature) {
        }
    }
}
