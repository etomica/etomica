/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.*;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.*;
import etomica.graph.operations.AllIsomorphs.AllIsomorphsParameters;
import etomica.graph.operations.BiComponentSubst.BiComponentSubstParameters;
import etomica.graph.operations.ComponentSubst.ComponentSubstParameters;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DifByConstant.DifByConstantParameters;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.ReeHoover.ReeHooverParameters;
import etomica.graph.operations.SplitOneBiconnected.SplitOneParametersBC;
import etomica.graph.property.*;
import etomica.graph.traversal.BCVisitor;
import etomica.graph.traversal.CVisitor;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.util.Arrays;
import etomica.virial.MayerFunction;
import etomica.virial.MayerFunctionNonAdditive;
import etomica.virial.cluster.CondenseExchange.CondenseExchangeParameters;
import etomica.virial.cluster.ExchangeSplit.ExchangeSplitParameters;

import java.util.*;

public class VirialDiagrams {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, disconnectedP;
    protected Set<Graph> minMultiP, fullMultiP, trueMultiP;
    protected Set<Graph> rho;
    protected Set<Graph> lnfXi, fullLnXi;
    protected Map<Graph,Graph> cancelMap;
    protected boolean doShortcut;
    protected boolean doMinimalMulti;
    protected boolean doMultiFromPair;
    protected boolean doMinimalBC;
    protected boolean doKeepEBonds;
    protected boolean doExchange;
    protected boolean doExchangeF;
    protected boolean doExchangeCondensing;
    protected boolean doDisconnectedMatching = true;
    protected boolean doNegativeExchange = false;
    protected boolean doHB;
    protected boolean flexCancelOnly = false; // this boolean makes me sad
    protected final char nodeColor = Metadata.COLOR_CODE_0;
    protected char[] flexColors;
    protected boolean allPermutations = false;
    public char fBond, bBond, eBond, excBond, mBond, mmBond, fmBond, efbcBond, ffBond, mxcBond, MxcBond;

    protected static int[][][] groupStart = new int[0][0][0];
    protected static int[][] myIds = new int[0][0];
    protected static int[][] tripletStart = new int[0][0];
    protected static int[][] quadStart = new int[0][0];
    protected static int[][] quintStart = new int[0][0];
    protected static int[][] sixStart = new int[0][0];

    public static void main(String[] args) {
        final int n = 4;
        boolean multibody = true;
        boolean flex = true;
        boolean doKeepEBonds = false;
        boolean doReeHoover = false;
        boolean doExchange = false;
        VirialDiagrams virialDiagrams = new VirialDiagrams(n, multibody, flex, true);
        virialDiagrams.setDoReeHoover(doReeHoover);
        virialDiagrams.setDoKeepEBonds(doKeepEBonds);
        virialDiagrams.setDoHB(false);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.setDoExchange(doExchange);
        if (doExchange) {
            virialDiagrams.setDoExchangeF(true);
            virialDiagrams.setDoExchangeCondensing(true);
        }
        if (multibody) {
            virialDiagrams.setDoMinimalMulti(true);
            virialDiagrams.setDoMultiFromPair(true);
        }
        if (doReeHoover && (flex || multibody)) {
            virialDiagrams.setDoMinimalBC(true);
        }
        virialDiagrams.makeVirialDiagrams();
    }

    public VirialDiagrams(int n, boolean multibody, boolean flex) {
        this(n, multibody, flex, false);
    }

    public VirialDiagrams(int n, boolean multibody, boolean flex, boolean interactive) {
        this.multibody = multibody;
        this.flex = flex;
        this.n = n;
        this.isInteractive = interactive;
        doReeHoover = true;
        doShortcut = false;
        init();
    }

    public void setDoReeHoover(boolean newDoReeHoover) {
        if (doKeepEBonds && newDoReeHoover) {
            throw new RuntimeException("can't have both Ree-Hoover and e-bond representation");
        }
        doReeHoover = newDoReeHoover;
    }

    public void setDoShortcut(boolean newDoShortcut) {
        doShortcut = newDoShortcut;
    }

    public void setDoKeepEBonds(boolean newDoKeepEBonds) {
        doKeepEBonds = newDoKeepEBonds;
        if (!flex && doKeepEBonds) {
            if (multibody) {
                throw new RuntimeException("keep-e-bonds doesn't work with multibody unless flex is on");
            }
            System.out.println("keep-e-bonds doesn't behave well with flex off");
        }
    }

    /**
     * Enable construction of diagrams from the Hellmann-Bich formulation
     */
    public void setDoHB(boolean newDoHB) {
        doHB = newDoHB;
        if (!doKeepEBonds && doHB) {
            throw new RuntimeException("Hellmann-Bich wants e-bonds");
        }
    }

    public void setDoExchange(boolean newDoExchange) {
        if (!doKeepEBonds && newDoExchange) {
            throw new RuntimeException("exchange only works with keep-eBonds");
        }
        doExchange = newDoExchange;
    }

    public void setDoDisconnectedMatching(boolean newDoDisconnectedMatching) {
        doDisconnectedMatching = newDoDisconnectedMatching;
    }

    /**
     * Turn on replacement of eBonds with bonds that represent a single fBond
     * between an exchange group and another molecule (or another exchange
     * group).
     */
    public void setDoExchangeF(boolean newDoExchangeF) {
        if (!doExchange) {
            throw new RuntimeException("this only makes sense with exchange on");
        }
        doExchangeF = newDoExchangeF;
    }
    
    /**
     * Turn on replacement of exchanged groups of nodes with a single node (of
     * a different color)
     */
    public void setDoExchangeCondensing(boolean newDoExchangeCondensing) {
        if (!doExchange) {
            throw new RuntimeException("this only makes sense with exchange on");
        }
        doExchangeCondensing = newDoExchangeCondensing;
        if (doExchangeCondensing && Metadata.COLORS.size() < 15) {
            // need more colors for points
            Metadata.COLORS.add("brown");
            Metadata.COLORS.add("greenyellow");
            Metadata.COLORS.add("powderblue");
            Metadata.COLORS.add("thistle");
            Metadata.COLORS.add("lightslategray");
            Metadata.COLORS.add("coral");
        }
    }
    
    public void setDoNegativeExchange(boolean newDoNegativeExchange) {
        if (!doExchange) {
            throw new RuntimeException("this only makes sense with exchange on");
        }
        doNegativeExchange = newDoNegativeExchange;
    }

    /**
     * Only include graphs with nonadditive bonds.
     */
    public void setDoMinimalMulti(boolean newDoMinimalMulti) {
        if (!multibody) {
            throw new RuntimeException("can't set minimalMulti without multi");
        }
        doMinimalMulti = newDoMinimalMulti;
    }

    public void setDoMultiFromPair(boolean newDoMultiFromPair) {
        if (!multibody) {
            throw new RuntimeException("can't set multi from pair without multi");
        }
        doMultiFromPair = newDoMultiFromPair;
    }

    public void setDoMinimalBC(boolean newDoMinimalBC) {
        if (!multibody && !flex) {
            throw new RuntimeException("you need multi or flex to do minimal bc");
        }
        doMinimalBC = newDoMinimalBC;
    }
    
    protected void init() {
        fBond = 'f';
        bBond = 'B';
        eBond = 'e';
        mBond = 'm';  // multi-body
        fmBond = 'z';  // full-multi-body (diagram represents the full virial coefficient)
        mmBond = 'M';  // Multi-body
        mxcBond = 'n';  // multi-body exchange (yes, of course 'n' is bad)
        MxcBond = 'N';
        efbcBond = 'b';
        excBond = 'x';
        ffBond = 'F';
    }

    public void setAllPermutations(boolean newAllPermutations) {
        if (flex && newAllPermutations) {
            System.err.println("all additive permutations and flex don't mix!");
        }
        allPermutations = newAllPermutations;
    }

    public Set<Graph> getVirialGraphs() {
        if (p == null) {
            makeVirialDiagrams();
        }
        return p;
    }

    /**
     * This method only effects the graphs returned by getMSMCGraphs
     */
    public void setFlexCancelOnly(boolean newFlexCancelOnly) {
        if (!flex) {
            throw new RuntimeException("this only makes sense with flex on");
        }
        flexCancelOnly = newFlexCancelOnly;
    }

    public Set<Graph> getMSMCGraphs(boolean connectedOnly, boolean getMultiGraphs) {
        if (p == null) {
            makeVirialDiagrams();
        }
        GraphList allP = makeGraphList();
        if (getMultiGraphs) {
            if (!multibody) throw new RuntimeException("oops");
            if (doMinimalMulti) {
                if (doMultiFromPair) {
                    for (Graph g : p) {
                        if (graphHasEdgeColor(g, mmBond)) {
                            if (flexCancelOnly && cancelMap.get(g) == null) continue;
                            allP.add(g);
                        }
                    }
                }
                else {
                    // we already constructed the exact set we want
                    allP.addAll(minMultiP);
                }
            }
            else {
                // grab the graphs with an mBond
                for (Graph g : p) {
                    if (graphHasEdgeColor(g, mBond)) {
                        allP.add(g);
                    }
                }
            }
        }
        else {
            for (Graph g : p) {
                if (!multibody || (doMinimalMulti && !graphHasEdgeColor(g, mmBond))
                               || (!doMinimalMulti && !graphHasEdgeColor(g, mBond))) {
                    if (flex && flexCancelOnly && cancelMap.get(g) == null) continue;
                    allP.add(g);
                    if (!connectedOnly && cancelMap != null) {
                        Graph c = cancelMap.get(g);
                        if (c != null) {
                            allP.add(c);
                        }
                    }
                }
            }
        }
        GraphList pn = makeGraphList();
        for (Graph g : allP) {
            int fieldCount = 0;
            for (Node node : g.nodes()) {
                if (node.getType() == 'F') {
                    fieldCount++;
                }
            }
            if (fieldCount == n) {
                pn.add(g);
            }
        }
        return pn;
    }
    
    public Map<Graph,Graph> getCancelMap() {
        if (p == null) {
            makeVirialDiagrams();
        }
        return cancelMap;
    }

    public ClusterSum makeVirialCluster(MayerFunction f) {
        return makeVirialCluster(f, null, false);
    }
    public ClusterSum makeVirialCluster(MayerFunction f, MayerFunctionNonAdditive fMulti) {
        return makeVirialCluster(f, fMulti, false);
    }
    public ClusterSum makeVirialCluster(MayerFunction f, MayerFunctionNonAdditive fMulti, boolean doTotal) {
        if (p == null) {
            makeVirialDiagrams();
        }
        boolean doMulti = fMulti != null;
        if (!multibody && doMulti) {
            throw new RuntimeException("can't make multi-body bonds without multi-body diagrams");
        }
        if (flex && multibody && doTotal) {
            throw new RuntimeException("seems like a bad combination: flex + multibody + doTotal");
        }
        if (!doMulti && doTotal) {
            throw new RuntimeException("you want total but don't have a multi-bond?");
        }
        Set<Graph> pn = getMSMCGraphs(false, doMulti);
        if (doTotal) {
            pn.addAll(getMSMCGraphs(false, !doMulti));
        }
        return makeVirialCluster(pn, f, fMulti);
    }
    
    public ClusterSum makeVirialCluster(Set<Graph> graphs, MayerFunction f, MayerFunctionNonAdditive fMulti) {
        
        boolean doMulti = fMulti != null;
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        for (Graph g : graphs) {
            if (multibody && flex && doMultiFromPair && NumRootNodes.value(g) > 1) continue;
            int nDiagrams = populateEFBonds(g, allBonds, weights, false);
            if (nDiagrams > 0 && flex && (!doMulti || doMultiFromPair)) {
                populateEFBonds(g, allBonds, weights, true);
            }
        }

        double[] w = new double[weights.size()];
        for (int i=0; i<w.length; i++) {
            w[i] = weights.get(i);
        }
        if (!doMulti) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f});
        }
        return new ClusterSumMultibody(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}, new MayerFunctionNonAdditive[]{fMulti});
    }
    public int[][] getFlipPointsforDiagram(String diagram) {
        if (diagram.equals("1")) {
            return new int[][]{{0, 1}};
        } else if (diagram.equals("5c")) {
            return new int[][]{{0, 1}, {0, 2}};

        } else if (diagram.equals("54c")) {
            return new int[][]{{0, 3}};
        } else if (diagram.equals("52c")) {
            return new int[][]{{0, 3}, {1, 2}, {0, 1, 2}};
        } else if (diagram.equals("38c")) {
            return new int[][]{{0, 1}, {0, 3}, {0, 2}};
        } else if (diagram.equals("954c")) {
            return new int[][]{{0, 4}};
        } else if (diagram.equals("946c")) {
            return new int[][]{{0, 4}};

        } else if (diagram.equals("952c")) {
            return new int[][]{{0, 4}};

        } else if (diagram.equals("882c")) {
            return new int[][]{};

        } else if (diagram.equals("936c")) {
            return new int[][]{{0, 4}, {0, 1, 2, 3}};

        } else if (diagram.equals("944c")) {
            return new int[][]{{0, 4}, {2, 3}};

        } else if (diagram.equals("930c")) {
            return new int[][]{{0, 4}};

        } else if (diagram.equals("818c")) {
            return new int[][]{{0, 3}, {0, 4}};

        } else if (diagram.equals("928c")) {
            return new int[][]{{0, 4}, {2, 3}, {0, 1, 2, 3}, {1, 2, 3}};

        } else if (diagram.equals("808c")) {
            return new int[][]{{1, 3}, {1, 2}, {0, 4}, {1, 0, 4}};

        } else if (diagram.equals("562c")) {
            return new int[][]{{0, 1}, {0, 2}, {0, 3}, {0, 4}};

        } else {
            throw new RuntimeException("unknown diagram " + diagram);
        }
    }
    public ClusterSum makeVirialCluster(Graph g, MayerFunction f){
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();

        int nDiagrams = populateEFBonds(g, allBonds, weights, false);
        if (nDiagrams > 0 && flex) {
            populateEFBonds(g, allBonds, weights, true);
        }
        if (flex && cancelMap.get(g) != null) {
            Graph cg = cancelMap.get(g);
            populateEFBonds(cg, allBonds, weights, false);
            populateEFBonds(cg, allBonds, weights, true);
        }
        double gCoef = g.coefficient().getValue();


        double[] w = new double[weights.size()];
        for (int i=0; i<w.length; i++) {
            w[i] = weights.get(i)/gCoef;
        }
        return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f});

    }
    public ClusterSum makeVirialClusterTempDeriv(MayerFunction f, MayerFunction e, MayerFunction dfdT) {
        if (p == null) {
            makeVirialDiagrams();
        }
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        Set<Graph> pn = getMSMCGraphs(false, false);
        char dfdTBond = 'a';
        DifByConstantParameters params = new DifByConstantParameters(fBond,dfdTBond);
        DifByConstant dif = new DifByConstant();
        Set<Graph> pn2 = dif.apply(pn, params);
        DifByConstantParameters params2 = new DifByConstantParameters(eBond,dfdTBond);
        Set<Graph> pn3 = dif.apply(pn2, params2);
        for (Graph g : pn3) {
        	populateEFdfdTBonds(g, allBonds, false);
            if (flex) {
            	populateEFdfdTBonds(g, allBonds, true);
            }
            double w = ((double)g.coefficient().getNumerator())/g.coefficient().getDenominator();
            if (flex) {
                w *= 0.5;
            }
            weights.add(w);
            if (flex) {
                weights.add(w);
            }
        }
        double[] w = new double[weights.size()];
        for (int i=0; i<w.length; i++) {
            w[i] = weights.get(i);
        }
        if (n > 3 && !flex && !multibody) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,e,dfdT});
        }
        else if (!multibody) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,dfdT});
        }
        return null;
    }
    
    
    public ClusterSum makeVirialClusterSecondTemperatureDerivative(MayerFunction f, MayerFunction e, MayerFunction dfdT, MayerFunction d2fdT2) {
        if (p == null) {
            makeVirialDiagrams();
        }
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        Set<Graph> pn = getMSMCGraphs(false, false);
        char dfdTBond = 'a';
        char df2dT2Bond = 'b';
        DifByConstantParameters params = new DifByConstantParameters(fBond,dfdTBond);
        DifByConstant dif = new DifByConstant();
        Set<Graph> pn2 = dif.apply(pn, params);
        DifByConstantParameters params2 = new DifByConstantParameters(eBond,dfdTBond);
        Set<Graph> pn3 = dif.apply(pn2, params2);
        
        
        
        DifByConstantParameters params3 = new DifByConstantParameters(fBond,dfdTBond);
        Set<Graph> pn4 = dif.apply(pn3, params3);
        DifByConstantParameters params4 = new DifByConstantParameters(eBond,dfdTBond);
        Set<Graph> pn5 = dif.apply(pn4, params4);
        
        DifByConstantParameters params5 = new DifByConstantParameters(dfdTBond,df2dT2Bond);
        Set<Graph> pn6 = dif.apply(pn3, params5);
        
        pn6.addAll(pn5);       
        
        for (Graph g : pn6) {
        	populateEFdfdTdf2dT2Bonds(g, allBonds, false);
            if (flex) {
            	populateEFdfdTdf2dT2Bonds(g, allBonds, true);
            }
            double w = ((double)g.coefficient().getNumerator())/g.coefficient().getDenominator();
            if (flex) {
                w *= 0.5;
            }
            weights.add(w);
            if (flex) {
                weights.add(w);
            }
        }
        double[] w = new double[weights.size()];
        for (int i=0; i<w.length; i++) {
            w[i] = weights.get(i);
        }
        if (n > 3 && !flex && !multibody) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,e,dfdT,d2fdT2});
        }
        else if (!multibody) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,dfdT,d2fdT2});
        }
        return null;
    }

    protected byte swap0n(byte i) {
        if (i == 0) return (byte)n;
        if (i == n) return 0;
        return i;
    }

    public int populateEFBonds(Graph g, List<ClusterBonds> allBonds, List<Double> weights, boolean swap) {
        boolean multiGraph = graphHasEdgeColor(g, mBond) || graphHasEdgeColor(g, mmBond);
        int rv = 0;
        if (allPermutations || multiGraph) {
            // multibody graph.  we need to generate all permutations in order
            // to get proper cancellation
            AllIsomorphs allIso = new AllIsomorphs();
            AllIsomorphsParameters allIsoParams = new AllIsomorphsParameters(true);
            Set<Graph> permutations = new HashSet<Graph>();
            Set<Graph> gSet = new HashSet<Graph>();
            gSet.add(g);
            if (graphHasEdgeColor(g, mmBond)) {
                gSet = substMinMulti(g);
            }
            if (NumRootNodes.value(g) > 0) {
                // we just want this graph and its cancelling graph
                permutations.addAll(gSet);
                Graph gCancel = cancelMap.get(g);
                if (gCancel != null) {
                    if (graphHasEdgeColor(gCancel, mmBond)) {
                        permutations.addAll(substMinMulti(gCancel));
                    }
                    else {
                        permutations.add(gCancel);
                    }
                }
            }
            else {
                // substitution generates lots of permtuations, there's no point
                // in doing all permutations of permutations
                IsoFree isofree = new IsoFree();
                gSet = isofree.apply(gSet, null);
                permutations = allIso.apply(gSet, allIsoParams);
            }

            rv = permutations.size();
            for (Graph gp : permutations) {
                ArrayList<int[]> fbonds = new ArrayList<int[]>();
                ArrayList<int[]> ebonds = new ArrayList<int[]>();
                int[][] mBonds = new int[n+1][0];
                List<List<Byte>> multiBiComponents = BCVisitor.getBiComponents(gp);
                boolean dupRoot = (!multiGraph || (multiGraph && doMultiFromPair)) && flex;
                int nn = dupRoot ? n+1 : n;
                for (List<Byte> comp : multiBiComponents) {
                    if (comp.size() == 1) continue;
                    if (!gp.hasEdge(comp.get(0), comp.get(1)) || gp.getEdge(comp.get(0), comp.get(1)).getColor() != mBond) {
                        // we've encountered e/f-bonds
                        // find all f-bonds within this component
                        for (int id1 = 0; id1<comp.size()-1; id1++) {
                            for (int id2 = id1+1; id2<comp.size(); id2++) {
                                byte n1 = comp.get(id1);
                                byte n2 = comp.get(id2);
                                if (!gp.hasEdge(n1, n2)) continue;
                                char edgeColor = gp.getEdge(n1, n2).getColor();
                                if (swap) n1 = swap0n(n1);
                                if (swap) n2 = swap0n(n2);
                                if (edgeColor == fBond) {
                                    fbonds.add(new int[]{n1,n2});
                                }
                                else if (edgeColor == eBond) {
                                    ebonds.add(new int[]{n1,n2});
                                }
                                else {
                                    throw new RuntimeException("oops, unknown bond "+edgeColor);
                                }
                            }
                        }
                        continue;
                    }
                    int groupID = -1;
                    // determine the groupID for this component
                    byte id0 = comp.get(0);
                    if (swap) id0 = swap0n(id0);
                    byte id1 = comp.get(1);
                    if (swap) id1 = swap0n(id1);
                    byte id2 = comp.get(2);
                    if (swap) id2 = swap0n(id2);
                    int size = comp.size();
                    if (comp.size() == 3) {
                        byte[] ids = new byte[]{id0, id1, id2};
                        java.util.Arrays.sort(ids);
                        groupID = tripletId(ids[0], ids[1], ids[2], nn);
                    }
                    else if (comp.size() == 4) {
                        byte id3 = comp.get(3);
                        if (swap) id3 = swap0n(id3);
                        byte[] ids = new byte[]{id0, id1, id2, id3};
                        java.util.Arrays.sort(ids);
                        groupID = quadId(ids[0], ids[1], ids[2], ids[3], nn);
                    }
                    else if (comp.size() == 5) {
                        byte id3 = comp.get(3);
                        if (swap) id3 = swap0n(id3);
                        byte id4 = comp.get(4);
                        if (swap) id4 = swap0n(id4);
                        byte[] ids = new byte[]{id0, id1, id2, id3, id4};
                        java.util.Arrays.sort(ids);
                        groupID = quintId(ids[0], ids[1], ids[2], ids[3], ids[4], nn);
                    }
                    else if (comp.size() == 6) {
                        byte id3 = comp.get(3);
                        if (swap) id3 = swap0n(id3);
                        byte id4 = comp.get(4);
                        if (swap) id4 = swap0n(id4);
                        byte id5 = comp.get(5);
                        if (swap) id5 = swap0n(id5);
                        byte[] ids = new byte[]{id0, id1, id2, id3, id4, id5};
                        java.util.Arrays.sort(ids);
                        groupID = sixId(ids[0], ids[1], ids[2], ids[3], ids[4], ids[5], nn);
                    }
                    int[] newGroups = new int[mBonds[size].length+1];
                    System.arraycopy(mBonds[size], 0, newGroups, 0, mBonds[size].length);
                    newGroups[newGroups.length-1] = groupID;
                    mBonds[size] = newGroups;
                }
                if (multiGraph) {
                    if (ebonds.size() > 0) {
                        allBonds.add(new ClusterBondsNonAdditive(nn, new int[][][]{fbonds.toArray(new int[0][0]), ebonds.toArray(new int[0][0])}, mBonds));
                    }
                    else {
                        allBonds.add(new ClusterBondsNonAdditive(nn, new int[][][]{fbonds.toArray(new int[0][0])}, mBonds));
                    }
                }
                else if (ebonds.size() > 0) {
                    allBonds.add(new ClusterBonds(nn, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0])}));
                }
                else {
                    allBonds.add(new ClusterBonds(nn, new int[][][]{fbonds.toArray(new int[0][0])}));
                }
                double w = gp.coefficient().getValue();
                if (dupRoot) {
                    // we'll get called twice.  once with swap, once without
                    w *= 0.5;
                }
                weights.add(w);
            }
        }
        else {
            ArrayList<int[]> ebonds = new ArrayList<int[]>();
            ArrayList<int[]> fbonds = new ArrayList<int[]>();
            rv = 1;
            for (Node node1 : g.nodes()) {
                for (Node node2 : g.nodes()) {
                    if (node1.getId() >= node2.getId()) continue;
                    if (g.hasEdge(node1.getId(), node2.getId())) {
                        byte n1 = node1.getId();
                        byte n2 = node2.getId();
                        char edgeColor = g.getEdge(n1, n2).getColor();
                        if (swap) n1 = swap0n(n1);
                        if (swap) n2 = swap0n(n2);
                        if (edgeColor == fBond) {
                            fbonds.add(new int[]{n1,n2});
                        }
                        else if (edgeColor == eBond) {
                            ebonds.add(new int[]{n1,n2});
                        }
                        else {
                            throw new RuntimeException("oops, unknown bond "+edgeColor);
                        }
                    }
                }
            }
            int gn = flex ? n+1 : n;
            if (ebonds.size() > 0) {
                allBonds.add(new ClusterBonds(gn, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0])}));
            }
            else {
                allBonds.add(new ClusterBonds(gn, new int[][][]{fbonds.toArray(new int[0][0])}));
            }
            double w = g.coefficient().getValue();
            if (flex) {
                w *= 0.5;
            }
            weights.add(w);
        }
        return rv;
    }
    
    public void populateEFdfdTBonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap) {
        ArrayList<int[]> fbonds = new ArrayList<int[]>();
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
        ArrayList<int[]> dfdTbonds = new ArrayList<int[]>();
        for (Node node1 : g.nodes()) {
            for (Node node2 : g.nodes()) {
                if (node1.getId() > node2.getId()) continue;
                if (g.hasEdge(node1.getId(), node2.getId())) {
                    byte n1 = node1.getId();
                    byte n2 = node2.getId();
                    if (swap) {
                        if (n1 == 0) n1 = (byte)n;
                        else if (n1 == n) n1 = (byte)0;
                        else if (n2 == 0) n2 = (byte)n;
                        else if (n2 == n) n2 = (byte)0;
                    }
                    char edgeColor = g.getEdge(node1.getId(), node2.getId()).getColor();
                    if (edgeColor == fBond) {
                        fbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == eBond) {
                        ebonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == 'a') {
                        dfdTbonds.add(new int[]{n1,n2});
                    }
                    else {
                        throw new RuntimeException("oops, unknown bond "+edgeColor);
                    }
                }
            }
        }
        if (!flex && n > 3) {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0]),dfdTbonds.toArray(new int[0][0])}));
        }
        else {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),dfdTbonds.toArray(new int[0][0])}));
        }

    }
    
    public void populateEFdfdTdf2dT2Bonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap) {
        ArrayList<int[]> fbonds = new ArrayList<int[]>();
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
        ArrayList<int[]> dfdTbonds = new ArrayList<int[]>();
        ArrayList<int[]> d2fdT2bonds = new ArrayList<int[]>();
        for (Node node1 : g.nodes()) {
            for (Node node2 : g.nodes()) {
                if (node1.getId() > node2.getId()) continue;
                if (g.hasEdge(node1.getId(), node2.getId())) {
                    byte n1 = node1.getId();
                    byte n2 = node2.getId();
                    if (swap) {
                        if (n1 == 0) n1 = (byte)n;
                        else if (n1 == n) n1 = (byte)0;
                        else if (n2 == 0) n2 = (byte)n;
                        else if (n2 == n) n2 = (byte)0;
                    }
                    char edgeColor = g.getEdge(node1.getId(), node2.getId()).getColor();
                    if (edgeColor == fBond) {
                        fbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == eBond) {
                        ebonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == 'a') {
                        dfdTbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == 'b') {
                        d2fdT2bonds.add(new int[]{n1,n2});
                    }
                    else {
                        throw new RuntimeException("oops, unknown bond "+edgeColor);
                    }
                }
            }
        }
        if (!flex && n > 3) {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0]),dfdTbonds.toArray(new int[0][0]),d2fdT2bonds.toArray(new int[0][0])}));
        }
        else {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),dfdTbonds.toArray(new int[0][0]),d2fdT2bonds.toArray(new int[0][0])}));
        }

    }

    public ClusterSumShell[] makeSingleVirialClusters(ClusterSum coreCluster, MayerFunction e, MayerFunction f) {
        if (p == null) {
            makeVirialDiagrams();
        }
        ArrayList<ClusterSumShell> allClusters = new ArrayList<ClusterSumShell>();
        Set<Graph> pn = getMSMCGraphs(true, false);
        IsBiconnected isBi = new IsBiconnected();
        if (doMinimalBC) {
            List<Double> weights = new ArrayList<Double>();
            // group all biconnected pairwise graphs into a single clusterSum
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            double leadingCoef = 0;
            for (Graph g : pn) {
                if (graphHasEdgeColor(g, mBond) || !isBi.check(g)) continue;
                if (!graphHasEdgeColor(g, eBond)) leadingCoef = g.coefficient().getValue();
                populateEFBonds(g, allBonds, weights, false);
                if (flex) {
                    populateEFBonds(g, allBonds, weights, true);
                }
            }
            double[] w = new double[weights.size()];
            for (int i=0; i<w.length; i++) {
                w[i] = weights.get(i) / leadingCoef;
            }
            if (n > 3 && !flex) {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,e}));
            }
            else if (allBonds.size() > 0) {
                // we might end up with nothing if we have flexOnlyCancel; that's OK
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}));
            }
        }
        for (Graph g : pn) {
            if (doMinimalBC && isBi.check(g)) continue;
            List<Double> weights = new ArrayList<Double>();
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            populateEFBonds(g, allBonds, weights, false);
            if (flex) {
                populateEFBonds(g, allBonds, weights, true);
            }
            if (flex && cancelMap.get(g) != null) {
                Graph cg = cancelMap.get(g);
                populateEFBonds(cg, allBonds, weights, false);
                populateEFBonds(cg, allBonds, weights, true);
            }
            double[] thisW = new double[weights.size()];
            double gCoef = g.coefficient().getValue();
            for (int i=0; i<thisW.length; i++) {
                thisW[i] = weights.get(i) / gCoef;
            }
            if (n > 3 && !flex) {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{f,e}));
            }
            else {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{f}));
            }
        }
        return allClusters.toArray(new ClusterSumShell[0]);
    }

    public Set<Graph> substMinMulti(Graph g) {
        if (!doMinimalMulti) {
            throw new RuntimeException("we need the minimal multi info to construct appropriate clusters");
        }
        Set<Graph>[] iMinMultiP = new Set[n+1];
        for (byte i=3; i<=n; i++) {
            iMinMultiP[i] = makeGraphList();
        }
        for (Graph gmm : minMultiP) {
            iMinMultiP[gmm.nodeCount()].add(gmm);
        }
        MulScalar mulScalar = new MulScalar();
        for (byte i=3; i<=n; i++) {
            iMinMultiP[i] = mulScalar.apply(iMinMultiP[i], new MulScalarParameters((int)SpecialFunctions.factorial(i),1-i));
        }
        BiComponentSubst biCompSubst = new BiComponentSubst();
        Set<Graph> gSubst = new HashSet<Graph>();
        gSubst.add(g);
        for (byte i=(byte)n; i>=3; i--) {
            Graph gComp = GraphFactory.createGraph(i);
            for (byte j=1; j<i; j++) {
                for (byte k=0; k<j; k++) {
                    gComp.putEdge(k,j);
                    gComp.getEdge(k,j).setColor(mmBond);
                }
            }
            gSubst = biCompSubst.apply(gSubst, new BiComponentSubstParameters(gComp, iMinMultiP[i], true, false));
        }
        return gSubst;
    }

    public ClusterSumMultibodyShell[] makeSingleVirialClustersMulti(ClusterSumMultibody coreCluster, MayerFunction f, MayerFunctionNonAdditive fMulti) {
        if (p == null) {
            makeVirialDiagrams();
        }
        List<ClusterSumMultibodyShell> allClusters = new ArrayList<ClusterSumMultibodyShell>();
        if (!doMinimalMulti) {
            throw new RuntimeException("we need the minimal multi info to construct appropriate clusters");
        }
        Set<Graph> pn = getMSMCGraphs(true, true);
        for (Graph g : pn) {
            if (multibody && flex && doMultiFromPair && NumRootNodes.value(g) > 1) continue;
            List<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            List<Double> weights = new ArrayList<Double>();
            populateEFBonds(g, allBonds, weights, false);
            if (doMultiFromPair && flex) {
                populateEFBonds(g, allBonds, weights, true);
            }
            double[] w = new double[weights.size()];
            double gCoef = g.coefficient().getValue();
            for (int i=0; i<w.length; i++) {
                w[i] = weights.get(i) / gCoef;
            }
            allClusters.add(new ClusterSumMultibodyShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}, new MayerFunctionNonAdditive[]{fMulti}));
        }
        return allClusters.toArray(new ClusterSumMultibodyShell[0]);
    }
    
    public Set<Graph> getExtraDisconnectedVirialGraphs() {
        if (p == null) {
            makeVirialDiagrams();
        }
        GraphList dpn = makeGraphList();
        for (Graph g : disconnectedP) {
            if (NumFieldNodes.value(g) == n) {
                dpn.add(g);
            }
        }
        if (multibody && doMultiFromPair) {
            for (Graph g : p) {
                if (NumFieldNodes.value(g) == n) {
                    Graph cancelGraph = cancelMap.get(g);
                    if (cancelGraph != null && NumRootNodes.value(g) > 0 && NumRootNodes.value(cancelGraph) > 0) {
                        // FIXME
                        // this is a diagram we list in p, but which we don't compute in MSMC...
                        // it's disconnectedP for our purposes.
                        dpn.add(g);
                    }
                }
            }
        }
        return dpn;
    }

    public Set<Graph> getSplitDisconnectedVirialGraphs(Graph g) {
        SplitGraph splitGraph = new SplitGraph();
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        Property happyArticulation = new ArticulatedAt0(doExchange, multibody ? mmBond : '0');
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), happyArticulation);
        // we want gSplit unsorted
        Set<Graph> gSplit = new GraphList(null);
        Set<Graph> gSplit1 = splitGraph.apply(g);
        for (Graph gs : gSplit1) {
            // the graph we get from splitting might not be in our preferred bonding arrangement
            Graph gsmax = maxIsomorph.apply(gs, mip);
            gSplit.add(gsmax);
        }
        return gSplit;
    }

    public static GraphList makeGraphList() {
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodesExchange());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList graphList = new GraphList(comp);
        return graphList;
    }

    public static HashMap<Character,Integer> initMetaDataComparator() {
        final HashMap<Character,Integer> colorOrderMap = new HashMap<Character,Integer>();
        if (MetadataImpl.metaDataComparator != null) {
            System.err.println("metaDataComparator already created");
            return null;
        }
        
        MetadataImpl.metaDataComparator = new Comparator<Metadata>() {
            public int compare(Metadata m1, Metadata m2) {
                Integer o1 = colorOrderMap.get(m1.getColor());
                Integer o2 = colorOrderMap.get(m2.getColor());
                if (o1 == o2) {
                    return m1.getColor() > m2.getColor() ? 1 : -1;
                }
                if (o1 == null) {
                    if (o2 == null) {
                        return m1.getColor() > m2.getColor() ? 1 : -1;
                    }
                    return -1;
                }
                if (o2 == null) {
                    return 1;
                }
                if (o1 != o2) {
                    return o1 > o2 ? 1 : -1;
                }
                if (m1.getType() != m2.getType()) {
                    return m1.getType() > m2.getType() ? 1 : -1;
                }
                return 0;
            }
        };
        return colorOrderMap;
    }
    
    public void makeRhoDiagrams() {
        flexColors = new char[0];
        if (flex) {
           flexColors = new char[]{nodeColor};
        }
        if (doExchange) {
            // all colors are flexible
            flexColors = new char[2*n];
            for (int i=0; i<n; i++) {
                flexColors[i] = (char)('A' + i);
                flexColors[n+i] = (char)('a' + i);
            }
        }

        final HashMap<Character,Integer> colorOrderMap = initMetaDataComparator();

        GraphList topSet = makeGraphList();

        char oneBond = 'o';
        Metadata.COLOR_MAP.put(eBond, "red");
        Metadata.COLOR_MAP.put(fBond, "green");
        Metadata.COLOR_MAP.put(mBond, "blue");
        Metadata.COLOR_MAP.put(mmBond, "orange");
        Metadata.COLOR_MAP.put(fmBond, "black");
        Metadata.COLOR_MAP.put(efbcBond, "fuchsia");
        Metadata.COLOR_MAP.put(excBond, "red");
        Metadata.COLOR_MAP.put(ffBond, "green");
        Metadata.COLOR_MAP.put(mxcBond, "blue");
        Metadata.DASH_MAP.put(excBond, 3);
        Metadata.DASH_MAP.put(ffBond, 3);
        Metadata.DASH_MAP.put(mxcBond, 3);
        Metadata.DASH_MAP.put(MxcBond, 3);
        if (colorOrderMap != null) {
            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(efbcBond, 2);
            colorOrderMap.put(fBond, 3);
            colorOrderMap.put(eBond, 4);
        }
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        
        if (doShortcut && !multibody && !doKeepEBonds) {
            // just take lnfXi to be the set of connected diagrams
            IsConnected isCon = new IsConnected();
            for (int i=1; i<n+1; i++) {
                GraphIterator iter = new PropertyFilter(new StoredIterator((byte)i), isCon);
                msp = new MulScalarParameters(1, (int)SpecialFunctions.factorial(i));
                while (iter.hasNext()) {
                    Graph g = iter.next();
                    for (Edge e : g.edges()) {
                        e.setColor(fBond);
                    }
                    lnfXi.add(mulScalar.apply(g, msp));
                }
            }
        }
        else {

            Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
            for (byte i=1; i<n+1; i++) {
                Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                for (Edge e : g.edges()) {
                    e.setColor(eBond);
                }
                g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                eXi.add(g);
    
                if (multibody && i>2 && !doKeepEBonds) {
                    g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                    g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                    for (Edge e : g.edges()) {
                        e.setColor(mBond);
                    }
                    eXi.add(g);
                }
            }
    
            if (isInteractive) {
                System.out.println("Xi");
                topSet.addAll(eXi);
                for (Graph g : topSet) {
                    System.out.println(g);
                }
            }
    
            Set<Graph> fXi;
            if (doShortcut && !doKeepEBonds) {
                fXi = new HashSet<Graph>();
                for (byte i=1; i<n+1; i++) {
                    GraphIterator iter = new StoredIterator(i);
                    msp = new MulScalarParameters(1, (int)SpecialFunctions.factorial(i));
                    while (iter.hasNext()) {
                        Graph g = iter.next();
                        for (Edge e : g.edges()) {
                            e.setColor(fBond);
                        }
                        fXi.add(mulScalar.apply(g, msp));
                    }
                    
                    if (multibody && i>2) {
                        Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                        g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                        for (Edge e : g.edges()) {
                            e.setColor(mBond);
                        }
                        fXi.add(g);
                    }
                }
            }
            else if (!doKeepEBonds) {
                Split split = new Split();
                SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
                Set<Graph> setOfSubstituted = split.apply(eXi, bonds);
        
                DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
                DeleteEdge deleteEdge = new DeleteEdge();
                //set of full star diagrams with f bonds
                fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
            }
            else {
                // fXi is the set we'll use from here on.  just take it to be
                // the eBond representation
                fXi = eXi;
            }

            if (isInteractive) {
                System.out.println("\nXi with f bonds");
                topSet.clear();
                topSet.addAll(fXi);
                for (Graph g : topSet) {
                    System.out.println(g);
                }
            }
            
            Set<Graph> fXipow = new HashSet<Graph>();
            fXipow.addAll(fXi);
            for (int i=1; i<n+1; i++) {

                lnfXi.addAll(fXipow);
                lnfXi = isoFree.apply(lnfXi, null);
                msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
                fXipow = isoFree.apply(mulScalar.apply(mulFlex.apply(fXipow, fXi, mfp), msp), null);
            }

            if (doHB) {
                if (isInteractive) {
                    topSet.clear();
                    topSet.addAll(lnfXi);
                    System.out.println("\nlnfXi before HB substitution");
                    for (Graph g : topSet) {
                        System.out.println(g);
                    }
                    ClusterViewer.createView("lnfXi before HB", topSet);
                }
                Set<Graph>[] lnXin = new Set[n+1];
                for  (int i=1; i<=n; i++) {
                    lnXin[i] = new HashSet<Graph>();
                }
                Coefficient[] coef = new Coefficient[n+1];
                Graph[] bigGraph = new Graph[n+1];
                for (Graph g : lnfXi) {
                    byte nc = g.nodeCount();
                    if (nc>1 && g.edgeCount() == nc*(nc-1)/2) {
                        coef[nc] = g.coefficient().copy();
                        bigGraph[nc] = g.copy();
                        g = g.copy();
                        bigGraph[nc].coefficient().divide(coef[nc]);
                        for (Edge e : g.edges()) {
                            e.setColor(bBond);
                        }
                        g.coefficient().multiply(new CoefficientImpl(-1));
                    }
                    lnXin[g.nodeCount()].add(g);
                }
                BiComponentSubst biCompSubst = new BiComponentSubst();
                for  (int i=2; i<n; i++) {
                    Coefficient iCoef = new CoefficientImpl(-1);
                    iCoef.divide(coef[i]);
                    lnXin[i] = mulScalar.apply(lnXin[i], new MulScalarParameters(iCoef));
                    BiComponentSubstParameters bcsp = new BiComponentSubstParameters(bigGraph[i], lnXin[i], true, false);
                    for (int j=i+1; j<=n; j++) {
                        lnXin[j] = isoFree.apply(biCompSubst.apply(lnXin[j], bcsp), null);
                    }
                }
                lnfXi.clear();
                lnfXi.addAll(lnXin[1]);
                for (int i=2; i<=n; i++) {
                    for (Graph g : lnXin[i]) {
                        byte nc = g.nodeCount();
                        if (g.edgeCount() == nc*(nc-1)/2) {
                            Graph gb = g.copy();
                            lnfXi.add(gb);
                            for (Edge e : g.edges()) {
                                e.setColor(eBond);
                            }
                            if (i==n) {
                                g.coefficient().multiply(new CoefficientImpl(-1));
                                gb.coefficient().multiply(new CoefficientImpl(-1));
                            }
                            else {
                                g.coefficient().multiply(coef[i]);
                                gb.coefficient().multiply(coef[i]);
                            }
                        }
                        else if (i<n) {
                            g.coefficient().multiply(new CoefficientImpl(-1));
                            g.coefficient().multiply(coef[i]);
                        }
                    }
                }
                fullLnXi = makeGraphList();
                for (int i=1; i<=n; i++) {
                    fullLnXi.addAll(lnXin[i]);
                }
                if (isInteractive) {
                    topSet.clear();
                    topSet.addAll(fullLnXi);
                    System.out.println("\nfull lnXi");
                    for (Graph g : topSet) {
                        System.out.println(g);
                    }
                    ClusterViewer.createView("full lnXi", topSet);
                }
                
            }
        }

        if (isInteractive) {
            topSet.clear();
            topSet.addAll(lnfXi);
            System.out.println("\nlnfXi");
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("lnfXi", topSet);
        }
    

        DifByNode opzdlnXidz = new DifByNode();
        DifParameters difParams = new DifParameters(nodeColor);
        rho = isoFree.apply(opzdlnXidz.apply(lnfXi, difParams), null);

        Relabel relabel = new Relabel();
        Set<Graph> rhoNew = new HashSet<Graph>();
        for (Graph g : rho) {
            if (g.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT) {
                rhoNew.add(g);
                continue;
            }
            
            // find the root node
            byte rootId = -1;
            for (byte i=1; i<g.nodeCount(); i++) {
                if (g.getNode(i).getType() == Metadata.TYPE_NODE_ROOT) {
                    rootId = i;
                    break;
                }
            }
            
            // the graph's root node is not at 0.
            // we want the diagram with the root node at 0, so permute
            byte[] newLabels = new byte[g.nodeCount()];
            for (int i=0; i<newLabels.length; i++) {
                newLabels[i] = (byte)i;
            }
            newLabels[0] = rootId;
            newLabels[rootId] = 0;
            RelabelParameters rp = new RelabelParameters(newLabels);
            rhoNew.add(relabel.apply(g, rp));
        }
        rho = rhoNew;

        if (isInteractive) {
            System.out.println("\nrho");
            topSet.clear();
            topSet.addAll(rho);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("rho", topSet);
        }
    }
    
    public void makeVirialDiagrams() {
        if (p != null) return;

        flexColors = new char[0];
        if (flex) {
           flexColors = new char[]{nodeColor};
        }
        if (doExchange) {
            // all colors are flexible
            flexColors = new char[2*n];
            for (int i=0; i<n; i++) {
                flexColors[i] = (char)('A' + i);
                flexColors[n+i] = (char)('a' + i);
            }
        }

        final HashMap<Character,Integer> colorOrderMap = new HashMap<Character,Integer>();
        MetadataImpl.metaDataComparator = new Comparator<Metadata>() {

            public int compare(Metadata m1, Metadata m2) {
                if (m1.getType() != m2.getType()) {
                    return m1.getType() > m2.getType() ? 1 : -1;
                }
                Integer o1 = colorOrderMap.get(m1.getColor());
                Integer o2 = colorOrderMap.get(m2.getColor());
                if (o1 == o2) {
                    return m1.getColor() > m2.getColor() ? 1 : -1;
                }
                if (o1 == null) {
                    if (o2 == null) {
                        return m1.getColor() > m2.getColor() ? 1 : -1;
                    }
                    return -1;
                }
                if (o2 == null) {
                    return 1;
                }
                return o1 > o2 ? 1 : (o1 < o2 ? -1 : 0);
            }
        };

        GraphList topSet = makeGraphList();

        char oneBond = 'o';
        IsoFree isoFree = new IsoFree();

        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(mBond, 1);
        colorOrderMap.put(mmBond, 2);
        colorOrderMap.put(eBond, 3);
        colorOrderMap.put(fBond, 4);
        colorOrderMap.put(efbcBond, 5);
        colorOrderMap.put(ffBond, 6);
        colorOrderMap.put(mxcBond, 7);
        colorOrderMap.put(MxcBond, 8);
        colorOrderMap.put(excBond, 9);

        Property happyArticulation = new ArticulatedAt0(doExchange, multibody ? mmBond : '0');

        if (doShortcut && !multibody && !flex) {
            lnfXi = new HashSet<Graph>();

            // skip directly to p diagrams
            p = new HashSet<Graph>();
            IsBiconnected isBi = new IsBiconnected();
            for (int i=1; i<n+1; i++) {
                GraphIterator iter = new PropertyFilter(new StoredIterator((byte)i), isBi);
                msp = new MulScalarParameters(i>1 ? 1-i : 1, (int)SpecialFunctions.factorial(i));
                while (iter.hasNext()) {
                    Graph g = iter.next();
                    for (Edge e : g.edges()) {
                        e.setColor(fBond);
                    }
                    p.add(mulScalar.apply(g, msp));
                }
            }
        }
        else {
            if (rho == null) {
                makeRhoDiagrams();
            }
            
            @SuppressWarnings("unchecked")
            HashSet<Graph>[] allRho = new HashSet[n+1];
            for (int i=0; i<n+1; i++) {
                allRho[i] = new HashSet<Graph>();
            }
            for (Graph g : rho) {
                allRho[g.nodeCount()].add(g);
            }
            
            Set<Graph> z = new HashSet<Graph>();
            
            // r = z + b*z^2 + c*z^3 + d*z^4 + e*z^5
            // z = r - b*r^2 + (2 b^2 - c) r^3 + (-5 b^3 + 5*b*c - d) r^4
            //     + (14 b^4 - 21 b^2*c + 3 c^2 + 6 b*d - e) r^5
            //     + (-42 b^5 + 84 b^3*c - 28 (b*c^2 + b^2*d) + 7 (b*e + c*d) - f) r^6
            // or
            // z = r - b*z^2 - c*z^3 - d*z^4 - e*z^5
            // z0 = r
            // z1 = r - b*z0^2
            //    = r - b*r^2
            // z2 = r - b*z1^2 - c*z1^3
            //    = r - b*(r-b*r^2)^2 - c*(r-b*r^2)^3
            //    = r - b*r^2 + 2*b^2*r^3 - c*r^3 + O[r^4]
            // etc

            MulFlexibleParameters mfpnm1 = MulFlexibleParameters.makeParameters(flexColors, (byte)(n-1));
            MulFlexible mulFlex = new MulFlexible();
            z.addAll(allRho[1]);
            for (int i=2; i<n+1; i++) {
                @SuppressWarnings("unchecked")
                Set<Graph>[] zPow = new HashSet[n+1];
                zPow[1] = new HashSet<Graph>();
                zPow[1].addAll(z);
                for (int j=2; j<i+1; j++) {
                    zPow[j] = new HashSet<Graph>();
                    zPow[j] = isoFree.apply(mulFlex.apply(zPow[j-1], z, mfpnm1), null);
                }
                z = new HashSet<Graph>();
                z.addAll(allRho[1]);
                msp = new MulScalarParameters(new CoefficientImpl(-1,1));
                for (int j=2; j<i+1; j++) {
                    z.addAll(mulScalar.apply(mulFlex.apply(allRho[j], zPow[j], mfpnm1), msp));
                }
            }

            z = isoFree.apply(z, null);
            if (isInteractive) {
                System.out.println("\nz");
                topSet.clear();
                topSet.addAll(z);
                for (Graph g : topSet) {
                    System.out.println(g);
                }
            }
            
            Decorate decorate = new Decorate();
            DecorateParameters dp = new DecorateParameters(nodeColor, mfp);
            p = decorate.apply(lnfXi, z, dp);
            p = isoFree.apply(p, null);

            cancelMap = new HashMap<Graph,Graph>();

//            if (isInteractive && multibody) {
//                System.out.println("\nP (full-multi)");
//                topSet.clear();
//                topSet.addAll(p);
//                for (Graph g : topSet) {
//                    System.out.println(g);
//                }
//                ClusterViewer.createView("P (full-multi)", topSet);
//            }
            
            if (doMinimalMulti && !doKeepEBonds) {
                minMultiP = makeGraphList();
                // pluck multibody graphs from lnfXi
                for (Graph g : lnfXi) {
                    if (graphHasEdgeColor(g, mBond)) {
                        // multibond graphs in lnfXi are part of min-multi
                        minMultiP.add(mulScalar.apply(g, new MulScalarParameters(1-g.nodeCount(),1)));
                    }
                }
                
                if (isInteractive) {
                    topSet.clear();
                    topSet.addAll(minMultiP);
                    ClusterViewer.createView("multiP", topSet);
                }
            }
            
            if (doMultiFromPair) {
                // we're going to use the pairwise additive diagrams to determine
                // appropriate singly-connected nonadditive diagrams.
                // we'll keep our disconnected nonadditive diagrams
                IsConnected isConnected = new IsConnected();
                Set<Graph> newP = makeGraphList();
                for (Graph g : p) {
                    if (!graphHasEdgeColor(g, mBond)) {
                        newP.add(g);
                    }
                }
                // pluck multibody graphs from lnfXi
                for (Graph g : lnfXi) {
                    if (graphHasEdgeColor(g, mBond) || !isConnected.check(g)) {
                        continue;
                    }
                    // a connected purely pairwise graph.
                    // we are interested if g is made up of fully connected bicomponents
                    List<List<Byte>> biComps = BCVisitor.getBiComponents(g);
                    List<List<Byte>> happyComps = new ArrayList<List<Byte>>();
                    for (List<Byte> biComp : biComps) {
                        if (biComp.size() > 2) {
                            boolean happy = true;
outer:                      for (int i=1; i<biComp.size(); i++) {
                                for (int j=0; j<i; j++) {
                                    if (!g.hasEdge(biComp.get(i), biComp.get(j))) {
                                        happy = false;
                                        break outer;
                                    }
                                }
                            }
                            if (happy) {
                                // this component is fully connected
                                happyComps.add(biComp);
                            }
                        }
                    }
                    if (happyComps.size() == 0) {
                        // we found no fully-connected components of size 3 or greater 
                        continue;
                    }
                    byte gn = g.nodeCount();
                    msp = new MulScalarParameters(gn-1,1);
                    if (g.edgeCount() == gn*(gn-1)/2) {
                        // a fully connected diagram.  take this one as min-multi
                        msp = new MulScalarParameters(1-gn,1);
                        g = mulScalar.apply(g, msp);
                        Graph g2 = g.copy();
                        for (Edge e : g2.edges()) {
                            e.setColor(mmBond);
                        }
                        newP.add(g2);
                        continue;
                    }
                    g = mulScalar.apply(g, msp);
                    // now replace each fully connected component (of size 3+) with a full-multi-component
                    Set<Graph> gSet = new HashSet<Graph>();
                    gSet.add(g);
                    for (List<Byte> biComp : happyComps) {
                        Set<Graph> newGSet = new HashSet<Graph>();
                        // we want one copy without replacement and one with
                        // at the end, we'll have one graph without replacement and
                        // throw it away (that graph is already in p)
                        newGSet.addAll(gSet);
                        for (Graph g2 : gSet) {
                            g2 = g2.copy();
                            for (int i=1; i<biComp.size(); i++) {
                                for (int j=0; j<i; j++) {
                                    g2.getEdge(biComp.get(j), biComp.get(i)).setColor(fmBond);
                                }
                            }
                            newGSet.add(g2);
                        }
                        gSet = newGSet;
                    }
                    Set<Graph> newGSet = new HashSet<Graph>();
                    for (Graph g2 : gSet) {
                        if (graphHasEdgeColor(g2, fmBond) || graphHasEdgeColor(g2, mmBond)) {
                            newGSet.add(g2);
                        }
                    }
                    newP.addAll(isoFree.apply(newGSet, null));
                }
                if (isInteractive) {
                    System.out.println("\nP (min-multi)");
                    topSet.clear();
                    topSet.addAll(minMultiP);
                    for (Graph g : topSet) {
                        System.out.println(g);
                    }
                    ClusterViewer.createView("P (min-multi)", topSet);
                }
                trueMultiP = makeGraphList();
                for (Graph g : p) {
                    if (graphHasEdgeColor(g, mBond)) {
                        trueMultiP.add(g);
                    }
                }
                p = newP;
                // now we have p in terms of min-multi fully connected graphs and
                // singly-connected graphs with full-multi components.  full-multi
                // components represent all non-additive diagrams in P at that order.
                // we want P in terms of min-multi graphs, and so we need to substitute.
                // our expression at 3rd order is correct (min = full) and so we can
                // substitute our ith order expression (starting at 3rd) to an (i+1)th
                // order expression of min-multi graphs.  and just continue up to nth
                // order.
                fullMultiP = makeGraphList();
                @SuppressWarnings("unchecked")
                Set<Graph>[] iFullMultiP = new Set[n+1];
                for (byte i=3; i<=n; i++) {
                    iFullMultiP[i] = new HashSet<Graph>();
                }
                for (Graph g : p) {
                    if (graphHasEdgeColor(g, mmBond) || graphHasEdgeColor(g, fmBond)) {
                        fullMultiP.add(g);
                        byte nc = g.nodeCount();
                        iFullMultiP[nc].add(g);
                    }
                }
                BiComponentSubst biCompSubst = new BiComponentSubst();
                for (byte i=3; i<=n; i++) {
                    Graph gComp = GraphFactory.createGraph(i);
                    for (byte j=1; j<i; j++) {
                        for (byte k=0; k<j; k++) {
                            gComp.putEdge(k,j);
                            gComp.getEdge(k,j).setColor(fmBond);
                        }
                    }
                    
                    Set<Graph> newFullMultiP = isoFree.apply(biCompSubst.apply(fullMultiP, new BiComponentSubstParameters(gComp, mulScalar.apply(iFullMultiP[i], new MulScalarParameters((int)SpecialFunctions.factorial(i), 1-i)), true, true)), null);
                    fullMultiP.clear();
                    fullMultiP.addAll(newFullMultiP);
                    for (byte j=i; j<=n; j++) {
                        iFullMultiP[j].clear();
                    }
                    for (Graph g : fullMultiP) {
                        byte nc = g.nodeCount();
                        if (nc >= i) {
                            iFullMultiP[nc].add(g);
                        }
                    }
                }
                MaxIsomorph maxIsomorph = new MaxIsomorph();
                MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), happyArticulation);
                Set<Graph> newFullMultiP = makeGraphList();
                newFullMultiP.addAll(maxIsomorph.apply(fullMultiP, mip));
                fullMultiP = newFullMultiP;
                
                if (isInteractive) {
                    System.out.println("\nP (full-multi)");
                    topSet.clear();
                    topSet.addAll(fullMultiP);
                    for (Graph g : topSet) {
                        System.out.println(g);
                    }
                    ClusterViewer.createView("P (full-multi)", topSet);
                }
                newP = makeGraphList();
                for (Graph g : p) {
                    if (!graphHasEdgeColor(g, mmBond) && !graphHasEdgeColor(g, fmBond)) {
                        // grab all the pairwise bonds
                        newP.add(g);
                    }
                }
                // now add in the multi bonds
                newP.addAll(fullMultiP);
                p = newP;

                if (flex) {
                    if (!doDisconnectedMatching) {
                        throw new RuntimeException("It doesn't make any sense to do multi-from-pair for a flex model without disconnected matching");
                    }
                    // now we need to add back in the disconnect P diagrams we originally had
                    for (byte i=(byte)n; i>=3; i--) {
                        Graph gComp = GraphFactory.createGraph(i);
                        for (byte j=1; j<i; j++) {
                            for (byte k=0; k<j; k++) {
                                gComp.putEdge(k,j);
                                gComp.getEdge(k,j).setColor(mmBond);
                            }
                        }
                        
                        // the original diagrams are m-bonds (sigh)
                        Set<Graph> mpSubst = new HashSet<Graph>();
                        for (Graph g : minMultiP) {
                            if (g.nodeCount() == i) {
                                g = mulScalar.apply(g, new MulScalarParameters((int)SpecialFunctions.factorial(i),1-i));
                                if (g.edgeCount() == i*(i-1)/2) {
                                    gComp = g.copy();
                                    for (byte j=1; j<i; j++) {
                                        for (byte k=0; k<j; k++) {
                                            g.getEdge(k,j).setColor(mmBond);
                                        }
                                    }
                                }
                                else {
                                    g = mulScalar.apply(g, new MulScalarParameters(-1,1));
                                }
                                mpSubst.add(g);
                            }
                        }
                        
                        trueMultiP = isoFree.apply(biCompSubst.apply(trueMultiP, new BiComponentSubstParameters(gComp, mpSubst, true, true)), null);
                    }
                    if (isInteractive) {
                        System.out.println("\nP (disconnected multi)");
                        topSet.clear();
                        topSet.addAll(trueMultiP);
                        for (Graph g : topSet) {
                            System.out.println(g);
                        }
                        ClusterViewer.createView("P (disconnected multi)", topSet);
                    }
                }
            }
            
            // clear these out -- we don't need them and (in extreme cases) we might need the memory
            lnfXi.clear();
            rho.clear();
            z.clear();
        }


        if (doKeepEBonds) {
            System.out.println("\nPe");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pe", topSet);
        }

        if (doExchange) {
            
            ExchangeSplit exchangeSplit = new ExchangeSplit();
            ExchangeSplitParameters esp = new ExchangeSplitParameters(mfp, (byte)n, eBond, excBond, doNegativeExchange);
            p = exchangeSplit.apply(p, esp);

            if (multibody) {
                // take every e+exc component of size >= 3.  replace each
                // component, c, with c+cm where cm has e=>m and exc=> mxc
                Set<Graph> newP = new HashSet<Graph>();
                for (Graph g : p) {
                    HashSet<Graph> gMultiSplit = new HashSet<Graph>();
                    gMultiSplit.add(g.copy());
                    List<List<Byte>> components = CVisitor.getComponents(g);
                    for (int i=0; i<components.size(); i++) {
                        List<Byte> comp = components.get(i);
                        if (comp.size() < 3) continue;
                        HashSet<Graph> gMultiNew = new HashSet<Graph>();
                        for (Graph gp : gMultiSplit) {
                            Graph gm = gp.copy();
                            for (byte inode1 = 0; inode1<comp.size()-1; inode1++) {
                                byte node1 = comp.get(inode1);
                                for (byte inode2 = (byte)(inode1+1); inode2<comp.size(); inode2++) {
                                    byte node2 = comp.get(inode2);
                                    Edge edge = gm.getEdge(node1, node2);
                                    if (edge.getColor() == eBond) {
                                        edge.setColor(mBond);
                                    }
                                    else {
                                        edge.setColor(mxcBond);
                                    }
                                }
                            }
                            gMultiNew.add(gm);
                        }
                        gMultiSplit.addAll(gMultiNew);
                    }
                    newP.addAll(gMultiSplit);
                }
                p = newP;
            }

            if (doExchangeF) {
                Set<Graph> newP = new HashSet<Graph>();
                for (Graph g : p) {
                    // for exchanged groups, we replace e1*e2*e3*e4 => f + 1
                    // where the e bonds are between the exchanged group and
                    // some other group.  The f bond is represented by a
                    // product of "F" bonds.
                    Graph gOnlyExc = g.copy();
                    for (Edge e : g.edges()) {
                        if (e.getColor() != excBond) {
                            gOnlyExc.deleteEdge(e.getId());
                        }
                    }
                    HashSet<Graph> gfSplit = new HashSet<Graph>();
                    gfSplit.add(g.copy());
                    List<List<Byte>> components = CVisitor.getComponents(gOnlyExc);
                    for (int i=0; i<components.size(); i++) {
                        List<Byte> comp = components.get(i);
                        if (comp.size() < 2) continue;
                        HashSet<Graph> gfSplitNew = null;
                        byte node0 = comp.get(0);
                        do {
                            gfSplitNew = new HashSet<Graph>();
                            for (Graph gf : gfSplit) {
                                // find all e-bonds between one of our exchange nodes and another component.
                                // create 2 graphs.
                                //   one replaces all e-bonds with any of our exchange nodes and that node with an ffBond
                                //   one deletes all e-bonds between any of our exchange nodes and that node
                                for (byte j=0; j<components.size(); j++) {
                                    if (i==j) continue;
                                    byte jNode0 = components.get(j).get(0);
                                    if (!gf.hasEdge(node0, jNode0) || gf.getEdge(node0,jNode0).getColor() == ffBond) continue;
                                    Graph gff = gf.copy();
                                    for (byte iNode : comp) {
                                        for (byte jNode : components.get(j)) {
                                            gff.getEdge(iNode,jNode).setColor(ffBond);
                                        }
                                    }
                                    gfSplitNew.add(gff);
                                    Graph gd = gf.copy();
                                    for (byte iNode : comp) {
                                        for (byte jNode : components.get(j)) {
                                            gd.deleteEdge(iNode,jNode);
                                        }
                                    }
                                    gfSplitNew.add(gd);
                                    break;
                                }
                            }
                            if (gfSplitNew.size() > 0) {
                                gfSplit = gfSplitNew;
                            }
                        } while (gfSplitNew.size() > 0);
                    }
                    newP.addAll(gfSplit);
                }
                p = newP;
            }
            
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), happyArticulation);
            p = maxIsomorph.apply(isoFree.apply(p, null), mip);
            
            System.out.println("\nPxc");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pxc", topSet);
        }
        else if (doKeepEBonds && multibody) {
            // take every e+exc component of size >= 3.  replace each
            // component, c, with c+cm where cm has e=>m and exc=> mxc
            Set<Graph> newP = new HashSet<Graph>();
            for (Graph g : p) {
                HashSet<Graph> gMultiSplit = new HashSet<Graph>();
                gMultiSplit.add(g.copy());
                List<List<Byte>> components = CVisitor.getComponents(g);
                for (int i=0; i<components.size(); i++) {
                    List<Byte> comp = components.get(i);
                    if (comp.size() < 3) continue;
                    HashSet<Graph> gMultiNew = new HashSet<Graph>();
                    for (Graph gp : gMultiSplit) {
                        Graph gm = gp.copy();
                        for (byte inode1 = 0; inode1<comp.size()-1; inode1++) {
                            byte node1 = comp.get(inode1);
                            for (byte inode2 = (byte)(inode1+1); inode2<comp.size(); inode2++) {
                                byte node2 = comp.get(inode2);
                                Edge edge = gm.getEdge(node1, node2);
                                if (edge.getColor() == eBond) {
                                    edge.setColor(mBond);
                                }
                            }
                        }
                        gMultiNew.add(gm);
                    }
                    gMultiSplit.addAll(gMultiNew);
                }
                newP.addAll(gMultiSplit);
            }
            p = newP;
        }

        if (doKeepEBonds) {
            // *now* switch the e-bonds back to f-bonds
            Split split = new Split();
            SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
            Set<Graph> newP = split.apply(p, bonds);
    
            DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
            DeleteEdge deleteEdge = new DeleteEdge();
            Unfactor unfactor = new Unfactor();
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), happyArticulation);

            p = maxIsomorph.apply(isoFree.apply(unfactor.apply(deleteEdge.apply(newP, deleteEdgeParameters), mfp), null), mip);
        }
        
        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), happyArticulation);
        if (doExchangeCondensing) {
            // we want to do this early so that we have simpler diagrams that are less likely to confuse
            // the rest of the operations.  However, the diagrams that result from this operation
            // have the root points (TYPE_NODE_ROOT) moved out of the exchange components and need
            // to stay that way, so we have to be careful that any operations that would move root points
            // happen before this.  Quite fragile.
            CondenseExchange condenser = new CondenseExchange();
            CondenseExchangeParameters condenseParams = new CondenseExchangeParameters(ffBond, excBond, fBond, 'A');
            p = condenser.apply(p, condenseParams);
            mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), happyArticulation);
            p = maxIsomorph.apply(isoFree.apply(p, null), mip);

            if (multibody) {
                condenseParams = new CondenseExchangeParameters(mBond, mxcBond, mBond, 'M');
                p = condenser.apply(p, condenseParams);
                p = maxIsomorph.apply(isoFree.apply(p, null), mip);
            }

            System.out.println("\nPc");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pc", topSet);
        }
        
        Set<Graph> newP = new HashSet<Graph>();
        // attempt to factor any graphs with an articulation point
        disconnectedP = new HashSet<Graph>();
        if (!flex) {

            if (doDisconnectedMatching) {
                Factor factor = new Factor();

                for (Graph g : p) {
                    if (doMultiFromPair && (graphHasEdgeColor(g, mmBond) || graphHasEdgeColor(g, fmBond))) {
                        newP.add(g);
                        continue;
                    }
                    boolean ap = hap.check(g);
                    boolean con = hap.isConnected();
                    if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                        Graph gf = factor.apply(g, mfp);
                        newP.add(gf);
                    }
                    else {
                        newP.add(g);
                    }
                }
                // we've split the graphs up; now combine them if possible
                p = isoFree.apply(newP, null);
                // now recondense them (only relevant for multibody interactions)
                Unfactor unfactor = new Unfactor();
                p = isoFree.apply(unfactor.apply(p, mfp), null);
            }
            // perform Ree-Hoover substitution (brute-force)
            if (doReeHoover) {
                if (doShortcut && !multibody) {
                    ReeHoover reeHoover = new ReeHoover();
                    p = reeHoover.apply(p, new ReeHooverParameters(eBond));
                }
                else {
                    char nfBond = 'F';
                    SplitOneParametersBC splitOneParameters = new SplitOneParametersBC(fBond, eBond, nfBond);
                    SplitOneBiconnected splitOneBC = new SplitOneBiconnected();
                    msp = new MulScalarParameters(-1, 1);
                    newP.clear();
                    for (Graph g : p) {
                        Set<Graph> gSet = splitOneBC.apply(g, splitOneParameters);
                        for (Graph g2 : gSet) {
                            boolean even = true;
                            for (Edge e : g2.edges()) {
                                if (e.getColor() == nfBond) {
                                    even = !even;
                                    e.setColor(fBond);
                                }
                            }
                            if (!even) {
                                g2 = mulScalar.apply(g2, msp);
                            }
                            newP.add(g2);
                        }
                    }
                    p = isoFree.apply(newP, null);
                    newP.clear();
                }
            }
            newP.clear();
            newP.addAll(maxIsomorph.apply(p, mip));
            p = newP;
        }
        else {

            // perform Ree-Hoover substitution (brute-force)
            if (doReeHoover) {
                char nfBond = 'Z';
                SplitOneParametersBC splitOneParameters = new SplitOneParametersBC(fBond, eBond, nfBond);
                SplitOneBiconnected splitOneBC = new SplitOneBiconnected();
                msp = new MulScalarParameters(-1, 1);
                newP.clear();
                for (Graph g : p) {
                    Set<Graph> gSet = splitOneBC.apply(g, splitOneParameters);
                    for (Graph g2 : gSet) {
                        boolean even = true;
                        for (Edge e : g2.edges()) {
                            if (e.getColor() == nfBond) {
                                even = !even;
                                e.setColor(fBond);
                            }
                        }
                        if (!even) {
                            g2 = mulScalar.apply(g2, msp);
                        }
                        newP.add(g2);
                    }
                }
                p = isoFree.apply(newP, null);
                newP.clear();
                if (multibody && doMultiFromPair) {
                    for (Graph g : trueMultiP) {
                        Set<Graph> gSet = splitOneBC.apply(g, splitOneParameters);
                        for (Graph g2 : gSet) {
                            boolean even = true;
                            for (Edge e : g2.edges()) {
                                if (e.getColor() == nfBond) {
                                    even = !even;
                                    e.setColor(fBond);
                                }
                            }
                            if (!even) {
                                g2 = mulScalar.apply(g2, msp);
                            }
                            newP.add(g2);
                        }
                    }
                    trueMultiP = isoFree.apply(newP, null);
                }
//                System.out.println("isofreeing on "+newP.size()+" Ree-Hooverish diagrams (from "+p.size()+" f-diagrams)");
            }

            FactorOnce factorOnce = new FactorOnce();
            FactorOnceParameters fop = new FactorOnceParameters((byte)0, false);
            newP.clear();
            newP.addAll(maxIsomorph.apply(p, mip));
            p.clear();
            p.addAll(newP);

            if (doDisconnectedMatching) {
                // match up singly-connected (in p) with disconnected diagrams.
                // we have to do this last so that our cancelMap remains valid.
                newP.clear();
                msp = new MulScalarParameters(-1, 1);

                if (multibody && doMultiFromPair) {
                    // first make passes through p to match up diagrams that are
                    // singly connected at a point with multi-bonds
                    happyArticulation = new ArticulatedAt0(doExchange, mmBond);
                    mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), happyArticulation);
                    Set<Graph> fullyDisconnectedMultiP = new HashSet<Graph>();
                    for (byte i=0; i<n-1; i++) {
                        for (Graph g : p) {
                            newP.add(g);
                            if (!graphHasEdgeColor(g, mmBond)) {
                                continue;
                            }
                            if (NumRootNodes.value(g) == i && hap.check(g)) {
                                // articulated and has (i+1) components
                                g = maxIsomorph.apply(g, mip);
                                boolean multiArticulated = false;
                                for (byte j=0; j<g.getOutDegree((byte)0); j++) {
                                    if (g.getEdge((byte)0,g.getOutNode((byte)0,j)).getColor() == mmBond) {
                                        multiArticulated = true;
                                        break;
                                    }
                                }
                                if (!multiArticulated) {
                                    fullyDisconnectedMultiP.add(g);
                                    continue;
                                }
                                Set<Graph> gfSet = factorOnce.apply(g, fop);
                                Graph gf = gfSet.iterator().next();
                                g = mulScalar.apply(g, new MulScalarParameters(-1,1));
                                newP.add(gf);
                                cancelMap.put(gf, g);
                                if (!hap.check(gf)) {
                                    // fully disconnected; no more articulation points.
                                    fullyDisconnectedMultiP.add(gf);
                                }
                            }
                        }
                        p = makeGraphList();
                        p.addAll(newP);
                        newP.clear();
                    }
                    // the not-cancelled disconnected parts of multi-body P are now in fullyDisconnectedP
                    // but we actually want our P expression to look like what we have in trueMultiP
                    // so take disconnected diagrams from trueMultiP and subtract fullyDisconnectedP
                    // this is then the remainder which we need to add to p
                    for (Graph g : trueMultiP) {
                        if (NumRootNodes.value(g) > 0) {
                            newP.add(g);
                        }
                    }
                    newP.addAll(mulScalar.apply(fullyDisconnectedMultiP, new MulScalarParameters(-1,1)));
                    newP = maxIsomorph.apply(isoFree.apply(newP, null), mip);
                    disconnectedP.addAll(newP);
                }
                
                for (Graph g : p) {
                    boolean ap = hap.check(g);
                    boolean con = hap.isConnected();
                    if (doMultiFromPair && graphHasEdgeColor(g, mmBond)) {
                        newP.add(g);
                        continue;
                    }
                    if (con && ap && hap.getArticulationPoints().contains((byte)0)) {
                        // newP will contain connected diagrams
                        g = g.copy();
                        newP.add(g);
                        Set<Graph> gfSet = factorOnce.apply(g, fop);
                        Graph gf = gfSet.iterator().next(); // we know we only have 1 iterate
                        disconnectedP.add(gf);
                        gf = mulScalar.apply(gf, msp);
                        cancelMap.put(g,gf);
                    }
                    else if (con) {
                        // this is a biconnected diagram;
                        newP.add(g.copy());
                    }
                    else if (graphHasEdgeColor(g, fBond)){
                        // this is a disconnected diagram;
                        disconnectedP.add(g.copy());
                    }
                    else {
                        newP.add(g.copy());
                    }
                }

                p = makeGraphList();
                p.addAll(newP);

                // we don't need to re-isofree p, we know that's still good.
                // some of our new disconnected diagrams might condense with the old ones
                disconnectedP = maxIsomorph.apply(isoFree.apply(disconnectedP, null), mip);
            }

            ComponentSubst compSubst = new ComponentSubst();
            List<ComponentSubstParameters> cpsList = new ArrayList<ComponentSubstParameters>();
//            ComponentSubstParameters[] csp = new ComponentSubstParameters[n+1];
            if (doMinimalBC) {
                // this only alters the diagrams in disconnectedP
                
                // group together all n-point diagrams (diagrams with n field nodes)
                // that evaluate to infinity (due to insufficient connectivity).
                // these are the diagrams that must be evaluated together during MSMC
                // at nth order.  all other n-point diagrams can be evaluated as
                // products of smaller diagrams

                Set<Graph> bcP = new HashSet<Graph>();
                IsBiconnected isBi = new IsBiconnected();
                for (Graph g : p) {
                    if (g.nodeCount() > 3 && !graphHasEdgeColor(g, mBond) && !graphHasEdgeColor(g, mmBond) && isBi.check(g) && 
                            (!doExchange || (!graphHasEdgeColor(g, excBond) && !graphHasEdgeColor(g, mxcBond)))) {
                        bcP.add(g);
                        continue;
                    }
                }
                
                for (byte i=4; i<=n; i++) {
                    HashSet<char[]> bicompSeen = new HashSet<char[]>();
                    char[] thisBicomp = null;
                    // find all multi graphs of size i without a root point
                    // then, Mi = mi1 + Mi2 + Mi3 + ...
                    // where mi1, Mi2, Mi3... are the diagrams we find here.
                    // mi1 is the fully connected diagram with mBonds, while Mi2, Mi3, etc.
                    //  will be diagrams containing smaller multibody components composed of MBonds
                    while (true) {
                        Coefficient mcoef = null;
                        Graph gbc = null;
                        Set<Graph> bcSubst = new HashSet<Graph>();
                        thisBicomp = null;
                        for (Graph g : bcP) {
                            if (g.nodeCount() == i && NumRootNodes.value(g) == 0) {
                                if (doExchange) {
                                    // sigh.  we have to keep things separate for each list of
                                    // exchange groups (AAAA, AAAB, AABB, AAAC, etc)
                                    char[] iBicomp = new char[i];
                                    for (byte j=0; j<i; j++) {
                                        iBicomp[j] = g.getNode(j).getColor();
                                    }
                                    java.util.Arrays.sort(iBicomp);
                                    if (thisBicomp == null) {
                                        boolean seen = false;
                                        for (char[] seenComp : bicompSeen) {
                                            if (java.util.Arrays.equals(seenComp,iBicomp)) {
                                                seen = true;
                                                break;
                                            }
                                        }
                                        if (seen) {
                                            // this is not the droid we are looking for
                                            continue;
                                        }
                                        thisBicomp = iBicomp;
                                        bicompSeen.add(thisBicomp);
                                    }
                                    else if (!java.util.Arrays.equals(iBicomp, thisBicomp)) {
                                        // this is not the droid we are looking for
                                        continue;
                                    }
                                }
                                        
                                g = g.copy();
                                if (graphHasEdgeColor(g, eBond)) {
                                    // this is a lower order (disconnected) diagram
                                    g.coefficient().multiply(new CoefficientImpl(-1));
                                    bcSubst.add(g);
                                }
                                else {
                                    gbc = g.copy();
                                    // this is the large diagram in this set.
                                    for (Edge e : g.edges()) {
                                        e.setColor(efbcBond);
                                    }
                                    bcSubst.add(g);
                                    mcoef = new CoefficientImpl(1);
                                    mcoef.divide(g.coefficient());
                                }
                            }
                        }
                        if (doExchange && thisBicomp == null) break;
                        bcSubst = mulScalar.apply(bcSubst, new MulScalarParameters(mcoef));
    
                        // replace all mBond groups of size i with MBond
                        // now mi1 = Mi - Mi2 - Mi3 - ...
                        // where mi1 is the fully connected diagram of size i
                        ComponentSubstParameters csp = new ComponentSubstParameters(gbc, bcSubst, mfp);
                        disconnectedP = maxIsomorph.apply(isoFree.apply(compSubst.apply(disconnectedP, csp), null), mip);
                        cpsList.add(csp);
                        if (!doExchange) break;
                    }
                }

            }

            if (doDisconnectedMatching) {
                // disconnected graphs with i-1 components
                @SuppressWarnings("unchecked")
                Set<Graph>[] newDisconnectedP = new HashSet[n+1];
                for (int i=0; i<n+1; i++) {
                    newDisconnectedP[i] = new HashSet<Graph>();
                }
                SplitGraph graphSplitter = new SplitGraph();
                for (Graph g : disconnectedP) {
                    Set<Graph> gSplit = graphSplitter.apply(g);
                    newDisconnectedP[gSplit.size()].add(g);
                }
                IsConnected isCon = new IsConnected();
                newP = makeGraphList();
                for (Graph g : p) {
                    if (graphHasEdgeColor(g, mmBond) && !isCon.check(g) && cancelMap.get(g) == null) {
                        // these are graphs that resulted from factoring multi-bond graphs
                        // they didn't make it into disconnectedP, but we now want to treat them
                        // the same as disconnectedP (continue factoring them)
                        Set<Graph> gSplit = graphSplitter.apply(g);
                        newDisconnectedP[gSplit.size()].add(g);
                    }
                    else {
                        newP.add(g);
                    }
                }
                p = newP;
                disconnectedP.clear();
                for (int i = 0; i<n; i++) {
                    // looking for graphs with i components
                    for (Graph g : newDisconnectedP[i]) {
                        boolean ap = hap.check(g);
                        if (ap && hap.getArticulationPoints().contains((byte)0)) {
                            // this is a disconnected diagram with a singly-connected component.
                            // we need to match it up with a more disconnected diagram
                            // all of the articulation points we're interested in will be at 0
                            // (because of happyArticulation.  Other articulation points will be
                            // from exchange groups, but we can't split those.
                            Set<Graph> gfSet = factorOnce.apply(g, fop);
                            Graph gf = gfSet.iterator().next();
                            gf = mulScalar.apply(gf, msp);
                            cancelMap.put(g,gf);

                            if (doMinimalBC && gf.nodeCount() > 5) {
                                // we've made a new diagram and we need to try to make our BC substitutions
                                for (ComponentSubstParameters csp : cpsList) {
                                    gfSet = compSubst.apply(gfSet, csp);
                                }
                            }
                            newDisconnectedP[i+1].addAll(gfSet);
                        }
                    }
    
                    newDisconnectedP[i+1] = maxIsomorph.apply(isoFree.apply(newDisconnectedP[i+1], null), mip);
    
                    disconnectedP.addAll(newDisconnectedP[i]);
                }
            }
        }

        
        if (isInteractive) {
            topSet.clear();
            topSet.addAll(p);
            topSet.addAll(disconnectedP);
            System.out.println("\nP");
            System.out.println(topSet.size()+" graphs");
            for (Graph g : topSet) {
                System.out.println(g);
                if (cancelMap != null) {
                    Graph cancelGraph = cancelMap.get(g);
                    if (cancelGraph != null) {
                        System.out.println("    "+cancelGraph);
                    }
                }
            }
            ClusterViewer.createView("P", topSet);
        }


        GraphList pFinal = makeGraphList();
        pFinal.addAll(p);
        p = pFinal;
        GraphList disconnectedPFinal = makeGraphList();
        disconnectedPFinal.addAll(disconnectedP);
        disconnectedP = disconnectedPFinal;

    }

    public static final class ArticulatedAt0 implements Property {
        protected final HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        protected final boolean doExchange;
        protected final char mBond;
        public ArticulatedAt0(boolean doExchange) {
            this(doExchange, '0');
        }
        public ArticulatedAt0(boolean doExchange, char mBond) {
            this.doExchange = doExchange;
            this.mBond = mBond;
        }
        public boolean check(Graph graph) {
            if (doExchange) {
                // exchange groups must be in order, smallest to largest
                char prevColor = graph.getNode((byte)0).getColor();
                if (prevColor >= 'a' && prevColor <= 'z') prevColor += 'A' - 'a';
                for (byte i=1; i<graph.nodeCount(); i++) {
                    char color = graph.getNode(i).getColor();
                    if (color >= 'a' && color <= 'z') color += 'A' - 'a';
                    if (color < prevColor) return false;
                    prevColor = color;
                }
            }
            List<List<Byte>> comps = CVisitor.getComponents(graph);
            if (comps.size() > 1 && !doExchange) {
                // disconnected diagram.  require contiguous components
                // exchange checks above might require that we have a not-contiguous component
                byte maxNode = -1;
                for (List<Byte> iComp : comps) {
                    maxNode += iComp.size();
                    for (byte iNode : iComp) {
                        if (iNode > maxNode) return false;
                    }
                }
            }
            boolean ap = hap.check(graph);
//            boolean con = hap.isConnected();
            if (!ap) return true;
            if (doExchange) {
                // it's OK if we're only articulated at exchange points
                boolean articulatedAtA = false;
                for (byte node : hap.getArticulationPoints()) {
                    if (node == 0) return true;
                    if (graph.getNode(node).getColor() == 'A') articulatedAtA = true;
                }
                if (!articulatedAtA) return true;
            }
            if (mBond != '0') {
                // we're OK if
                // 1. not articulated (not satisfied)
                // or
                // 2. articulation point at 0
                // and
                // 2a. 0 has an mBond
                // or
                // 2b. no articulation point has an mBond
                boolean multiArticulation = false;
                boolean articulatedAt0 = false;
                for (byte node : hap.getArticulationPoints()) {
                    if (node == 0) {
                        articulatedAt0 = true;
                    }
                    for (byte i=0; i<graph.getOutDegree(node); i++) {
                        byte j = graph.getOutNode(node, i);
                        if (graph.getEdge(node,j).getColor() == mBond) {
                            if (node == 0) return true; // satisfied 2a
                            multiArticulation = true;
                        }
                    }
                }
                // 
                return articulatedAt0 && !multiArticulation;
            }
            // 0 must be an articulation point
            for (byte node : hap.getArticulationPoints()) {
                if (node == 0) return true;
            }
            return false;
        }
    }

    public static class ComparatorNumFieldNodesExchange implements Comparator<Graph> {
        public int compare(Graph g1, Graph g2) {
            int fieldCount1 = NumFieldNodes.value(g1);
            int fieldCountEx1 = fieldCount1;
            int multi1 = 0;
            for (Node node : g1.nodes()) {
                char c = node.getColor();
                if (c >= 'A' && c <= 'Z') {
                    if (c >= 'M') {
                        c -= 'M' - 'A';
                        multi1++;
                    }
                    fieldCountEx1 += c - 'A';
                }
                else {
                    if (c >= 'm') {
                        c -= 'm' - 'a';
                        multi1++;
                    }
                    fieldCountEx1 += c - 'b';
                }
            }
            int fieldCount2 = NumFieldNodes.value(g2);
            int fieldCountEx2 = fieldCount2;
            int multi2 = 0;
            for (Node node : g2.nodes()) {
                char c = node.getColor();
                if (c >= 'A' && c <= 'Z') {
                    if (c >= 'M') {
                        c -= 'M' - 'A';
                        multi2++;
                    }
                    fieldCountEx2 += c - 'A';
                }
                else {
                    if (c >= 'm') {
                        multi2++;
                        c -= 'm' - 'a';
                    }
                    fieldCountEx2 += c - 'b';
                }
            }
            if (fieldCountEx1 != fieldCountEx2) return fieldCountEx1 - fieldCountEx2;
            int d = (fieldCountEx1-fieldCount1) - (fieldCountEx2 - fieldCount2);
            if (d != 0) return d;
            return multi1 - multi2;
        }
    }

    public static boolean graphHasEdgeColor(Graph g, char color) {
        for (Edge edge : g.edges()) {
            if (edge.getColor() == color) {
                return true;
            }
        }
        return false;
    }

    public static int getGroupID(int[] ids, int n) {
        int idl = ids.length;
        if (idl == 1) return ids[0];

        while (groupStart.length <= idl) {
            groupStart = (int[][][])Arrays.addObject(groupStart, new int[0][0]);
        }
        while (groupStart[idl].length <= n) {
            groupStart[idl] = (int[][])Arrays.addObject(groupStart[idl], new int[0]);
        }
        while (myIds.length <= idl) {
            myIds = (int[][])Arrays.addObject(myIds, new int[myIds.length]);
        }
        
        if (groupStart[idl][n].length == 0) {
            int[] nGroupStart = new int[n-idl+1];
            int nGroups = 0;
            int num = n-1;
            int den1 = n-1-(idl-1);
            int den2 = idl-1;
            int g = (int)(SpecialFunctions.factorial(num)/(SpecialFunctions.factorial(den1)*SpecialFunctions.factorial(den2)));
            for (int i=0; i<n-(idl-1); i++) {
                nGroupStart[i] = nGroups;
                nGroups += g;
                g *= (den1-i);
                g /= (num-i);
            }
            groupStart[idl][n] = nGroupStart;
        }
        for (int j=0; j<idl-1; j++) {
            myIds[idl-1][j] = ids[j+1]-ids[0]-1;
        }
        return groupStart[idl][n][ids[0]] + getGroupID(myIds[idl-1], n-ids[0]-1);
    }

    public static int singletId(int id0, int n) {
        return id0;
    }

    /**
     * Returns the pair ID for the given indices (id0, id1, id2)
     */
    public static int pairId(int id0, int id1, int n) {
        return (2*n-id0-1)*id0/2 + (id1-id0-1);
    }
    
    /**
     * Returns the triplet ID for the given indices (id0, id1, id2) and given
     * number of molecules (n).  Triplets are ordered as
     * (0,1,2), (0,1,3), (0,1,4)... (0,1,n-1)
     * (0,2,3), (0,2,4), (0,2,5)... (0,2,n-1)
     * ...
     * (0,n-3,n-2),(0,n-3,n-1)
     * (0,n-2,n-1)
     * 
     * Indices must be given in order.  If idx < 0 or idx >= n, the id returned
     * will be nonsense (or there will be an exception).
     */
    public static int tripletId(int id0, int id1, int id2, int n) {
        while (tripletStart.length <= n) {
            tripletStart = (int[][])Arrays.addObject(tripletStart, new int[0]);
        }
        if (tripletStart[n].length == 0) {
            int[] nTripletStart = new int[n-2];
            int nTriplets = 0;
            for (int i=1; i<n-1; i++) {
                nTripletStart[i-1] = nTriplets;
                nTriplets += (n-i)*(n-i-1)/2;
            }
            tripletStart[n] = nTripletStart;
        }
        return tripletStart[n][id0] + (2*n-id0-id1-2)*(id1-id0-1)/2 + (id2-id1-1);
    }

    /**
     * Returns the quad ID for the given indices (id0, id1, id2, id3) and given
     * number of molecules (n).  Quads are ordered the same as triplets.
     */
    public static int quadId(int id0, int id1, int id2, int id3, int n) {
        while (quadStart.length <= n) {
            quadStart = (int[][])Arrays.addObject(quadStart, new int[0]);
        }
        if (quadStart[n].length == 0) {
            int[] nQuadStart = new int[n-3];
            int nQuads = 0;
            for (int i=1; i<n-2; i++) {
                nQuadStart[i-1] = nQuads;
                nQuads += tripletId(n-i-3, n-i-2, n-i-1, n-i)+1;
            }
            quadStart[n] = nQuadStart;
        }
        return quadStart[n][id0] + tripletId(id1-id0-1, id2-id0-1, id3-id0-1, n-id0-1);
    }

    /**
     * Returns the quint ID for the given indices (id0, id1, id2, id3, id4) and given
     * number of molecules (n).  Quints are ordered the same as triplets.
     */
    public static int quintId(int id0, int id1, int id2, int id3, int id4, int n) {
        while (quintStart.length <= n) {
            quintStart = (int[][])Arrays.addObject(quintStart, new int[0]);
        }
        if (quintStart[n].length == 0) {
            int[] nQuintStart = new int[n-4];
            int nQuints = 0;
            for (int i=1; i<n-3; i++) {
                nQuintStart[i-1] = nQuints;
                nQuints += quadId(n-i-4, n-i-3, n-i-2, n-i-1, n-i)+1;
            }
            quintStart[n] = nQuintStart;
        }
        return quintStart[n][id0] + quadId(id1-id0-1, id2-id0-1, id3-id0-1, id4-id0-1, n-id0-1);
    }

    /**
     * Returns the six-molecule ID for the given indices (id0, id1, id2, id3, id4, id5)
     * and given number of molecules (n).  Sixes are ordered the same as triplets.
     */
    public static int sixId(int id0, int id1, int id2, int id3, int id4, int id5, int n) {
        while (sixStart.length <= n) {
            sixStart = (int[][])Arrays.addObject(sixStart, new int[0]);
        }
        if (sixStart[n].length == 0) {
            int[] nSixStart = new int[n-5];
            int nSixes = 0;
            for (int i=1; i<n-4; i++) {
                nSixStart[i-1] = nSixes;
                nSixes += quintId(n-i-5, n-i-4, n-i-3, n-i-2, n-i-1, n-i)+1;
            }
            sixStart[n] = nSixStart;
        }
        return sixStart[n][id0] + quintId(id1-id0-1, id2-id0-1, id3-id0-1, id4-id0-1, id5-id0-1, n-id0-1);
    }
}
