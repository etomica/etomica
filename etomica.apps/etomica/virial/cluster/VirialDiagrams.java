package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IdenticalGraphFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.AllIsomorphs;
import etomica.graph.operations.AllIsomorphs.AllIsomorphsParameters;
import etomica.graph.operations.ComponentSubst;
import etomica.graph.operations.ComponentSubst.ComponentSubstParameters;
import etomica.graph.operations.Decorate;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByConstant;
import etomica.graph.operations.DifByConstant.DifByConstantParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.Factor.BCVisitor;
import etomica.graph.operations.FactorOnce;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.PCopy;
import etomica.graph.operations.ReeHoover;
import etomica.graph.operations.ReeHoover.ReeHooverParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.SplitOne;
import etomica.graph.operations.SplitGraph.CVisitor;
import etomica.graph.operations.SplitOne.SplitOneParameters;
import etomica.graph.operations.SplitOneBiconnected;
import etomica.graph.operations.SplitParameters;
import etomica.graph.operations.Unfactor;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.property.NumRootNodes;
import etomica.graph.property.Property;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.DepthFirst;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.util.Arrays;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterBondsNonAdditive;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumMultibody;
import etomica.virial.ClusterSumShell;
import etomica.virial.MayerFunction;
import etomica.virial.MayerFunctionNonAdditive;

public class VirialDiagrams {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Set<Graph> multiP;
    protected Set<Graph> rho;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Graph> cancelMap;
    protected boolean doShortcut;
    protected boolean doMinimalMulti;
    protected boolean doMinimalBC;
    protected boolean doKeepEBonds;
    protected boolean doExchange;
    protected boolean doExchangeNoF;
    protected final char nodeColor = Metadata.COLOR_CODE_0;
    protected char[] flexColors;
    public char fBond, eBond, excBond, mBond, MBond, efbcBond, ffBond, mxcBond;

    protected static int[][] tripletStart = new int[0][0];
    protected static int[][] quadStart = new int[0][0];
    protected static int[][] quintStart = new int[0][0];
    protected static int[][] sixStart = new int[0][0];

    public static void main(String[] args) {
        final int n = 4;
        boolean multibody = true;
        boolean flex = true;
        boolean doKeepEBonds = false;
        boolean doReeHoover = true;
        boolean doExchange = false;
        VirialDiagrams virialDiagrams = new VirialDiagrams(n, multibody, flex, true);
        virialDiagrams.setDoReeHoover(doReeHoover);
        virialDiagrams.setDoKeepEBonds(doKeepEBonds);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.setDoExchange(doExchange);
        if (doExchange) {
            virialDiagrams.setDoExchangeNoF(false);
        }
        if (multibody) {
            virialDiagrams.setDoMinimalMulti(true);
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
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
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
    }

    public void setDoExchange(boolean newDoExchange) {
        if (!doKeepEBonds && newDoExchange) {
            throw new RuntimeException("exchange only works with keep-eBonds");
        }
        doExchange = newDoExchange;
    }

    public void setDoExchangeNoF(boolean newDoExchangeNoF) {
        if (!doExchange) {
            throw new RuntimeException("this only makes sense with exchange on");
        }
        doExchangeNoF = newDoExchangeNoF;
    }

    public void setDoMinimalMulti(boolean newDoMinimalMulti) {
        if (!multibody) {
            throw new RuntimeException("can't set minimalMulti without multi");
        }
        doMinimalMulti = newDoMinimalMulti;
    }

    public void setDoMinimalBC(boolean newDoMinimalBC) {
        if (!multibody && !flex) {
            throw new RuntimeException("you need multi or flex to do minimal bc");
        }
        doMinimalBC = newDoMinimalBC;
    }

    public Set<Graph> getVirialGraphs() {
        if (p == null) {
            makeVirialDiagrams();
        }
        return p;
    }

    public Set<Graph> getMSMCGraphs(boolean connectedOnly, boolean getMultiGraphs) {
        if (p == null) {
            makeVirialDiagrams();
        }
        GraphList<Graph> allP = makeGraphList();
        if (getMultiGraphs) {
            if (!multibody) throw new RuntimeException("oops");
            if (doMinimalMulti) {
                // we already constructed the exact set we want
                allP.addAll(multiP);
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
                if (!multibody || (doMinimalMulti && !graphHasEdgeColor(g, MBond))
                               || (!doMinimalMulti && !graphHasEdgeColor(g, mBond))) {
                    allP.add(g);
                    if (!connectedOnly) {
                        Graph c = cancelMap.get(g);
                        if (c != null) {
                            allP.add(c);
                        }
                    }
                }
            }
        }
        GraphList<Graph> pn = makeGraphList();
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
        return makeVirialCluster(f, null);
    }
    public ClusterSum makeVirialCluster(MayerFunction f, MayerFunctionNonAdditive fMulti) {
        if (p == null) {
            makeVirialDiagrams();
        }
        boolean doMulti = fMulti != null;
        if (!multibody && doMulti) {
            throw new RuntimeException("can't make multi-body bonds without multi-body diagrams");
        }
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        Set<Graph> pn = getMSMCGraphs(false, fMulti!=null);
        for (Graph g : pn) {
            int nDiagrams = populateEFBonds(g, allBonds, false, !doMulti);
            if (nDiagrams == 0) continue;
            if (flex && !doMulti) {
                populateEFBonds(g, allBonds, true, true);
                nDiagrams *= 2;
            }
            double w = g.coefficient().getValue()/nDiagrams;
            for (int i=0; i<nDiagrams; i++) {
                weights.add(w);
            }
            if (weights.size() != allBonds.size()) {
                throw new RuntimeException("oops");
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
    
    public int populateEFBonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap, boolean pairOnly) {
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
        boolean multiGraph = graphHasEdgeColor(g, mBond);
        if (multiGraph && pairOnly) return 0;
        if (!pairOnly && !multiGraph) return 0;
        int rv = 0;
        if (multiGraph) {
            // multibody graph.  we need to generate all permutations in order
            // to get proper cancellation
            AllIsomorphs allIso = new AllIsomorphs();
            AllIsomorphsParameters allIsoParams = new AllIsomorphsParameters(true);
            Set<Graph> permutations = allIso.apply(g, allIsoParams);
            List<List<Byte>> multiBiComponents = new ArrayList<List<Byte>>();
            int nPoints = g.nodeCount();

            rv = permutations.size();
            for (Graph gp : permutations) {
                ArrayList<int[]> fbonds = new ArrayList<int[]>();
                int[][] mBonds = new int[nPoints+1][0];
                multiBiComponents.clear();
                BCVisitor bcv = new BCVisitor(multiBiComponents);
                new Biconnected().traverseAll(gp, bcv);
                for (List<Byte> comp : multiBiComponents) {
                    if (comp.size() == 1) continue;
                    if (!gp.hasEdge(comp.get(0), comp.get(1)) || gp.getEdge(comp.get(0), comp.get(1)).getColor() != mBond) {
                        // we're doing multi graphs, but encountered f-bonds
                        // find all f-bonds within this component
                        for (int id1 = 0; id1<comp.size()-1; id1++) {
                            byte n1 = comp.get(id1);
                            for (int id2 = id1+1; id2<comp.size(); id2++) {
                                byte n2 = comp.get(id2);
                                if (!gp.hasEdge(n1, n2)) continue;
                                if (gp.getEdge(n1, n2).getColor() != fBond) {
                                    throw new RuntimeException("oops");
                                }
                                fbonds.add(new int[]{n1, n2});
                            }
                        }
                        continue;
                    }
                    int groupID = -1;
                    // determine the groupID for this component
                    byte id0 = comp.get(0);
                    byte id1 = comp.get(1);
                    byte id2 = comp.get(2);
                    int size = comp.size();
                    if (comp.size() == 3) {
                        groupID = tripletId(id0, id1, id2, nPoints);
                    }
                    else if (comp.size() == 4) {
                        byte id3 = comp.get(3);
                        groupID = quadId(id0, id1, id2, id3, nPoints);
                    }
                    else if (comp.size() == 5) {
                        byte id3 = comp.get(3);
                        byte id4 = comp.get(4);
                        groupID = quintId(id0, id1, id2, id3, id4, nPoints);
                    }
                    else if (comp.size() == 6) {
                        byte id3 = comp.get(3);
                        byte id4 = comp.get(4);
                        byte id5 = comp.get(5);
                        groupID = sixId(id0, id1, id2, id3, id4, id5, nPoints);
                    }
                    int[] newGroups = new int[mBonds[size].length+1];
                    System.arraycopy(mBonds[size], 0, newGroups, 0, mBonds[size].length);
                    newGroups[newGroups.length-1] = groupID;
                    mBonds[size] = newGroups;
                }
                allBonds.add(new ClusterBondsNonAdditive(nPoints, new int[][][]{fbonds.toArray(new int[0][0])}, mBonds));
            }
        }
        else {
            ArrayList<int[]> fbonds = new ArrayList<int[]>();
            rv = 1;
            for (Node node1 : g.nodes()) {
                for (Node node2 : g.nodes()) {
                    if (node1.getId() >= node2.getId()) continue;
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
                        else {
                            throw new RuntimeException("oops, unknown bond "+edgeColor);
                        }
                    }
                }
            }
            if (ebonds.size() > 0) {
                allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0])}));
            }
            else {
                allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0])}));
            }
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
        ArrayList<Double> weights = new ArrayList<Double>();
        if (doMinimalBC) {
            // group all biconnected pairwise graphs into a single clusterSum
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            double leadingCoef = 0;
            for (Graph g : pn) {
                if (graphHasEdgeColor(g, mBond) || !isBi.check(g)) continue;
                populateEFBonds(g, allBonds, false, true);
                double coef = g.coefficient().getValue();
                if (leadingCoef == 0) {
                    leadingCoef = coef;
                    coef = 1;
                }
                else {
                    coef /= leadingCoef;
                }
                if (flex) {
                    populateEFBonds(g, allBonds, true, true);
                    weights.add(0.5*coef);
                    weights.add(0.5*coef);
                }
                else {
                    weights.add(coef);
                }
            }
            double[] w = new double[weights.size()];
            for (int i=0; i<w.length; i++) {
                w[i] = weights.get(i);
            }
            if (n > 3 && !flex) {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f,e}));
            }
            else {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}));
            }
        }
        double[] w1 = new double[]{1};
        if (flex) {
            w1 = new double[]{0.5,0.5};
        }
        for (Graph g : pn) {
            if (graphHasEdgeColor(g, mBond)) continue;
            if (doMinimalBC && isBi.check(g)) continue;
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            populateEFBonds(g, allBonds, false, true);
            double[] thisW = w1;
            if (flex) {
                populateEFBonds(g, allBonds, true, true);
            }
            if (flex && cancelMap.get(g) != null) {
                Graph cg = cancelMap.get(g);
                populateEFBonds(cg, allBonds, false, true);
                populateEFBonds(cg, allBonds, true, true);
                thisW = new double[2+2];
                thisW[0] = thisW[1] = 0.5;
                for (int i=2; i<thisW.length; i++) {
                    thisW[i] = -0.5/1;
                }
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
    
    public Set<Graph> getExtraDisconnectedVirialGraphs() {
        if (p == null) {
            makeVirialDiagrams();
        }
        GraphList<Graph> dpn = makeGraphList();
        for (Graph g : disconnectedP) {
            int fieldCount = 0;
            for (Node node : g.nodes()) {
              if (node.getType() == 'F') {
                fieldCount++;
              }
            }
            if (fieldCount == n) {
                dpn.add(g);
            }
        }
        return dpn;
    }

    public HashMap<Graph,Set<Graph>> getSplitDisconnectedVirialGraphs(Set<Graph> disconnectedGraphs) {
        HashMap<Graph,Set<Graph>> map = new HashMap<Graph,Set<Graph>>();
        SplitGraph splitGraph = new SplitGraph();
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        Relabel relabel = new Relabel();
        for (Graph g : disconnectedGraphs) {
            Set<Graph> gSplit = makeGraphList();
            Set<Graph> gSplit1 = splitGraph.apply(g);
            for (Graph gs : gSplit1) {
                // the graph we get from splitting might not be in our preferred bonding arrangement
                Graph gsmax = maxIsomorph.apply(gs, MaxIsomorph.PARAM_ALL);
                HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
                if (hap.check(gsmax)) {
                    if (!hap.getArticulationPoints().contains(0)) {
                        // the graph is articulated but not with the articulation point at 0.
                        // we calculate the diagram with the articulation point at 0, so permute
                        // now so that we can easily identify which graph this is
                        byte[] permutations = new byte[gsmax.nodeCount()];
                        for (int i=0; i<permutations.length; i++) {
                            permutations[i] = (byte)i;
                        }
                        permutations[0] = hap.getArticulationPoints().get(0);
                        permutations[hap.getArticulationPoints().get(0)] = 0;
                        RelabelParameters rp = new RelabelParameters(permutations);
                        gsmax = relabel.apply(gsmax, rp);
                    }
                }
                gSplit.add(gsmax);
            }
            map.put(g, gSplit);
        }
        return map;
    }

    public GraphList<Graph> makeGraphList() {
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList<Graph> graphList = new GraphList<Graph>(comp);
        return graphList;
    }        
    
    public void makeRhoDiagrams() {
        flexColors = new char[0];
        if (flex) {
           flexColors = new char[]{nodeColor};
        }

        final HashMap<Character,Integer> colorOrderMap = new HashMap<Character,Integer>();
        if (MetadataImpl.metaDataComparator == null) {
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
        }

        GraphList<Graph> topSet = makeGraphList();

        char oneBond = 'o';
        fBond = 'f';
        eBond = 'e';
        mBond = 'm';  // multi-body
        MBond = 'M';  // Multi-body
        efbcBond = 'b';
        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(mBond, 1);
        colorOrderMap.put(efbcBond, 2);
        colorOrderMap.put(fBond, 3);
        colorOrderMap.put(eBond, 4);
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = new MulFlexibleParameters(flexColors, (byte)n);
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

        }
        Metadata.COLOR_MAP.put(eBond, "red");
        Metadata.COLOR_MAP.put(fBond, "green");
        Metadata.COLOR_MAP.put(mBond, "blue");
        Metadata.COLOR_MAP.put(MBond, "orange");
        Metadata.COLOR_MAP.put(efbcBond, "fuchsia");
        Metadata.COLOR_MAP.put(excBond, "red");
        Metadata.COLOR_MAP.put(ffBond, "green");
        Metadata.COLOR_MAP.put(mxcBond, "blue");
        Metadata.DASH_MAP.put(excBond, 3);
        Metadata.DASH_MAP.put(ffBond, 3);
        Metadata.DASH_MAP.put(mxcBond, 3);

        if (isInteractive) {
            topSet.clear();
            topSet.addAll(lnfXi);
            System.out.println("\nlnfXi");
            for (Graph g : topSet) {
                System.out.println(g);
            }
//            ClusterViewer.createView("lnfXi", topSet);
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

        GraphList<Graph> topSet = makeGraphList();

        char oneBond = 'o';
        fBond = 'f';
        eBond = 'e';
        mBond = 'm';  // multi-body
        MBond = 'M';  // Multi-body
        mxcBond = 'n';  // multi-body exchange (yes, of course 'n' is bad)
        efbcBond = 'b';
        excBond = 'x';
        ffBond = 'F';
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = new MulFlexibleParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(mBond, 1);
        colorOrderMap.put(MBond, 2);
        colorOrderMap.put(eBond, 3);
        colorOrderMap.put(fBond, 4);
        colorOrderMap.put(efbcBond, 5);
        colorOrderMap.put(ffBond, 6);
        colorOrderMap.put(excBond, 7);
        colorOrderMap.put(mxcBond, 8);

        if (doShortcut && !multibody && !flex) {

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
    
            MulFlexibleParameters mfpnm1 = new MulFlexibleParameters(flexColors, (byte)(n-1));
            z.addAll(allRho[1]);
            for (int i=2; i<n+1; i++) {
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
            
            // clear these out -- we don't need them and (in extreme cases) we might need the memory
            lnfXi.clear();
            rho.clear();
            z.clear();
        }

        

        if (doExchange) {
            System.out.println("\nPe");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pe", topSet);
            // we substitute excBonds (exchange bonds, between atoms involved the exchange)
            // if exchanging atoms are connected to other atoms (with e-Bonds), we change those
            // to ff-Bonds (what a horrible name).  Our Graph algebra can't convert the e-Bonds
            // to f-bonds (because f = e*e-1) between a exchanging double-ring and a single atom,
            // so we do it manually here.  The ff-bonds are those new f bonds. 
            HashSet<Graph> pxc = new HashSet<Graph>();
            char efBond = doExchangeNoF ? eBond : ffBond;
            for (Graph g : p) {
                List<List<Byte>> components = new ArrayList<List<Byte>>();
                CVisitor v = new CVisitor(components);
                new DepthFirst().traverseAll(g, v);
                HashSet<Graph> gPermuted = new HashSet<Graph>();
                gPermuted.add(g.copy());
                for (int i=0; i<components.size(); i++) {
                    List<Byte> comp = components.get(i);
                    if (comp.size() == 1) continue;
                    HashSet<Graph> gPermutedNew = new HashSet<Graph>();
                    for (Graph gp : gPermuted) {
                        gPermutedNew.add(gp.copy());
                        if (comp.size() == 2) {
                            Graph gCopy = gp.copy();
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gPermutedNew.add(gCopy);
                        }
                        else if (comp.size() == 3) {
                            if (multibody) {
                                Graph gCopy = gp.copy();
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mBond);
                                gCopy.getEdge(comp.get(0), comp.get(2)).setColor(mBond);
                                gCopy.getEdge(comp.get(1), comp.get(2)).setColor(mBond);
                                gPermutedNew.add(gCopy);
                            }

                            Graph gCopy = gp.copy();
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gCopy.getEdge(comp.get(0), comp.get(2)).setColor(efBond);
                            gCopy.getEdge(comp.get(1), comp.get(2)).setColor(efBond);
                            gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(3,1)));

                            if (multibody) {
                                gCopy = gCopy.copy();
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(0), comp.get(2)).setColor(mBond);
                                gCopy.getEdge(comp.get(1), comp.get(2)).setColor(mBond);
                                gCopy = mulScalar.apply(gCopy, new MulScalarParameters(3,1));
                                gPermutedNew.add(gCopy);
                            }

                            if (!doExchangeNoF) {
                                gCopy = gp.copy();
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                                gCopy.deleteEdge(comp.get(0), comp.get(2));
                                gCopy.deleteEdge(comp.get(1), comp.get(2));
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(3,1)));
                            }

                            gCopy = gp.copy();
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gCopy.getEdge(comp.get(0), comp.get(2)).setColor(excBond);
                            gCopy.getEdge(comp.get(1), comp.get(2)).setColor(excBond);
                            gCopy = mulScalar.apply(gCopy, new MulScalarParameters(2,1));
                            gPermutedNew.add(gCopy);
                            if (multibody) {
                                gCopy = gCopy.copy();
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(0), comp.get(2)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(1), comp.get(2)).setColor(mxcBond);
                                gCopy = mulScalar.apply(gCopy, new MulScalarParameters(2,1));
                                gPermutedNew.add(gCopy);
                            }
                        }
                        else if (comp.size() == 4) {
                            if (gp.getEdge(comp.get(0), comp.get(1)).getColor() == mBond) continue;

                            //2
                            Graph gCopy = gp.copy();
                            for (byte j=0; j<3; j++) {
                                for (byte k=(byte)(j+1); k<4; k++) {
                                    gCopy.getEdge(comp.get(j), comp.get(k)).setColor(efBond);
                                }
                            }
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gCopy.getEdge(comp.get(2), comp.get(3)).setColor(eBond);
                            gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(6,1)));
                            if (multibody) {
                                gCopy = gp.copy();
                                for (byte j=0; j<3; j++) {
                                    for (byte k=(byte)(j+1); k<4; k++) {
                                        gCopy.getEdge(comp.get(j), comp.get(k)).setColor(mBond);
                                    }
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mxcBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(6,1)));
                            }
                            if (!doExchangeNoF) {
                                gCopy = gp.copy();
                                for (byte j=0; j<2; j++) {
                                    gCopy.deleteEdge(comp.get(j), comp.get(2));
                                    gCopy.getEdge(comp.get(j), comp.get(3)).setColor(ffBond);
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(12,1)));
                                gCopy = gp.copy();
                                for (byte j=0; j<2; j++) {
                                    for (byte k=2; k<4; k++) {
                                        gCopy.deleteEdge(comp.get(j), comp.get(k));
                                    }
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(6,1)));
                            }

                            //2+2
                            gCopy = gp.copy();
                            for (byte j=0; j<3; j++) {
                                for (byte k=(byte)(j+1); k<4; k++) {
                                    gCopy.getEdge(comp.get(j), comp.get(k)).setColor(efBond);
                                }
                            }
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gCopy.getEdge(comp.get(2), comp.get(3)).setColor(excBond);
                            gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(3,1)));
                            if (multibody) {
                                gCopy = gp.copy();
                                for (byte j=0; j<3; j++) {
                                    for (byte k=(byte)(j+1); k<4; k++) {
                                        gCopy.getEdge(comp.get(j), comp.get(k)).setColor(mBond);
                                    }
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(2), comp.get(3)).setColor(mxcBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(3,1)));
                            }
                            if (!doExchangeNoF) {
                                gCopy = gp.copy();
                                for (byte j=0; j<2; j++) {
                                    for (byte k=2; k<4; k++) {
                                        gCopy.deleteEdge(comp.get(j), comp.get(k));
                                    }
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                                gCopy.getEdge(comp.get(2), comp.get(3)).setColor(excBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(3,1)));
                            }

                            //3
                            gCopy = gp.copy();
                            for (byte j=0; j<3; j++) {
                                gCopy.getEdge(comp.get(j), comp.get(3)).setColor(efBond);
                            }
                            gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                            gCopy.getEdge(comp.get(0), comp.get(2)).setColor(excBond);
                            gCopy.getEdge(comp.get(1), comp.get(2)).setColor(excBond);
                            gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(8,1)));
                            if (multibody) {
                                gCopy = gp.copy();
                                for (byte j=0; j<3; j++) {
                                    gCopy.getEdge(comp.get(j), comp.get(3)).setColor(mBond);
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(0), comp.get(2)).setColor(mxcBond);
                                gCopy.getEdge(comp.get(1), comp.get(2)).setColor(mxcBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(8,1)));
                            }
                            if (!doExchangeNoF) {
                                gCopy = gp.copy();
                                for (byte j=0; j<3; j++) {
                                    gCopy.deleteEdge(comp.get(j), comp.get(3));
                                }
                                gCopy.getEdge(comp.get(0), comp.get(1)).setColor(excBond);
                                gCopy.getEdge(comp.get(0), comp.get(2)).setColor(excBond);
                                gCopy.getEdge(comp.get(1), comp.get(2)).setColor(excBond);
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(8,1)));
                            }

                            //4
                            gCopy = gp.copy();
                            for (byte j=0; j<3; j++) {
                                for (byte k=(byte)(j+1); k<4; k++) {
                                    gCopy.getEdge(comp.get(j), comp.get(k)).setColor(excBond);
                                }
                            }
                            gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(6,1)));
                            if (multibody) {
                                gCopy = gp.copy();
                                for (byte j=0; j<3; j++) {
                                    for (byte k=(byte)(j+1); k<4; k++) {
                                        gCopy.getEdge(comp.get(j), comp.get(k)).setColor(mxcBond);
                                    }
                                }
                                gPermutedNew.add(mulScalar.apply(gCopy, new MulScalarParameters(6,1)));
                            }
                        }
                        else {
                            throw new RuntimeException("don't have code for 5th-order yet.");
                        }
                    }
                    gPermuted = gPermutedNew;
                }
                pxc.addAll(gPermuted);
            }
            Property happyArticulation = new ArticulatedAt0();
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(happyArticulation);
            p.clear();
            p.addAll(maxIsomorph.apply(isoFree.apply(pxc, null), mip));
            
            System.out.println("\nPxc");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pxc", topSet);
        }

        if (doKeepEBonds) {
            Split split = new Split();
            SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
            Set<Graph> newP = split.apply(p, bonds);
    
            DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
            DeleteEdge deleteEdge = new DeleteEdge();
            p.clear();
            Unfactor unfactor = new Unfactor();
            Property happyArticulation = new ArticulatedAt0();
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(happyArticulation);
            p.addAll(maxIsomorph.apply(isoFree.apply(unfactor.apply(deleteEdge.apply(newP, deleteEdgeParameters), mfp), null), mip));
        }
        
        
        if (doMinimalMulti) {
            doMinimalMulti();
        }

        Set<Graph> newP = new HashSet<Graph>();

        // attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Graph>();
        cancelP = new GraphList<Graph>();
        disconnectedP = new HashSet<Graph>();
        if (!flex) {

            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Factor factor = new Factor();
            newP.clear();

            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newP.add(gf);
                }
                else {
                    newP.add(g.copy());
                }
            }
            // we've split the graphs up; now combine them if possible
            p = isoFree.apply(newP, null);
            // now recondense them (only relevant for multibody interactions)
            Unfactor unfactor = new Unfactor();
            p = isoFree.apply(unfactor.apply(p, mfp), null);
            // perform Ree-Hoover substitution (brute-force)
            if (doReeHoover) {
                if (doShortcut && !multibody) {
                    ReeHoover reeHoover = new ReeHoover();
                    p = reeHoover.apply(p, new ReeHooverParameters(eBond));
                }
                else {
                    char nfBond = 'F';
                    SplitOneParameters splitOneParameters = new SplitOneParameters(eBond, nfBond);
                    SplitOne splitOne = new SplitOne();
                    msp = new MulScalarParameters(-1, 1);
                    newP.clear();
                    for (Graph g : p) {
                        Set<Graph> gSet = splitOne.apply(g, splitOneParameters);
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
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            newP.addAll(maxIsomorph.apply(p, MaxIsomorph.PARAM_ALL));
            p = newP;
        }
        else {

            // perform Ree-Hoover substitution (brute-force)
            if (doReeHoover) {
                char nfBond = 'Z';
                SplitOneParameters splitOneParameters = new SplitOneParameters(eBond, nfBond);
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
//                System.out.println("isofreeing on "+newP.size()+" Ree-Hooverish diagrams (from "+p.size()+" f-diagrams)");
                p = isoFree.apply(newP, null);
                newP.clear();
            }

            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            FactorOnce factorOnce = new FactorOnce();
            boolean allPermutations = false;
            FactorOnceParameters fop = new FactorOnceParameters((byte)0, new char[0], allPermutations);
            newP.clear();
            Property happyArticulation = new ArticulatedAt0();
            MaxIsomorph maxIsomorph = new MaxIsomorph();
            MaxIsomorphParameters mip = new MaxIsomorphParameters(happyArticulation);
            newP.addAll(maxIsomorph.apply(p, mip));
            p.clear();
            p.addAll(newP);
            newP.clear();
            msp = new MulScalarParameters(-1, 1);
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if (con && ap) {
                    // newP will contain connected diagrams
                    g = g.copy();
                    newP.add(g);
                    Set<Graph> gfSet = factorOnce.apply(g, fop);
                    Graph gf = gfSet.iterator().next(); // we know we only have 1 iterate
                    disconnectedP.add(gf);
                    gf = mulScalar.apply(gf, msp);
                    cancelP.add(gf);
                    cancelMap.put(g,gf);
                }
                else if (con) {
                    // this is a biconnected diagram;
                    newP.add(g.copy());
                }
                else {
                    // this is a disconnected diagram;
                    disconnectedP.add(g.copy());
                }
            }
            p = newP;
            // we don't need to re-isofree p, we know that's still good.
            // some of our new disconnected diagrams might condense with the old ones
            disconnectedP = isoFree.apply(disconnectedP, null);

            Set<Graph>[] bcSubst = new HashSet[n+1];
            ComponentSubst compSubst = new ComponentSubst();
            ComponentSubstParameters[] csp = new ComponentSubstParameters[n+1];
            if (doMinimalBC) {
                // group together all n-point diagrams (diagrams with n field nodes)
                // that evaluate to infinity (due to insufficient connectivity).
                // these are the diagrams that must be evaluated together during MSMC
                // at nth order.  all other n-point diagrams can be evaluated as
                // products of smaller diagrams
                for (int i=4; i<=n; i++) {
                    bcSubst[i] = new HashSet<Graph>();
                }

                Set<Graph> bcP = new HashSet<Graph>();
                IsBiconnected isBi = new IsBiconnected();
                for (Graph g : p) {
                    if (!graphHasEdgeColor(g, mBond) && !graphHasEdgeColor(g, MBond) && isBi.check(g)) {
                        bcP.add(g);
                        continue;
                    }
                }
                MulFlexibleParameters mfpc = new MulFlexibleParameters(new char[]{nodeColor}, (byte)n);
                for (int i=4; i<=n; i++) {
                    // find all multi graphs of size i without a root point
                    // then, Mi = mi1 + Mi2 + Mi3 + ...
                    // where mi1, Mi2, Mi3... are the diagrams we find here.
                    // mi1 is the fully connected diagram with mBonds, while Mi2, Mi3, etc.
                    //  will be diagrams containing smaller multibody components composed of MBonds
                    Coefficient mcoef = null;
                    Graph gbc = null;
                    for (Graph g : bcP) {
                        if (g.nodeCount() == i && NumRootNodes.value(g) == 0) {
                            g = g.copy();
                            if (graphHasEdgeColor(g, eBond)) {
                                // this is a lower order (disconnected) diagram
                                g.coefficient().multiply(new CoefficientImpl(-1));
                                bcSubst[i].add(g);
                            }
                            else {
                                gbc = g.copy();
                                // this is the large diagram in this set.
                                for (Edge e : g.edges()) {
                                    e.setColor(efbcBond);
                                }
                                bcSubst[i].add(g);
                                mcoef = new CoefficientImpl(1);
                                mcoef.divide(g.coefficient());
                            }
                        }
                    }
                    bcSubst[i] = mulScalar.apply(bcSubst[i], new MulScalarParameters(mcoef));

                    // replace all mBond groups of size i with MBond
                    // now mi1 = Mi - Mi2 - Mi3 - ...
                    // where mi1 is the fully connected diagram of size i
                    csp[i] = new ComponentSubstParameters(gbc, bcSubst[i], mfpc);
//                    p = isoFree.apply(compSubst.apply(p, csp), null);
                    disconnectedP = isoFree.apply(compSubst.apply(disconnectedP, csp[i]), null);
                }

//                biconP = new HashSet<Graph>();
                // we have mi1 = Mi + Mi2 + Mi3
                // we want to reconstruct Mi = mi1 + mi2 + mi3
                // first we just do  Mi = mi - Mi2 - Mi3
                for (int i=4; i<n+1; i++) {
                    if (true) break;
                    Graph gm = null;
                    for (Graph g : bcSubst[i]) {
                        if (g.edgeCount() == i*(i-1)/2) {
                            // this was Mi, turn it back into mi1
                            gm = g.copy(); // remember this so we can substitute for it elsewhere
                            // invert sign here because we'll uninvert below
                            g.coefficient().multiply(new CoefficientImpl(-1));
                            for (Edge e : g.edges()) {
                                e.setColor(fBond);
                            }
                            break;
                        }
                    }
                    // flip sign of every diagram
                    bcSubst[i] = mulScalar.apply(bcSubst[i], new MulScalarParameters(new CoefficientImpl(-1)));
                    // We have Mi = mi1 + mi2 + mi3
                    // now substitute that back for j>i
                    ComponentSubstParameters cspBack = new ComponentSubstParameters(gm, bcSubst[i], mfpc);
                    for (int j=i+1; j<n+1; j++) {
                        bcSubst[j] = isoFree.apply(compSubst.apply(bcSubst[j], cspBack), null);
                    }
                    // now the leading coefficient is (1-i)/i!
//                    biconP.addAll(mulScalar.apply(bcSubst[i], new MulScalarParameters(new CoefficientImpl(1-i, (int)SpecialFunctions.factorial(i)))));
                }

//                p.clear();
//                p.addAll(nbcP);
//                p.addAll(bcP);
            }
            

            // disconnected graphs with i-1 components
            Set<Graph>[] newDisconnectedP = new HashSet[n+1];
            for (int i=0; i<n+1; i++) {
                newDisconnectedP[i] = new HashSet<Graph>();
            }
            SplitGraph graphSplitter = new SplitGraph();
            for (Graph g : disconnectedP) {
                Set<Graph> gSplit = graphSplitter.apply(g);
                newDisconnectedP[gSplit.size()].add(g);
            }
            disconnectedP.clear();
            for (int i = 0; i<n; i++) {
                // looking for graphs with i components
                for (Graph g : newDisconnectedP[i]) {
                    boolean ap = hap.check(g);
                    if (ap) {
                        byte apid = hap.getArticulationPoints().get(0);
                        Set<Graph> gfSet = factorOnce.apply(g, new FactorOnceParameters(apid, new char[0], false));
                        Graph gf = gfSet.iterator().next();
                        gf = mulScalar.apply(gf, msp);
                        cancelP.add(gf);
                        cancelMap.put(g,gf);

                        if (doMinimalBC && gf.nodeCount() > 5) {
                            // we've made a new diagram and we need to try to make our BC substitutions
                            for (int j=4; j<gf.nodeCount()-1 && j<n; j++) {
                                gfSet = compSubst.apply(gfSet, csp[j]);
                            }
                        }
                        newDisconnectedP[i+1].addAll(gfSet);
                    }
                }
                newDisconnectedP[i+1] = isoFree.apply(newDisconnectedP[i+1], null);
                disconnectedP.addAll(newDisconnectedP[i]);
            }

            Set<Graph> pNoRoot = new HashSet<Graph>();
            for (Graph g : p) {
                g = g.copy();
                for (Node node : g.nodes()) {
                    node.setType(Metadata.TYPE_NODE_FIELD);
                }
                pNoRoot.add(g);
            }
            
            // we want to condense cancelP (in case multiple diagrams were factored into the
            // same one -- is that even possible?), but want to be careful not to permute bonds.
            PCopy pcopy = new PCopy();
            IteratorWrapper wrapper = new IteratorWrapper(pcopy.apply(cancelP, null).iterator());
            GraphIterator isomorphs = new IdenticalGraphFilter(wrapper);
            cancelP = new GraphList<Graph>();
            while (isomorphs.hasNext()) {
                cancelP.add(isomorphs.next());
            }
        }

        if (isInteractive) {
            topSet.clear();
            topSet.addAll(p);
            topSet.addAll(disconnectedP);
            System.out.println("\nP");
            for (Graph g : topSet) {
                System.out.println(g);
                Graph cancelGraph = cancelMap.get(g);
                if (cancelGraph != null) {
                    System.out.println("    "+cancelGraph);
                }
            }
            ClusterViewer.createView("P", topSet);
        }


        GraphList<Graph> pFinal = makeGraphList();
        pFinal.addAll(p);
        p = pFinal;
        GraphList<Graph> disconnectedPFinal = makeGraphList();
        disconnectedPFinal.addAll(disconnectedP);
        disconnectedP = disconnectedPFinal;

    }

    protected void doMinimalMulti() {
        MulFlexibleParameters mfp = new MulFlexibleParameters(flexColors, (byte)n);
        // group together all n-point diagrams (diagrams with n field nodes)
        // that evaluate to infinity (due to insufficient connectivity).
        // these are the diagrams that must be evaluated together during MSMC
        // at nth order.  all other n-point diagrams can be evaluated as
        // products of smaller diagrams
        Factor factor = new Factor();
        MulScalar mulScalar = new MulScalar();
        IsoFree isoFree = new IsoFree();
        Set<Graph>[] multiPSubst = new HashSet[n+1];
        for (int i=3; i<=n; i++) {
            multiPSubst[i] = new HashSet<Graph>();
        }

        Set<Graph> mP = new HashSet<Graph>();
        Set<Graph> pairP = new HashSet<Graph>();
        for (Graph g : p) {
            g = factor.apply(g,mfp);
            if (!graphHasEdgeColor(g, mBond)) {
                pairP.add(g);
                continue;
            }
            mP.add(g);
        }
        MulFlexibleParameters mfpc = new MulFlexibleParameters(new char[]{nodeColor}, (byte)n);
        ComponentSubst compSubst = new ComponentSubst();
        for (int i=3; i<=n; i++) {
            // find all multi graphs of size i without a root point
            // then, Mi = mi1 + Mi2 + Mi3 + ...
            // where mi1, Mi2, Mi3... are the diagrams we find here.
            // mi1 is the fully connected diagram with mBonds, while Mi2, Mi3, etc.
            //  will be diagrams containing smaller multibody components composed of MBonds
            Coefficient mcoef = null;
            Graph gm = null;
            for (Graph g : mP) {
                if (g.nodeCount() == i && NumRootNodes.value(g) == 0) {
                    g = g.copy();
                    if (g.edgeCount() < i*(i-1)/2) {
                        // this is a lower order (disconnected) diagram
                        g.coefficient().multiply(new CoefficientImpl(-1));
                        multiPSubst[i].add(g);
                    }
                    else {
                        gm = g.copy();
                        // this is the large diagram in this set.
                        for (Edge e : g.edges()) {
                            e.setColor(MBond);
                        }
                        multiPSubst[i].add(g);
                        mcoef = new CoefficientImpl(1);
                        mcoef.divide(g.coefficient());
                    }
                }
            }
            multiPSubst[i] = mulScalar.apply(multiPSubst[i], new MulScalarParameters(mcoef));

            // replace all mBond groups of size i with MBond
            // now mi1 = Mi - Mi2 - Mi3 - ...
            // where mi1 is the fully connected diagram of size i
            ComponentSubstParameters csp = new ComponentSubstParameters(gm, multiPSubst[i], mfpc);
            mP = isoFree.apply(compSubst.apply(mP, csp), null);
        }

        multiP = new HashSet<Graph>();
        // we have mi1 = Mi + Mi2 + Mi3
        // we want to reconstruct Mi = mi1 + mi2 + mi3
        // first we just do  Mi = mi - Mi2 - Mi3
        for (int i=3; i<n+1; i++) {
            Graph gm = null;
            for (Graph g : multiPSubst[i]) {
                if (g.edgeCount() == i*(i-1)/2) {
                    // this was Mi, turn it back into mi1
                    gm = g.copy(); // remember this so we can substitute for it elsewhere
                    // invert sign here because we'll uninvert below
                    g.coefficient().multiply(new CoefficientImpl(-1));
                    for (Edge e : g.edges()) {
                        e.setColor(mBond);
                    }
                    break;
                }
            }
            // flip sign of every diagram
            multiPSubst[i] = mulScalar.apply(multiPSubst[i], new MulScalarParameters(new CoefficientImpl(-1)));
            // We have Mi = mi1 + mi2 + mi3
            // now substitute that back for j>i
            ComponentSubstParameters csp = new ComponentSubstParameters(gm, multiPSubst[i], mfpc);
            for (int j=i+1; j<n+1; j++) {
                multiPSubst[j] = isoFree.apply(compSubst.apply(multiPSubst[j], csp), null);
            }
            // now the leading coefficient is (1-i)/i!
            multiP.addAll(mulScalar.apply(multiPSubst[i], new MulScalarParameters(new CoefficientImpl(1-i, (int)SpecialFunctions.factorial(i)))));
        }

        multiP = isoFree.apply(multiP, null);
        if (isInteractive) {
            Set<Graph> topSet = makeGraphList();
            topSet.addAll(multiP);
            ClusterViewer.createView("multiP", topSet);
        }

        p.clear();
        p.addAll(pairP);
        p.addAll(mP);
    }

    public static final class ArticulatedAt0 implements Property {
        protected final IsBiconnected isBi = new IsBiconnected();
        protected final HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        public boolean check(Graph graph) {
            if (isBi.check(graph)) return true;
            boolean ap = hap.check(graph);
//            boolean con = hap.isConnected();
            if (!ap) return true;
            return hap.getArticulationPoints().contains((byte)0); 
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
            quintStart[n] = nSixStart;
        }
        return sixStart[n][id0] + quintId(id1-id0-1, id2-id0-1, id3-id0-1, id4-id0-1, id5-id0-1, n-id0-1);
    }
}
