package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.StoredIterator;
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
import etomica.graph.operations.FactorOnce;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.GraphOp;
import etomica.graph.operations.GraphOpMaxRoot;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.ReeHoover;
import etomica.graph.operations.ReeHoover.ReeHooverParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.SplitOneBiconnected;
import etomica.graph.operations.SplitOneBiconnected.SplitOneParametersBC;
import etomica.graph.operations.SplitParameters;
import etomica.graph.operations.Unfactor;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.property.NumFieldNodes;
import etomica.graph.property.NumRootNodes;
import etomica.graph.property.Property;
import etomica.graph.traversal.BCVisitor;
import etomica.graph.traversal.CVisitor;
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
import etomica.virial.cluster.CondenseExchange.CondenseExchangeParameters;
import etomica.virial.cluster.ExchangeSplit.ExchangeSplitParameters;

public class VirialDiagrams {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, disconnectedP;
    protected Set<Graph> multiP;
    protected Set<Graph> rho;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Graph> cancelMap;
    protected boolean doShortcut;
    protected boolean doMinimalMulti;
    protected boolean doMinimalBC;
    protected boolean doKeepEBonds;
    protected boolean doExchange;
    protected boolean doExchangeF;
    protected boolean doExchangeCondensing;
    protected boolean doDisconnectedMatching = true;
    protected boolean doNegativeExchange = false;
    protected final char nodeColor = Metadata.COLOR_CODE_0;
    protected char[] flexColors;
    public char fBond, eBond, excBond, mBond, MBond, efbcBond, ffBond, mxcBond, MxcBond;

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
            virialDiagrams.setDoExchangeF(true);
            virialDiagrams.setDoExchangeCondensing(true);
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
            int nDiagrams = populateEFBonds(g, allBonds, false);
            if (nDiagrams > 0) {
	            if (flex && !doMulti) {
	                populateEFBonds(g, allBonds, true);
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

    public int populateEFBonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap) {
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
        boolean multiGraph = graphHasEdgeColor(g, mBond);
        int rv = 0;
        if (multiGraph) {
            // multibody graph.  we need to generate all permutations in order
            // to get proper cancellation
            AllIsomorphs allIso = new AllIsomorphs();
            AllIsomorphsParameters allIsoParams = new AllIsomorphsParameters(true);
            Set<Graph> permutations = allIso.apply(g, allIsoParams);
            int nPoints = g.nodeCount();
            if (nPoints != n) {
                throw new RuntimeException("I would rather have these be equal");
            }

            rv = permutations.size();
            for (Graph gp : permutations) {
                ArrayList<int[]> fbonds = new ArrayList<int[]>();
                int[][] mBonds = new int[nPoints+1][0];
                List<List<Byte>> multiBiComponents = BCVisitor.getBiComponents(gp);
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
                allBonds.add(new ClusterBondsNonAdditive(flex ? nPoints+1 : nPoints, new int[][][]{fbonds.toArray(new int[0][0])}, mBonds));
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
                populateEFBonds(g, allBonds, false);
                double coef = g.coefficient().getValue();
                if (leadingCoef == 0) {
                    leadingCoef = coef;
                    coef = 1;
                }
                else {
                    coef /= leadingCoef;
                }
                if (flex) {
                    populateEFBonds(g, allBonds, true);
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
            populateEFBonds(g, allBonds, false);
            double[] thisW = w1;
            if (flex) {
                populateEFBonds(g, allBonds, true);
            }
            if (flex && cancelMap.get(g) != null) {
                Graph cg = cancelMap.get(g);
                populateEFBonds(cg, allBonds, false);
                populateEFBonds(cg, allBonds, true);
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
        comp.addComparator(new ComparatorNumFieldNodesExchange());
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
        if (doExchange) {
            // all colors are flexible
            flexColors = new char[2*n];
            for (int i=0; i<n; i++) {
                flexColors[i] = (char)('A' + i);
                flexColors[n+i] = (char)('a' + i);
            }
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
        Metadata.DASH_MAP.put(MxcBond, 3);

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

        GraphList<Graph> topSet = makeGraphList();

        char oneBond = 'o';
        fBond = 'f';
        eBond = 'e';
        mBond = 'm';  // multi-body
        MBond = 'M';  // Multi-body
        mxcBond = 'n';  // multi-body exchange (yes, of course 'n' is bad)
        MxcBond = 'N';
        efbcBond = 'b';
        excBond = 'x';
        ffBond = 'F';
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(mBond, 1);
        colorOrderMap.put(MBond, 2);
        colorOrderMap.put(eBond, 3);
        colorOrderMap.put(fBond, 4);
        colorOrderMap.put(efbcBond, 5);
        colorOrderMap.put(ffBond, 6);
        colorOrderMap.put(mxcBond, 7);
        colorOrderMap.put(MxcBond, 8);
        colorOrderMap.put(excBond, 9);

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
            
            // clear these out -- we don't need them and (in extreme cases) we might need the memory
            lnfXi.clear();
            rho.clear();
            z.clear();
        }

        

        Property happyArticulation = new ArticulatedAt0(doExchange);
        if (doExchange) {
            System.out.println("\nPe");
            topSet.clear();
            topSet.addAll(p);
            for (Graph g : topSet) {
                System.out.println(g);
            }
            ClusterViewer.createView("Pe", topSet);
            
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
        
        if (doMinimalMulti) {
            doMinimalMulti();
        }

        Set<Graph> newP = new HashSet<Graph>();

        // attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Graph>();
        disconnectedP = new HashSet<Graph>();
        if (!flex) {

            if (doDisconnectedMatching) {
                Factor factor = new Factor();

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
            newP.addAll(maxIsomorph.apply(p, MaxIsomorph.PARAM_ALL));
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
//                System.out.println("isofreeing on "+newP.size()+" Ree-Hooverish diagrams (from "+p.size()+" f-diagrams)");
                p = isoFree.apply(newP, null);
                newP.clear();
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
                for (Graph g : p) {
                    boolean ap = hap.check(g);
                    boolean con = hap.isConnected();
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
                p = newP;
    
                // we don't need to re-isofree p, we know that's still good.
                // some of our new disconnected diagrams might condense with the old ones
                disconnectedP = isoFree.apply(maxIsomorph.apply(disconnectedP, mip), null);
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
                    if (g.nodeCount() > 3 && !graphHasEdgeColor(g, mBond) && !graphHasEdgeColor(g, MBond) && isBi.check(g) && 
                            (!doExchange || (!graphHasEdgeColor(g, excBond) && !graphHasEdgeColor(g, mxcBond)))) {
                        bcP.add(g);
                        continue;
                    }
                }
                Set<Graph> bcSubst = new HashSet<Graph>();
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
                        bcSubst.clear();
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
//                        p = isoFree.apply(compSubst.apply(p, csp), null);
                        disconnectedP = isoFree.apply(maxIsomorph.apply(compSubst.apply(disconnectedP, csp), mip), null);
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
    
                    newDisconnectedP[i+1] = isoFree.apply(maxIsomorph.apply(newDisconnectedP[i+1], mip), null);
    
                    disconnectedP.addAll(newDisconnectedP[i]);
                }
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
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        // group together all n-point diagrams (diagrams with n field nodes)
        // that evaluate to infinity (due to insufficient connectivity).
        // these are the diagrams that must be evaluated together during MSMC
        // at nth order.  all other n-point diagrams can be evaluated as
        // products of smaller diagrams
        Factor factor = new Factor();
        MulScalar mulScalar = new MulScalar();
        IsoFree isoFree = new IsoFree();

        Set<Graph> mP = new HashSet<Graph>();
        Set<Graph> pairP = new HashSet<Graph>();
        for (Graph g : p) {
            g = factor.apply(g,mfp);
            if (!graphHasEdgeColor(g, mBond) && !graphHasEdgeColor(g, mxcBond)) {
                if (doExchange) {
                    // we need to look for single points that are internally multibonded
                    boolean isMulti = false;
                    for (Node node : g.nodes()) {
                        if (node.getColor() >= 'M' && node.getColor() <= 'Z') {
                            isMulti = true;
                            break;
                        }
                    }
                    if (!isMulti) {
                        pairP.add(g);
                        continue;
                    }
                }
                else {
                    pairP.add(g);
                    continue;
                }
            }
            mP.add(g);
        }
        MulFlexibleParameters mfpc = MulFlexibleParameters.makeParameters(new char[]{nodeColor}, (byte)n);
        if (flex) {
            // we just want to force no-superimposing here.
            mfpc = mfp;
        }
        ComponentSubst compSubst = new ComponentSubst();
        // we need to start with i=1 so we catch fully exchanged diagrams
        HashSet<char[]> bicompSeen = new HashSet<char[]>();
        Set<Graph> mSubst = new HashSet<Graph>();
        // now we make a second pass where we actually make substitutions
        multiP = new HashSet<Graph>();
        multiP.addAll(mP);
        List<ComponentSubstParameters> cspUnList = new ArrayList<ComponentSubstParameters>();
        for (int i=1; i<=n; i++) {
            if (!doExchange && i<3) continue;
            // find all multi graphs of size i without a root point
            // then, Mi = mi1 + Mi2 + Mi3 + ...
            // where mi1, Mi2, Mi3... are the diagrams we find here.
            // mi1 is the fully connected diagram with mBonds, while Mi2, Mi3, etc.
            //  will be diagrams containing smaller multibody components composed of MBonds
            bicompSeen.clear();
            while (true) {
                char[] thisBicomp = null;
                Graph gm = null;
                Graph gM = null;
                Coefficient mcoef = null;
                mSubst.clear();
                Set<Graph> mUnSubst = new HashSet<Graph>();
                for (Graph g : mP) {
                    if (g.nodeCount() != i || NumRootNodes.value(g) > 0) {
                        continue;
                    }
                    if (doExchange) {
                        char[] iBicomp = new char[i];
                        boolean hasRoot = false;
                        for (byte j=0; j<i; j++) {
                            iBicomp[j] = g.getNode(j).getColor();
                            if (iBicomp[j] >= 'a' && iBicomp[j] <= 'z') {
                                // root point in disguise
                                hasRoot = true;
                                break;
                            }
                            if (iBicomp[j] >= 'M') {
                                iBicomp[j] -= 'M' - 'A';
                            }
                        }
                        if (hasRoot) continue;
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
                    if (g.edgeCount() < i*(i-1)/2) {
                        // this is a lower order (disconnected) diagram
                        mUnSubst.add(g.copy());
                        g.coefficient().multiply(new CoefficientImpl(-1));
                        mSubst.add(g);
                    }
                    else {
                        gm = g.copy();
                        // this is the large diagram in this set.
                        for (Edge e : g.edges()) {
                            if (e.getColor() == mBond) {
                                e.setColor(MBond);
                            }
                            else if (e.getColor() == mxcBond) {
                                e.setColor(MxcBond);
                            }
                            else {
                                throw new RuntimeException("oops");
                            }
                        }
                        gM = g;
                        mSubst.add(g);
                        mUnSubst.add(gm.copy());
                        mcoef = new CoefficientImpl(1);
                        mcoef.divide(g.coefficient());
                    }
                }
                if (doExchange && thisBicomp == null) break;
                MulScalarParameters msp = new MulScalarParameters(mcoef);
                mSubst = mulScalar.apply(mSubst, msp);
                mUnSubst = mulScalar.apply(mUnSubst, msp);
    
                // replace all mBond groups of size i with MBond
                // now mi1 = Mi - Mi2 - Mi3 - ...
                // where mi1 is the fully connected diagram of size i
                ComponentSubstParameters csp = new ComponentSubstParameters(gm, mSubst, mfpc);
                mP = isoFree.apply(compSubst.apply(mP, csp), null);

                // do a more targeted replacement for multiP
                Set<Graph> newMultiP = new HashSet<Graph>();
                for (Graph gmp : multiP) {
                    if (gmp.nodeCount() > i) {
                        newMultiP.addAll(compSubst.apply(gmp, csp));
                    }
                    else {
                        newMultiP.add(gmp);
                    }
                }
                multiP = isoFree.apply(newMultiP, null);
                // multiP will now contain M-bond and m-graph
                // we will remove any M-bond graphs that are finite
                cspUnList.add(new ComponentSubstParameters(gM, mUnSubst, mfpc));
                
                if (!doExchange) break;
            }
        }
        Set<Graph> newMultiP = new HashSet<Graph>();
        for (Graph gmp : multiP) {
            if (NumRootNodes.value(gmp) == 0) {
                newMultiP.add(gmp);
            }
        }

        // multiP is now the set of multi-body diagrams, but with the parts
        // that are products of Multi-bond graphs removed.  now substitute back
        // M => m
        multiP = newMultiP;
        for (int i=cspUnList.size()-1; i>=0; i--) {
            multiP = isoFree.apply(compSubst.apply(multiP, cspUnList.get(i)), null);
        }

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
        protected final boolean doExchange;
        public ArticulatedAt0(boolean doExchange) {
            this.doExchange = doExchange;
        }
        public boolean check(Graph graph) {
            if (doExchange) {
                char prevColor = graph.getNode((byte)0).getColor();
                if (prevColor >= 'a' && prevColor <= 'z') prevColor += 'A' - 'a';
                for (byte i=1; i<graph.nodeCount(); i++) {
                    char color = graph.getNode(i).getColor();
                    if (color >= 'a' && color <= 'z') color += 'A' - 'a';
                    if (color < prevColor) return false;
                    prevColor = color;
                }
            }
            if (isBi.check(graph)) return true;
            boolean ap = hap.check(graph);
//            boolean con = hap.isConnected();
            if (!ap) return true;
            boolean articulatedAtA = false;
            for (byte node : hap.getArticulationPoints()) {
                if (node == 0) return true;
                if (graph.getNode(node).getColor() == 'A') articulatedAtA = true;
            }
            return !articulatedAtA; 
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
