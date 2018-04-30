/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IdenticalGraphFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.BitmapFactory;
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
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.Decorate;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DecorateWertheim2SiteRho.DecorateWertheimParameters2Site;
import etomica.graph.operations.DecorateWertheim2SiteRho;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.FactorOnce;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.PCopy;
import etomica.graph.operations.ReeHoover;
import etomica.graph.operations.ReeHoover.ReeHooverParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.RelabelWertheim;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.SplitOne;
import etomica.graph.operations.SplitOne.SplitOneParameters;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterSumShell;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.BondConnectedSubstitution.BondConnectedSubstitutionParameters2Site;

/**
 * Wertheim diagram generator for 2 site model
 *
 * @author Hye Min Kim
 */
public class WertheimDiagrams2SiteRho {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean bondDecomp;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP,pWertheim2Site;
    protected Set<Graph> rho, rho0, rhoA,rhoB,rhoAB,rhoAB0;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Set<Graph>> cancelMap;
    protected boolean doShortcut;
    protected char fBond, eBond, mBond, capFBond, fRBond, capFABBond,capFBABond;

    public static void main(String[] args) {
    	long t1 = System.currentTimeMillis();
    	GraphImpl.useReverseEdges = true;
    	MetadataImpl.rootPointsSpecial = true;
        final int n = 4;
        boolean multibody = false;
        boolean bondDecomp = true;
        boolean flex = false;
        WertheimDiagrams2SiteRho virialDiagrams = new WertheimDiagrams2SiteRho(n, multibody, flex, bondDecomp, true);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.makeWertheimDiagrams();
        long t2 = System.currentTimeMillis();
        System.out.println("time taken "+(t2-t1)/1000.0);
    }
    
    public WertheimDiagrams2SiteRho(int n, boolean multibody, boolean flex, int numSite, boolean bondDecomp) {
        this(n, multibody, flex, bondDecomp, false);
    }
    
    public WertheimDiagrams2SiteRho(int n, boolean multibody, boolean flex, boolean bondDecomp, boolean interactive) {
        this.multibody = multibody;
        this.flex = flex;
        this.bondDecomp = bondDecomp;
        this.n = n;
        this.isInteractive = interactive;
        doShortcut = false;
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
    }
    
    public void setDoShortcut(boolean newDoShortcut) {
        doShortcut = newDoShortcut;
        if (multibody && doShortcut) {
            throw new RuntimeException("shortcut only works for pairwise models");
        }
    }

    public Set<Graph> getVirialGraphs() {
        if (p == null) {
            makeWertheimDiagrams();
        }
        return p;
    }

    public Set<Graph> getMSMCGraphs(boolean connectedOnly) {
        if (p == null) {
            makeWertheimDiagrams();
        }
        GraphList pn = makeGraphList();
        GraphList allP = makeGraphList();
        allP.addAll(p);
        if (!connectedOnly) {
            allP.addAll(cancelP);
        }
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
    
    public Map<Graph,Set<Graph>> getCancelMap() {
        if (p == null) {
            makeWertheimDiagrams();
        }
        return cancelMap;
    }

    public ClusterSum makeWertheimCluster(MayerFunction f, MayerFunction e) {//Mayer sampling
        if (p == null) {
            makeWertheimDiagrams();
        }
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        Set<Graph> pn = getMSMCGraphs(false);
        for (Graph g : pn) {
            populateEFBonds(g, allBonds, false);
            if (flex) {
                populateEFBonds(g, allBonds, true);
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
            return new ClusterSumEF(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{e});
        }
        else if (!multibody) {
            return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f});
        }
        return null;
    }
    
    public void populateEFBonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap) {//Mayer sampling
        ArrayList<int[]> fbonds = new ArrayList<int[]>();
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
        ArrayList<int[]> fRbonds = new ArrayList<int[]>();
        ArrayList<int[]> capFbonds = new ArrayList<int[]>();
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
                    else if (edgeColor == fRBond) {
                        fRbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFBond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFABBond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFBABond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else {
                        throw new RuntimeException("oops, unknown bond "+edgeColor);
                    }
                }
            }
        }
        if (!flex && n > 3) {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0]),ebonds.toArray(new int[0][0])}));
        }
        else {
            allBonds.add(new ClusterBonds(flex ? n+1 : n, new int[][][]{fbonds.toArray(new int[0][0])}));
        }

    }

    public ClusterSumShell[] makeSingleVirialClusters(ClusterSum coreCluster, MayerFunction e, MayerFunction f) {//Mayer sampling
        if (p == null) {
            makeWertheimDiagrams();
        }
        ArrayList<ClusterSumShell> allClusters = new ArrayList<ClusterSumShell>();
        double[] w = new double[]{1};
        if (flex) {
            w = new double[]{0.5,0.5};
        }
        Set<Graph> pn = getMSMCGraphs(true);
        for (Graph g : pn) {
            double[] thisW = w;
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            populateEFBonds(g, allBonds, false);
            if (flex) {
                populateEFBonds(g, allBonds, true);
            }
            if (flex && cancelMap.get(g) != null) {
                Set<Graph> cgSet = cancelMap.get(g);
                for (Graph cg : cgSet) {
                    populateEFBonds(cg, allBonds, false);
                    populateEFBonds(cg, allBonds, true);
                }
                thisW = new double[2+2*cgSet.size()];
                thisW[0] = thisW[1] = 0.5;
                for (int i=2; i<thisW.length; i++) {
                    thisW[i] = -0.5/cgSet.size();
                }
            }
            if (n > 3 && !flex && !multibody) {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{e,f}));
            }
            else if (!multibody) {
                allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{f}));
            }
        }
        return allClusters.toArray(new ClusterSumShell[0]);
    }
    
    public Set<Graph> getExtraDisconnectedVirialGraphs() {//Mayer sampling
        if (p == null) {
            makeWertheimDiagrams();
        }
        GraphList dpn = makeGraphList();
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

    public static GraphList makeGraphList() {
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList graphList = new GraphList(comp);
        return graphList;
    }        
    
    public void makeRhoDiagrams() {
        char nodeColor = COLOR_CODE_0;
        char[] flexColors = new char[0];
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

        GraphList topSet = makeGraphList();

        char oneBond = 'o';
        mBond = 'm';  // multi-body
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        fRBond = 'B';
        capFBond = 'U';
        capFABBond = 'F';
        capFBABond = 'G';
        ArrayList<Character> associationBond = new ArrayList<Character>();
        associationBond.add(capFABBond);
        associationBond.add(capFBABond);
        MetadataImpl.edgeColorPairs.add(associationBond);
        
        if (doShortcut) {
            // just take lnfXi to be the set of connected diagrams
            fBond = 'A';
            eBond = 'e';

            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(eBond, 2);
            colorOrderMap.put(fBond, 3);
            colorOrderMap.put(fRBond, 4);
            colorOrderMap.put(capFABBond, 5);
            colorOrderMap.put(capFBABond, 6);
            colorOrderMap.put(capFBond, 7);

            IsConnected isCon = new IsConnected();
            for (int i=1; i<n+1; i++) {
                GraphIterator iter = new PropertyFilter(new StoredIterator((byte)i), isCon);
                msp = new MulScalarParameters(1, (int)SpecialFunctions.factorial(i));
                while (iter.hasNext()) {
                    lnfXi.add(mulScalar.apply(iter.next(), msp));
                }
            }
        }
        else {
            
            eBond = COLOR_CODE_0;//color of edge
            fBond = 'f';

            Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(eBond, 2);
            colorOrderMap.put(fBond, 3);
            colorOrderMap.put(fRBond, 4);
            colorOrderMap.put(capFABBond, 5);
            colorOrderMap.put(capFBABond, 6);
            colorOrderMap.put(capFBond, 7);
            
            for (byte i=1; i<n+1; i++) {
                Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                eXi.add(g);
    
                if (multibody && i>2) {
                    g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                    g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                    for (Edge e : g.edges()) {
                        e.setColor(mBond);
                    }
                    eXi.add(g);
                }
            }
    
            if (isInteractive) {
                System.out.println("Xi with e bonds");
                topSet.addAll(eXi);
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                System.out.println("eXi has "+eXi.size()+" graphs");
                //ClusterViewer.createView("Xi with e bonds", topSet);
            }
    
            Split split = new Split();
            SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
            Set<Graph> setOfSubstituted = split.apply(eXi, bonds);         

    
            DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
            DeleteEdge deleteEdge = new DeleteEdge();
            //set of full star diagrams with f bonds
            Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
            setOfSubstituted.clear();

            if (isInteractive) {
                System.out.println("\nXi with f bonds");
                topSet.clear();
                topSet.addAll(fXi);
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                System.out.println("fXi has "+fXi.size()+" graphs");
               //ClusterViewer.createView("fXi", topSet);
            }
            
            Set<Graph> fXipow = new HashSet<Graph>();
            fXipow.addAll(fXi);
            
            for (int i=1; i<n+1; i++) {
                lnfXi.addAll(fXipow);
                lnfXi = isoFree.apply(lnfXi, null);
                msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
                fXipow = mulScalar.apply(mulFlex.apply(fXipow, fXi, mfp), msp);
            }
            fXipow = null;
            fXi.clear();
        }
        if (isInteractive) {
            topSet.clear();
            topSet.addAll(lnfXi);
            System.out.println("\nlnfXi");
            for (Graph g : topSet) {
                System.out.println(g);
            }
            System.out.println("lnfXi has "+lnfXi.size()+" graphs");
            //ClusterViewer.createView("lnfXi", topSet);
        }
    	for (Graph g:lnfXi){
        	IsConnected isConnected = new IsConnected();
        	if (!isConnected.check(g)){
        		System.out.println(" g fXipow "+g);
        		throw new RuntimeException();
        	}
    	}
    	//f = fR+FAB+FBA
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(fBond, fRBond, capFBond);
        lnfXi = split.apply(lnfXi, bonds);
        bonds = new SplitParameters(capFBond, capFABBond, capFBABond);
        lnfXi = split.apply(lnfXi, bonds);   

        DifByNode opzdlnXidz = new DifByNode();//z*dlnXi/dz
        DifParameters difParams = new DifParameters(nodeColor);
        rho = opzdlnXidz.apply(lnfXi, difParams);
        
        rho = isoFree.apply(rho, null);
        rho0 = new HashSet<Graph>();
        rhoA = new HashSet<Graph>();
        rhoB = new HashSet<Graph>();
        rhoAB = new HashSet<Graph>();
        rhoAB0 = new HashSet<Graph>();
        
        for (Graph g : rho) {
        	boolean hascapFABBond = false;
        	boolean hascapFBABond = false;
        	byte rootNode = 0;
        	for (byte i = 0; i<g.nodeCount(); i++){
        		if (g.getNode(i).getType()== Metadata.TYPE_NODE_ROOT){
        			rootNode = i;
        			break;
        		}
        	}

        	for (int i = 0;i<g.nodeCount(); i++){
        		if (i == rootNode || !g.hasEdge(rootNode, (byte)i))continue;
        		if (g.getEdge(rootNode, (byte)i).getColor() == capFABBond)hascapFABBond = true;
        		if (g.getEdge(rootNode, (byte)i).getColor() == capFBABond)hascapFBABond = true;
        	}
        	if (hascapFABBond&&hascapFBABond)
        		rhoAB.add(g);
        	else if (hascapFABBond)
        		rhoA.add(g);
        	else if (hascapFBABond)
        		rhoB.add(g);
        	else 
        		rho0.add(g);
        }
        
        Set<Graph> newRhoA = new HashSet<Graph>();
        for (Graph g : rhoA){
        	byte rootId = 0;
        	for (Node node:g.nodes()){
        		if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
        	}
        	if (rootId == 0){
        		newRhoA.add(g);
        		continue;
        	}
            RelabelWertheim relabel = new RelabelWertheim();
            byte[] permutations = new byte[g.nodeCount()];
            for (int i=0; i<permutations.length; i++) {
                permutations[i] = (byte)i;
            }
            permutations[0] = rootId;
            permutations[rootId] = 0;
            RelabelParameters rp = new RelabelParameters(permutations);
            g = relabel.apply(g, rp);
            newRhoA.add(g);
        }
        rhoA.clear();
        rhoA.addAll(newRhoA);
        
        Set<Graph> newRhoB = new HashSet<Graph>();
        for (Graph g : rhoB){
        	byte rootId = 0;
        	for (Node node:g.nodes()){
        		if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
        	}
        	if (rootId == 0){
        		newRhoB.add(g);
        		continue;
        	}
            RelabelWertheim relabel = new RelabelWertheim();
            byte[] permutations = new byte[g.nodeCount()];
            for (int i=0; i<permutations.length; i++) {
                permutations[i] = (byte)i;
            }
            permutations[0] = rootId;
            permutations[rootId] = 0;
            RelabelParameters rp = new RelabelParameters(permutations);
            g = relabel.apply(g, rp);
            newRhoB.add(g);
        }
        rhoB.clear();
        rhoB.addAll(newRhoB);
        
        Set<Graph> newRhoAB = new HashSet<Graph>();
        for (Graph g : rhoAB){
        	byte rootId = 0;
        	for (Node node:g.nodes()){
        		if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
        	}
        	if (rootId == 0){
        		newRhoAB.add(g);
        		continue;
        	}
            RelabelWertheim relabel = new RelabelWertheim();
            byte[] permutations = new byte[g.nodeCount()];
            for (int i=0; i<permutations.length; i++) {
                permutations[i] = (byte)i;
            }
            permutations[0] = rootId;
            permutations[rootId] = 0;
            RelabelParameters rp = new RelabelParameters(permutations);
            g = relabel.apply(g, rp);
            newRhoAB.add(g);
        }
        rhoAB.clear();
        rhoAB.addAll(newRhoAB);
        
        Set<Graph> rho0pow = new HashSet<Graph>();
        Set<Graph> rho0m1 = new HashSet<Graph>();

        for (Graph g : rho0){
        	if (g.nodeCount() > 1)rho0m1.add(g);
        }

        rho0pow.addAll(rho0m1);

        Set<Graph> rho0m1pow = new HashSet<Graph>();

        msp = new MulScalarParameters(new CoefficientImpl(-1,1));
        rho0m1 = mulScalar.apply(rho0m1, msp);
        rho0m1pow.addAll(rho0m1);
    }
    
    public void makeWertheimDiagrams() {
        
        char nodeColor = COLOR_CODE_0;
        char[] flexColors = new char[0];
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

        GraphList topSet = makeGraphList();

        char oneBond = 'o';
        mBond = 'm';  // multi-body
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();

        if (doShortcut && !flex) {
            // just take lnfXi to be the set of connected diagrams
            fBond = COLOR_CODE_0;
            eBond = 'e';

            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(eBond, 2);
            colorOrderMap.put(fBond, 3);

            // skip directly to p diagrams
            p = new HashSet<Graph>();
            IsBiconnected isBi = new IsBiconnected();
            for (int i=1; i<n+1; i++) {
                GraphIterator iter = new PropertyFilter(new StoredIterator((byte)i), isBi);
                msp = new MulScalarParameters(1-i, (int)SpecialFunctions.factorial(i));
                while (iter.hasNext()) {
                    p.add(mulScalar.apply(iter.next(), msp));
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
    
            MulFlexibleParameters mfpnm1 = MulFlexibleParameters.makeParameters(flexColors, (byte)(n-1));//n-1 field point
            MulFlexibleParameters mfpnm1zWertheim = MulFlexibleParameters.makeParametersOnlyRootPt(flexColors, (byte)(n-1));
            z.addAll(allRho[1]);
            for (int i=2; i<n+1; i++) {
                Set<Graph>[] zPow = new Set[n+1];
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

            HashSet<Graph>[] allRho0 = new HashSet[n+1];
            for (int i=0; i<n+1; i++) {
                allRho0[i] = new HashSet<Graph>();
            }
            for (Graph g : rho0) {
                allRho0[g.nodeCount()].add(g);
            }
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoAB);
                System.out.println("\nrhoAB(z)");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoAB(z)", topSet);
                System.out.println("rhoAB(z) has "+rhoAB.size()+" graphs");
            }
            Set<Graph> zWertheim = new HashSet<Graph>();

            zWertheim.addAll(allRho0[1]);//add all of 1 point diagram in rho0
            msp = new MulScalarParameters(new CoefficientImpl(-1,1));
            for (int i=2; i<n+1; i++) {
                Set<Graph>[] zPow = new Set[n+1];
                zPow[1] = new HashSet<Graph>();
                zPow[1].addAll(zWertheim);
                for (int j=2; j<i+1; j++) {
                    zPow[j] = new HashSet<Graph>();
                    zPow[j] = isoFree.apply(mulFlex.apply(zPow[j-1], zWertheim, mfpnm1zWertheim), null);
                }
                zWertheim = new HashSet<Graph>();
                zWertheim.addAll(allRho0[1]);
                for (int j=2; j<i+1; j++) {
                    zWertheim.addAll(mulScalar.apply(mulFlex.apply(allRho0[j], zPow[j], mfpnm1zWertheim), msp));
                }
            }
    
            zWertheim = isoFree.apply(zWertheim, null);
            
            Decorate decorate = new Decorate();
            DecorateParameters dp = new DecorateParameters(nodeColor, mfp);
            DecorateParameters dpnm1 = new DecorateParameters(0,mfpnm1zWertheim);

        	p = isoFree.apply(decorate.apply(lnfXi, zWertheim, dp),null);
            lnfXi.clear();

	       for(Graph g:rho){
	    	   g.setNumFactors(1);
	    	   g.addFactors(new int[]{g.nodeCount()});
	       } 
           for(Graph g:rhoA){
        	   g.setNumFactors(1);
        	   g.addFactors(new int[]{g.nodeCount()});
           }
           for(Graph g:rhoB){
        	   g.setNumFactors(1);
        	   g.addFactors(new int[]{g.nodeCount()});
           }
           for(Graph g:rhoAB){
        	   g.setNumFactors(1);
        	   g.addFactors(new int[]{g.nodeCount()});
           }
           for(Graph g:zWertheim){
        	   g.setNumFactors(1);
           }

            rho = isoFree.apply(decorate.apply(rho, zWertheim, dpnm1),null);
            rhoA = isoFree.apply(decorate.apply(rhoA, zWertheim, dpnm1),null);
            rhoB = isoFree.apply(decorate.apply(rhoB, zWertheim, dpnm1),null);
            rhoAB = isoFree.apply(decorate.apply(rhoAB, zWertheim, dpnm1),null);
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rho0);
                System.out.println("\nrho0(z)");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                System.out.println("rho0(z) has "+rho0.size()+" graphs");
                //ClusterViewer.createView("rho0(z)", topSet);
            }
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(zWertheim);
                System.out.println("\nzWertheim");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("zWertheim", topSet);
                System.out.println("zWertheim has "+zWertheim.size()+" graphs");
            }
            
            zWertheim.clear();
            
	        for(Graph g:rho){
	        	g.setNumFactors(2);
	        	g.addFactors(new int[]{g.nodeCount(),0});
	         }          
	         for(Graph g:rhoA){
	        	 g.setNumFactors(2);
	        	 g.addFactors(new int[]{g.nodeCount(),0});
	         }
	         for(Graph g:rhoB){
	        	 g.setNumFactors(2);
	        	 g.addFactors(new int[]{g.nodeCount(),0});
	         }
	         for(Graph g:rhoAB){
	        	 g.setNumFactors(2);
	        	 g.addFactors(new int[]{g.nodeCount(),0});
	         }

            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Factor factor = new Factor();
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoA);
                System.out.println("\nrhoA(rho0)");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoA(rho0)", topSet);
                System.out.println("rhoA(rho0) has "+rhoA.size()+" graphs");
            }
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoAB);
                System.out.println("\nrhoAB(rho0)");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoAB(rho0)", topSet);
                System.out.println("rhoAB(rho0) has "+rhoAB.size()+" graphs");
            }

            DecorateWertheim2SiteRho decorateWertheim2Site = new DecorateWertheim2SiteRho();
            DecorateWertheimParameters2Site dpWertheim2Site = new DecorateWertheimParameters2Site(mfp, capFABBond, capFBABond, rhoA, rhoB, rhoAB);//decorate rho point with rho0C1
            DecorateWertheimParameters2Site drhoWertheim2Site = new DecorateWertheimParameters2Site(mfpnm1, capFABBond, capFBABond, rhoA, rhoB, rhoAB);
            pWertheim2Site = decorateWertheim2Site.apply(p, dpWertheim2Site);
            rhoA = decorateWertheim2Site.apply(rhoA,drhoWertheim2Site);
            rhoB = decorateWertheim2Site.apply(rhoB,drhoWertheim2Site);
            rhoAB = decorateWertheim2Site.apply(rhoAB,drhoWertheim2Site);
            rho = decorateWertheim2Site.apply(rho,drhoWertheim2Site);
            
            
            pWertheim2Site = isoFree.apply(pWertheim2Site, null);
            rhoA = isoFree.apply(rhoA,null);
            rhoB = isoFree.apply(rhoB,null);
            rhoAB = isoFree.apply(rhoAB,null);
            rho = isoFree.apply(rho,null);
            
            BondConnectedSubstitution bondConnected = new BondConnectedSubstitution();
            BondConnectedSubstitutionParameters2Site bcp2Site = new BondConnectedSubstitutionParameters2Site(fRBond, eBond);
            pWertheim2Site = bondConnected.apply(pWertheim2Site,bcp2Site);
            rhoA = bondConnected.apply(rhoA,bcp2Site);
            rhoB = bondConnected.apply(rhoB,bcp2Site);
            rhoAB = bondConnected.apply(rhoAB,bcp2Site);
            rho = bondConnected.apply(rho,bcp2Site);

            Set<Graph> newpWertheim2Site = new HashSet<Graph>();
            for (Graph g : pWertheim2Site) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newpWertheim2Site.add(gf);
                }
                else {
                	newpWertheim2Site.add(g.copy());
                }
            }
            pWertheim2Site = newpWertheim2Site;
            MetadataImpl.rootPointsSpecial = false;
            pWertheim2Site = isoFree.apply(pWertheim2Site, null);
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(pWertheim2Site);
                System.out.println("\npWertheim2Site");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                ClusterViewer.createView("pWertheim2Site", topSet);
                System.out.println("pwertheim has "+pWertheim2Site.size()+" graphs");
            }
            
            Set<Graph> newRhoA = new HashSet<Graph>();
            for (Graph g : rhoA) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newRhoA.add(gf);
                }
                else {
                	newRhoA.add(g.copy());
                }
            }
            rhoA = newRhoA;
            MetadataImpl.rootPointsSpecial = false;
            rhoA = isoFree.apply(rhoA, null);
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoA);
                System.out.println("\nrhoA");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoA", topSet);
                System.out.println("rhoA has "+rhoA.size()+" graphs");
            }
            
            Set<Graph> newRhoB = new HashSet<Graph>();
            for (Graph g : rhoB) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newRhoB.add(gf);
                }
                else {
                	newRhoB.add(g.copy());
                }
            }
            rhoB = newRhoB;
            MetadataImpl.rootPointsSpecial = false;
            rhoB = isoFree.apply(rhoB, null);
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoB);
                System.out.println("\nrhoB");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoB", topSet);
                System.out.println("rhoB has "+rhoB.size()+" graphs");
            }
            
            Set<Graph> newRhoAB = new HashSet<Graph>();
            for (Graph g : rhoAB) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newRhoAB.add(gf);
                }
                else {
                	newRhoAB.add(g.copy());
                }
            }
            rhoAB = newRhoAB;
            MetadataImpl.rootPointsSpecial = false;
            rhoAB = isoFree.apply(rhoAB, null);
            
            rhoAB0.addAll(rhoAB);
            MulFlexibleParameters mfpCi = MulFlexibleParameters.makeParametersOnlyRootPt(new char[]{nodeColor}, (byte)(n-1));//truncated in nth order, avoid superimposing every point
            Set<Graph> rhoArhoB = mulFlex.apply(rhoA, rhoB, mfpCi);   
            int[] rhoArhoBFactor = new int[]{-1,0,0,0};
            for (Graph g:rhoArhoB){
            	g.addFactors(rhoArhoBFactor);
            }
            rhoAB0.addAll((mulScalar.apply(rhoArhoB, msp)));
            rhoAB0 = isoFree.apply(rhoAB0, null);
            
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoAB);
                System.out.println("\nrhoAB");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoAB", topSet);
                System.out.println("rhoAB has "+rhoAB.size()+" graphs");
            }
            
            if (isInteractive) {
                topSet.clear();
                topSet.addAll(rhoAB0);
                System.out.println("\nrhoAB0");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("rhoAB0", topSet);
                System.out.println("rhoAB0 has "+rhoAB0.size()+" graphs");
            }
            
            Set<Graph> newRho = new HashSet<Graph>();
            for (Graph g : rho) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);
                    newRho.add(gf);
                }
                else {
                	newRho.add(g.copy());
                }
            }
            rho = newRho;
            MetadataImpl.rootPointsSpecial = false;
            rho = isoFree.apply(rho, null);

            // clear these out -- we don't need them and (in extreme cases) we might need the memory

                
            lnfXi.clear();
            rho.clear();
            z.clear();
        }

        Set<Graph> newP = new HashSet<Graph>();

        // attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Set<Graph>>();
        cancelP = new GraphList();
        disconnectedP = new HashSet<Graph>();
        if (!flex) {
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Factor factor = new Factor();
            newP.clear();
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);//
                    newP.add(gf);
                }
                else {
                    newP.add(g.copy());
                }
            }
            //p = isoFree.apply(newP, null);

            // perform Ree-Hoover substitution (brute-force)
            if (doReeHoover) {
                if (doShortcut) {
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
        }
        else {
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Relabel relabel = new Relabel();
            FactorOnce factor = new FactorOnce();
            // we can take all permutations (all ways to factor the diagram) here
            // it can make the MSMC a bit more efficient, but (in practice) the difference
            // seems to be negligible
            boolean allPermutations = false;
            FactorOnceParameters fop = new FactorOnceParameters((byte)0, allPermutations);
            newP.clear();
            msp = new MulScalarParameters(-1, 1);
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if (con && ap) {
                    if (!hap.getArticulationPoints().contains(0)) {
                        byte[] permutations = new byte[g.nodeCount()];
                        for (int i=0; i<permutations.length; i++) {
                            permutations[i] = (byte)i;
                        }
                        permutations[0] = hap.getArticulationPoints().get(0);
                        permutations[hap.getArticulationPoints().get(0)] = 0;
                        RelabelParameters rp = new RelabelParameters(permutations);
                        g = relabel.apply(g, rp);
                    }
                    // newP will contain connected diagrams
                    g = g.copy();
                    newP.add(g);
                    Set<Graph> gf = factor.apply(g, fop);
                    disconnectedP.addAll(gf);
                    gf = mulScalar.apply(gf, msp);
                    cancelP.addAll(gf);
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

            // we want to condense cancelP (in case multiple diagrams were factored into the
            // same one -- is that even possible?), but want to be careful not to permute bonds.
            PCopy pcopy = new PCopy();
            IteratorWrapper wrapper = new IteratorWrapper(pcopy.apply(cancelP, null).iterator());
            GraphIterator isomorphs = new IdenticalGraphFilter(wrapper);
            cancelP = new GraphList();
            while (isomorphs.hasNext()) {
                cancelP.add(isomorphs.next());
            }
        }


        GraphList pFinal = makeGraphList();
        pFinal.addAll(p);
        p = pFinal;
        GraphList disconnectedPFinal = makeGraphList();
        disconnectedPFinal.addAll(disconnectedP);
        disconnectedP = disconnectedPFinal;

    }
}
