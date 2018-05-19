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

import etomica.graph.iterators.StoredIterator;
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
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.Unfactor;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
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

/**
 * Mayer's theory on ionic solution
 * f-bond is splitted to k-bond, (kg+g),(kg2+g2)
 * diagrams with g chain(s) are reduced to lower order of diagrams(here they are just simply discarded)
 * gi is substituted by (qi - kgi) 
 * qi represent modified q bonds, gi represent REAL q bonds
 * even number of q bonds , the final result contains q1,q2,kg1,kg2bonds
 * 
 * 
 * @author shu
 * Date: April 20 2012
 */
public class MayerIonicDiagram2qbonds {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean bondDecomp;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Set<Graph> s;
    protected Set<Graph> rho;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Set<Graph>> cancelMap;
    protected boolean doShortcut;
    protected char fBond, eBond, mBond, kBond, GBond, G1Bond, G2Bond,g1Bond, g2Bond,kg1Bond, kg2Bond;
    protected char q1Bond,q2Bond;//modified q bond, actually
    protected char kg1BondMinus, kg2BondMinus;
    // Let's begin!
    public static void main(String[] args) {
    	GraphImpl.useReverseEdges = true;
    	MetadataImpl.rootPointsSpecial = true;
        final int n = 3;
        boolean multibody = false;
        boolean bondDecomp = true;
        boolean flex = false;
        MayerIonicDiagram2qbonds virialDiagrams = new MayerIonicDiagram2qbonds(n, multibody, flex, bondDecomp, true);
        virialDiagrams.setDoReeHoover(false);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.makeVirialDiagrams();
    }
    // constructor
    public MayerIonicDiagram2qbonds(int n, boolean multibody, boolean flex, int numSite, boolean bondDecomp) {
        this(n, multibody, flex, bondDecomp, false);
    }
    // the other constructor
    public MayerIonicDiagram2qbonds(int n, boolean multibody, boolean flex, boolean bondDecomp, boolean interactive) {
        this.multibody = multibody;
        this.flex = flex;
        this.bondDecomp = bondDecomp;
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
        doReeHoover = newDoReeHoover;
    }
    
    public void setDoShortcut(boolean newDoShortcut) {
        doShortcut = newDoShortcut;
        if (multibody && doShortcut) {
            throw new RuntimeException("shortcut only works for pairwise models");
        }
    }

    public Set<Graph> getVirialGraphs() {
        if (p == null) {
            makeVirialDiagrams();
        }
        return p;
    }

    public Set<Graph> getMSMCGraphs(boolean connectedOnly) {
        if (p == null) {
            makeVirialDiagrams();
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
            makeVirialDiagrams();
        }
        return cancelMap;
    }

    public ClusterSum makeVirialCluster(MayerFunction f, MayerFunction e) {//Mayer sampling
        if (p == null) {
            makeVirialDiagrams();
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
        ArrayList<int[]> kBonds = new ArrayList<int[]>();
        ArrayList<int[]> gBonds = new ArrayList<int[]>();
      
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
                    else if (edgeColor == kBond) {
                    	kBonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == GBond) {
                    	gBonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == g1Bond) {
                    	gBonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == g2Bond) {
                    	gBonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == kg1Bond) {
                    	gBonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == kg2Bond) {
                    	gBonds.add(new int[]{n1,n2});
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
            makeVirialDiagrams();
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
            makeVirialDiagrams();
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
        comp.addComparator(new ComparatorNumFieldNodes());// sort according to number of nodes
        comp.addComparator(new ComparatorBiConnected());// bi-connected diagrams first
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
        //these will show up in console after simulation
        GBond = 'X';
        G1Bond = 'x';
        G2Bond = 'Y';
  
        kBond = 'k';
        g1Bond = 'A';
        g2Bond = 'B';
              
        kg1Bond = 'D';
        kg2Bond = 'E';
        
        q1Bond = 'H';
        q2Bond = 'I';
        
        kg1BondMinus = 'L';
        kg2BondMinus = 'M';

        if (doShortcut) {
            // just take lnfXi to be the set of connected diagrams
            fBond = 'Z';
            eBond = 'e';
            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(eBond, 2);
            colorOrderMap.put(fBond, 3);
            colorOrderMap.put(kBond, 4);
            colorOrderMap.put(GBond, 5);
            colorOrderMap.put(G1Bond, 6);
            colorOrderMap.put(G2Bond, 7);
         
            colorOrderMap.put(g1Bond, 10);
            colorOrderMap.put(g2Bond, 11);
            colorOrderMap.put(kg1Bond, 13);
            colorOrderMap.put(kg2Bond, 14);
            colorOrderMap.put(q1Bond, 16);
            colorOrderMap.put(q2Bond, 17);
            colorOrderMap.put(kg1BondMinus, 19);
            colorOrderMap.put(kg2BondMinus, 20);

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
            colorOrderMap.put(kBond, 4);
            colorOrderMap.put(GBond, 5);
            colorOrderMap.put(G1Bond, 6);
            colorOrderMap.put(G2Bond, 7);
            
            colorOrderMap.put(g1Bond, 10);
            colorOrderMap.put(g2Bond, 11);
            colorOrderMap.put(kg1Bond,13);
            colorOrderMap.put(kg2Bond,14);
            colorOrderMap.put(q1Bond, 16);
            colorOrderMap.put(q2Bond, 17);
            colorOrderMap.put(kg1BondMinus, 19);
            colorOrderMap.put(kg2BondMinus, 20);
        
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
                ClusterViewer.createView("Xi with e bonds", topSet);
            }
    
            Split split = new Split();// e = 1 + f
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
        Metadata.COLOR_MAP.put(kBond, "black");
        Metadata.COLOR_MAP.put(g1Bond,"red");
        Metadata.COLOR_MAP.put(q1Bond,"red");
        Metadata.COLOR_MAP.put(kg1Bond,"red");
        Metadata.DASH_MAP.put(kg1Bond, 3);
     

        Metadata.COLOR_MAP.put(g2Bond, "green");
        Metadata.COLOR_MAP.put(q2Bond, "green");
        Metadata.COLOR_MAP.put(kg2Bond,"green");
        Metadata.DASH_MAP.put(kg2Bond, 3);
             
        Metadata.COLOR_MAP.put(kg1BondMinus,"yellow");
        Metadata.COLOR_MAP.put(kg2BondMinus,"yellow");
        DifByNode opzdlnXidz = new DifByNode();//rho = z*dlnXi/dz, functional derivative
        DifParameters difParams = new DifParameters(nodeColor);
        rho = isoFree.apply(opzdlnXidz.apply(lnfXi, difParams), null);
    }
    
    public void makeVirialDiagrams() {
        
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
            
            rho.clear();
            Set<Graph> z = new HashSet<Graph>();
            // r is rho , density
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
            // get p 
        	p = decorate.apply(lnfXi, z, dp);
            p = isoFree.apply(p, null);
            
            lnfXi.clear();
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Factor factor = new Factor();
            Set<Graph> newP = new HashSet<Graph>();
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
            System.out.println("\np");
            for (Graph g : p) {
                System.out.println(g);
            }
            ClusterViewer.createView("p", p);
            // get s 
            s = new HashSet<Graph>();
            for (Graph g : p) {
            	g = g.copy();
            	g.coefficient().multiply(new CoefficientImpl(1,1-g.nodeCount()));
           	  	s.add(g);
            }
            System.out.println("\ns from p");
            for (Graph g :s) {
                System.out.println(g);
            }
            ClusterViewer.createView("s from p", s);
            // split f bonds in s 
            Split split = new Split();
            SplitParameters bonds = new SplitParameters(fBond, kBond, GBond);
            s = split.apply(s, bonds);
            bonds = new SplitParameters(GBond, G1Bond, G2Bond);
            s = split.apply(s, bonds);  
            bonds = new SplitParameters(G1Bond, g1Bond, kg1Bond);
            s = split.apply(s, bonds); 
            bonds = new SplitParameters(G2Bond, g2Bond, kg2Bond);
            s = split.apply(s, bonds);  

            System.out.println("\ns after first splitting");
            for (Graph g :s) {
                System.out.println(g);
            }
            ClusterViewer.createView("s after first splitting", s);
            
            //find g bond chains and delete diagrams containing them
            Set<Graph> sTemp = new HashSet<Graph>();
            for (Graph g : s) {//loop over all graphs in s
            	boolean hasgChain = false; 
            	// find 2pt cycle diagram first
            	if (g.nodeCount() == 2) {// 2 pt cycle diagram:actually it is just g2, throw this away!
            		for (Edge e:g.edges()){
            			if ( e.getColor() == g2Bond ){
            				hasgChain = true;// label node that has two g bonds
            			}
            		}
            	}
            	
            	for (byte j = 0;j<(g.nodeCount());j++){// loop over all nodes in a given graph
            		Node node = g.getNode(j);
            		if ( g.getOutDegree(node.getId())==2 ) {// find out all nodes that have two edges 
            			// getOutNode :get the bonds 
            			byte index0 = g.getOutNode(node.getId(), (byte)0);
            			byte index1 = g.getOutNode(node.getId(), (byte)1);
            			char color0 = g.getEdge(node.getId(), index0).getColor();
            			char color1 = g.getEdge(node.getId(), index1).getColor();
            			if((color0 == g1Bond)&& (color1 == g1Bond)){ 
            				hasgChain = true;// label node that has two g bonds
            			} 
            		}
            	}
            	// if there is no g chain, then add such diagram in sTemp, discard otherwise
            	if(!hasgChain){
            		sTemp.add(g);
            	}
            }
            s = isoFree.apply(sTemp, null);//add diagrams in sTemp to s 
            sTemp.clear();// so that sTemp is clear and can be used later again
            
            System.out.println("\ns after deleting all diagrams that contain g bond chains");
            for (Graph g :s) {
                System.out.println(g);
            }
            ClusterViewer.createView("s after deleting all diagrams that contain g bond chains ", s);
            
            // substitute g with q.(NOW g represents exp(-kappa*r)/4*PI*r instead of exp(-alpha*r)/4*PI*r. g + kg = q ==> g = q - kg
            bonds = new SplitParameters(g1Bond, q1Bond, kg1BondMinus);
            s = split.apply(s, bonds);
            bonds = new SplitParameters(g2Bond, q2Bond, kg2BondMinus);
            s = split.apply(s, bonds); 
         
            // (-1) * diagrams with kg(i)BondMinus bond
            for (Graph g : s) {
    			g = g.copy();
    			for (Edge e:g.edges()){// loop over all edges to find kg1BondMinus OR kg2BondMinus OR kg3BondMinus OR kg4BondMinus
    				char color = e.getColor();
    				if((color == kg1BondMinus) || (color == kg2BondMinus) ){ 
    					g.coefficient().multiply(new CoefficientImpl(-1));    		
    				}
    			}
    			sTemp.add(g);
            }
            s = isoFree.apply(sTemp, null);
            sTemp.clear();
          
            // kg1Minus--> kg1, kg2Minus--> kg2, kg3Minus--> kg3, kg4Minus--> kg4
            for (Graph g : s) {
    			g = g.copy();
    			for (Edge e:g.edges()){  // loop over all edges in every diagram to find these three types of bonds
    				char color = e.getColor();
    				if (color == kg1BondMinus){
    					e.setColor(kg1Bond);  // kg1Minus--> kg1
    				} 
    				if (color == kg2BondMinus){ 
    					e.setColor(kg2Bond); // kg2Minus--> kg2
    				}
    			}
    			sTemp.add(g);
            }
            s = isoFree.apply(sTemp, null);
            sTemp.clear();
            System.out.println("\ns after second splitting, with q bonds and kg bonds");
            for (Graph g :s) {
                System.out.println(g);
            }
            ClusterViewer.createView("s after second splitting, with q bonds and kg bonds", s);
            
            //delete diagrams with odd number of qbonds, loop over all edges in every diagram to find
            for (Graph g : s) {
            	int qBond = 0;
            	g = g.copy();
            	for (Edge e:g.edges()){
    				if ( (e.getColor() == kg1Bond ) || (e.getColor() == q1Bond ) ){
    					qBond +=1;
    				} else if ( (e.getColor() == kg2Bond ) || (e.getColor() == q2Bond ) ){
    					qBond +=2;
    				} else if ( e.getColor() == kBond  ){
    					qBond += 0;
    				} else {
    		        	throw new RuntimeException("This is strange");
    		        }
            	}
            	if ( qBond % 2 == 0 ) {// qBond number is even
            		sTemp.add(g);
            	}
            }
            
        	s = isoFree.apply(sTemp, null);
            sTemp.clear();
        	System.out.println("\ns with only even number of q bonds, unordered");
            for (Graph g :s) {
            	System.out.println(g);
            }
            ClusterViewer.createView("s with only even number of q bonds, unordered", s);
            // topset is an ordered set of graphs
            if (isInteractive) {
            	System.out.println("\ns with only even number of q bonds, ordered");
                topSet.clear();
                topSet.addAll(s);
                for (Graph g : topSet) {
                	System.out.println(g);
                }
                ClusterViewer.createView("s with only even number of bonds, ordered", topSet);
            }
            
        }
    }
}
 