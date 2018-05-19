/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.io.FileWriter;
import java.io.IOException;
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
import etomica.graph.operations.CombineABSite.CombineABSiteParameters;
import etomica.graph.operations.Decorate;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DecorateWertheim3SiteRho;
import etomica.graph.operations.DecorateWertheim3SiteRho.DecorateWertheimParameters3Site;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorphWertheim;
import etomica.graph.operations.MaxIsomorphWertheim.MaxIsomorphWertheimParameters;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.RelabelWertheim;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.math.SpecialFunctions;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterSumShell;
import etomica.virial.MayerFunction;
import etomica.virial.cluster.BondConnectedSubstitution3Site.BondConnectedSubstitutionParameters3Site;

/**
 * Wertheim diagram generator for GCPM water (3 site model)
 *
 * @author Hye Min Kim
 */
public class WertheimDiagrams3SiteRho {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected final boolean bondDecomp;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Set<Graph> rho, rho0, rho1, rhoA,rhoB,rhoC,rhoAB,rhoAC,rhoBC,rhoABC;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Set<Graph>> cancelMap;
    protected boolean doShortcut;
    protected char fBond, eBond, mBond, mERBond, capFBond,mCapFBond, capF1Bond, mCapF1Bond, capF2Bond, mCapF2Bond, fRBond, capFACBond,mCapFACBond, capFBCBond,mCapFBCBond, capFCABond,mCapFCABond, capFCBBond,mCapFCBBond;

    public static void main(String[] args) {
        long t1 = System.currentTimeMillis();
    	GraphImpl.useReverseEdges = true;
    	MetadataImpl.rootPointsSpecial = true;
        final int n = 4;
        boolean multibody = true;
        boolean bondDecomp = true;
        boolean flex = false;
        WertheimDiagrams3SiteRho virialDiagrams = new WertheimDiagrams3SiteRho(n, multibody, flex, bondDecomp, true);
        virialDiagrams.setDoReeHoover(false);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.makeVirialDiagrams();
        long t2 = System.currentTimeMillis();
        System.out.println("time taken "+(t2-t1)/1000.0);
    }
   
    public WertheimDiagrams3SiteRho(int n, boolean multibody, boolean flex, boolean bondDecomp) {
        this(n, multibody, flex, bondDecomp, false);
    }
   
    public WertheimDiagrams3SiteRho(int n, boolean multibody, boolean flex, boolean bondDecomp, boolean interactive) {
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
        ArrayList<int[]> fRbonds = new ArrayList<int[]>();
        ArrayList<int[]> capFbonds = new ArrayList<int[]>();
        ArrayList<int[]> mCapFbonds = new ArrayList<int[]>();
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
                    else if (edgeColor == mCapFBond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capF1Bond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapF1Bond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capF2Bond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapF2Bond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFACBond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapFACBond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFBCBond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapFBCBond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFCABond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapFCABond) {
                        mCapFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == capFCBBond) {
                        capFbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == mCapFCBBond) {
                        mCapFbonds.add(new int[]{n1,n2});
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
        mBond = 'm';// multi-body
        lnfXi = new HashSet<Graph>();
        IsoFree isoFree = new IsoFree();

        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();
        
        Metadata.COLORS.add("salmon");
        Metadata.COLORS.add("moccasin");
        Metadata.COLORS.add("steelblue");
        Metadata.COLORS.add("springgreen");
        Metadata.COLORS.add("lime");
        Metadata.COLORS.add("deeppink");
        Metadata.COLORS.add("gold");
        Metadata.COLORS.add("seagreen");
        Metadata.COLORS.add("mintcream");

        fRBond = 'B';
        capFBond = 'U';
        capF1Bond = 'W';
        capF2Bond = 'V';
        capFACBond = 'F';
        capFCABond = 'G';
        capFBCBond = 'I';
        capFCBBond = 'K';
        if (multibody){
        	mERBond = 'P';
	        mCapFBond = 'X';
	        mCapF1Bond = 'Y';
	        mCapF2Bond = 'Z';
	        mCapFACBond = 'L';
	        mCapFCABond = 'M';
	        mCapFBCBond = 'N';
	        mCapFCBBond = 'O';
        }
       
        ArrayList<Character> associationBondAC = new ArrayList<Character>();

        associationBondAC.add(capFACBond);
        associationBondAC.add(capFCABond);
        MetadataImpl.edgeColorPairs.add(associationBondAC);
        ArrayList<Character> associationBondBC = new ArrayList<Character>();
        associationBondBC.add(capFBCBond);
        associationBondBC.add(capFCBBond);
        MetadataImpl.edgeColorPairs.add(associationBondBC);

        
        if (multibody){
            ArrayList<Character> mAssociationBondAC = new ArrayList<Character>();
            mAssociationBondAC.add(mCapFACBond);
            mAssociationBondAC.add(mCapFCABond);
            MetadataImpl.edgeColorPairs.add(mAssociationBondAC);
            ArrayList<Character> mAssociationBondBC = new ArrayList<Character>();
            mAssociationBondBC.add(mCapFBCBond);
            mAssociationBondBC.add(mCapFCBBond);
            MetadataImpl.edgeColorPairs.add(mAssociationBondBC);
        }
       
        if (doShortcut) {
            // just take lnfXi to be the set of connected diagrams
            fBond = 'A';
            eBond = 'e';

            colorOrderMap.put(oneBond, 0);
            colorOrderMap.put(mBond, 1);
            colorOrderMap.put(eBond, 2);
            colorOrderMap.put(fBond, 3);
            colorOrderMap.put(fRBond, 4);
            colorOrderMap.put(capFACBond, 5);
            colorOrderMap.put(capFCABond, 6);
            colorOrderMap.put(capFBCBond, 7);
            colorOrderMap.put(capFCBBond, 8);
            colorOrderMap.put(capFBond, 9);
            colorOrderMap.put(capF1Bond, 10);
            colorOrderMap.put(capF2Bond, 11);
            if (multibody){
	            colorOrderMap.put(mCapFACBond, 12);
	            colorOrderMap.put(mCapFCABond, 13);
	            colorOrderMap.put(mCapFBCBond, 14);
	            colorOrderMap.put(mCapFCBBond, 15);
	            colorOrderMap.put(mCapFBond, 16);
	            colorOrderMap.put(mCapF1Bond, 17);
	            colorOrderMap.put(mCapF2Bond, 18);
            }

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
            colorOrderMap.put(capFACBond, 5);
            colorOrderMap.put(capFCABond, 6);
            colorOrderMap.put(capFBCBond, 7);
            colorOrderMap.put(capFCBBond, 8);
            colorOrderMap.put(capFBond, 9);
            colorOrderMap.put(capF1Bond, 10);
            colorOrderMap.put(capF2Bond, 11);
            if (multibody){
	            colorOrderMap.put(mCapFACBond, 12);
	            colorOrderMap.put(mCapFCABond, 13);
	            colorOrderMap.put(mCapFBCBond, 14);
	            colorOrderMap.put(mCapFCBBond, 15);
	            colorOrderMap.put(mCapFBond, 16);
	            colorOrderMap.put(mCapF1Bond, 17);
	            colorOrderMap.put(mCapF2Bond, 18);
            }
            for (byte i=1; i<n+1; i++) {//1 means 1 point
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
//                System.out.println("Xi with e bonds");
//                topSet.addAll(eXi);
//                for (Graph g : topSet) {
//                    System.out.println(g);
//                }
                System.out.println("eXi has "+eXi.size()+" graphs");
//                ClusterViewer.createView("Xi with e bonds", topSet);
            }
   
           
            Split split = new Split();
            SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
            Set<Graph> setOfSubstituted = split.apply(eXi, bonds);
            eXi.clear();

            DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
            DeleteEdge deleteEdge = new DeleteEdge();
            //set of full star diagrams with f bonds
            Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
            setOfSubstituted.clear();

            if (isInteractive) {
//                topSet.clear();
//                topSet.addAll(fXi);
//                System.out.println("\nXi with f bonds");
//                for (Graph g : topSet) {
//                    System.out.println(g);
//                }
//                topSet.clear();
                System.out.println("fXi has "+fXi.size()+" graphs");
                //ClusterViewer.createView("fXi", topSet);
            }

            Set<Graph> fXipow = new HashSet<Graph>();
            fXipow.addAll(fXi);
            for (int i=1; i<n+1; i++) {
                lnfXi.addAll(fXipow);
                if (i>1) lnfXi = isoFree.apply(lnfXi, null);
                msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
                fXipow = mulScalar.apply(mulFlex.apply(fXipow, fXi, mfp), msp);
            }
            fXipow = null;
            fXi.clear();
        }
        if (isInteractive) {
            topSet.clear();
            topSet.addAll(lnfXi);
//            System.out.println("\nlnfXi");
//            for (Graph g : topSet) {
//                System.out.println(g);
//            }
            System.out.println("lnfXi has "+lnfXi.size()+" graphs");
            //ClusterViewer.createView("lnfXi", topSet);
        }
        //for GCPM water. f = fR+FAC+FBC+FCA+FCB
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(fBond, fRBond, capFBond);
        lnfXi = split.apply(lnfXi, bonds);
        bonds = new SplitParameters(capFBond, capF1Bond, capF2Bond);
        lnfXi = split.apply(lnfXi, bonds);
        bonds = new SplitParameters(capF1Bond, capFACBond, capFCABond);
        lnfXi = split.apply(lnfXi, bonds);
        bonds = new SplitParameters(capF2Bond, capFBCBond, capFCBBond);
        lnfXi = split.apply(lnfXi, bonds);

        if (multibody){//m-bond decomposition
        	bonds = new SplitParameters(mBond,mERBond,mCapFBond);
        	lnfXi = split.apply(lnfXi, bonds);
            bonds = new SplitParameters(mCapFBond,mCapF1Bond,mCapF2Bond);
            lnfXi = split.apply(lnfXi, bonds);
            bonds = new SplitParameters(mCapF1Bond,mCapFACBond,mCapFCABond);
            lnfXi = split.apply(lnfXi, bonds);
            bonds = new SplitParameters(mCapF2Bond,mCapFBCBond,mCapFCBBond);
            lnfXi = split.apply(lnfXi, bonds);
        }
        
        DifByNode opzdlnXidz = new DifByNode();//z*dlnXi/dz
        DifParameters difParams = new DifParameters(nodeColor);
        rho = opzdlnXidz.apply(lnfXi, difParams);
       

        rho0 = new HashSet<Graph>();
        rhoA = new HashSet<Graph>();
        rhoB = new HashSet<Graph>();
        rhoC = new HashSet<Graph>();
        rhoAB = new HashSet<Graph>();
        rhoAC = new HashSet<Graph>();
        rhoBC = new HashSet<Graph>();
        rhoABC = new HashSet<Graph>();
        for (Graph g : rho) {
         boolean hascapFACBond = false;
         boolean hascapFBCBond = false;
         boolean hascapFCABond = false;
         boolean hascapFCBBond = false;
         byte rootNode = 0;
         for (byte i = 0; i<g.nodeCount(); i++){
          if (g.getNode(i).getType()== Metadata.TYPE_NODE_ROOT){
           rootNode = i;
           break;
          }
          
         }

         for (int i = 0;i<g.nodeCount(); i++){
          if (i == rootNode || !g.hasEdge(rootNode, (byte)i))continue;
          if ((g.getEdge(rootNode, (byte)i).getColor() == capFACBond)||(g.getEdge(rootNode, (byte)i).getColor() == mCapFACBond))hascapFACBond = true;
          if ((g.getEdge(rootNode, (byte)i).getColor() == capFBCBond)||(g.getEdge(rootNode, (byte)i).getColor() == mCapFBCBond))hascapFBCBond = true;
          if ((g.getEdge(rootNode, (byte)i).getColor() == capFCABond)||(g.getEdge(rootNode, (byte)i).getColor() == mCapFCABond))hascapFCABond = true;
          if ((g.getEdge(rootNode, (byte)i).getColor() == capFCBBond)||(g.getEdge(rootNode, (byte)i).getColor() == mCapFCBBond))hascapFCBBond = true;
         }
         if (hascapFACBond&&hascapFBCBond&&(hascapFCABond||hascapFCBBond))
          rhoABC.add(g);
         else if (hascapFACBond&&(hascapFCABond||hascapFCBBond))
          rhoAC.add(g);
         else if (hascapFBCBond&&(hascapFCABond||hascapFCBBond))
          rhoBC.add(g);
         else if (hascapFACBond&&hascapFBCBond)
          rhoAB.add(g);
         else if (hascapFACBond)
          rhoA.add(g);
         else if (hascapFBCBond)
          rhoB.add(g);
         else if (hascapFCABond||hascapFCBBond)
          rhoC.add(g);
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
       
        Set<Graph> newRhoC = new HashSet<Graph>();
        for (Graph g : rhoC){
         byte rootId = 0;
         for (Node node:g.nodes()){
          if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
         }
         if (rootId == 0){
	          newRhoC.add(g);
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
            newRhoC.add(g);
        }
        rhoC.clear();
        rhoC.addAll(newRhoC);
       
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
       
        Set<Graph> newRhoAC = new HashSet<Graph>();
        for (Graph g : rhoAC){
         byte rootId = 0;
         for (Node node:g.nodes()){
          if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
         }
         if (rootId == 0){
	          newRhoAC.add(g);
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
            newRhoAC.add(g);
        }
        rhoAC.clear();
        rhoAC.addAll(newRhoAC);
       
        Set<Graph> newRhoBC = new HashSet<Graph>();
        for (Graph g : rhoBC){
         byte rootId = 0;
         for (Node node:g.nodes()){
          if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
         }
         if (rootId == 0){
	          newRhoBC.add(g);
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
            newRhoBC.add(g);
        }
        rhoBC.clear();
        rhoBC.addAll(newRhoBC);
       
        Set<Graph> newRhoABC = new HashSet<Graph>();
        for (Graph g : rhoABC){
         byte rootId = 0;
         for (Node node:g.nodes()){
          if (node.getType() == Metadata.TYPE_NODE_ROOT) rootId = node.getId();
         }
         if (rootId == 0){
	          newRhoABC.add(g);
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
            newRhoABC.add(g);
        }
        rhoABC.clear();
        rhoABC.addAll(newRhoABC);
       
        Set<Graph> rho0pow = new HashSet<Graph>();
        Set<Graph> rho0m1 = new HashSet<Graph>();

        for (Graph g : rho0){//every graph in rho0
	         if (g.nodeCount() > 1)rho0m1.add(g);//rho0m1 = rho0 - single pt
        }

        rho0pow.addAll(rho0m1);
        mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)(n-1));
       
        Set<Graph> rho0m1pow = new HashSet<Graph>();

        msp = new MulScalarParameters(new CoefficientImpl(-1,1));
        rho0m1 = mulScalar.apply(rho0m1, msp);
        rho0m1pow.addAll(rho0m1);
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
            if (isInteractive) {
//                topSet.clear();
//                topSet.addAll(zWertheim);
//                System.out.println("\nzWertheim");
//                for (Graph g : topSet) {
//                    System.out.println(g);
//                }
                System.out.println("zWertheim has "+zWertheim.size()+" graphs");
                //ClusterViewer.createView("p(rho0)", topSet);
            }
           
            Decorate decorate = new Decorate();
            DecorateParameters dp = new DecorateParameters(nodeColor, mfp);
            DecorateParameters dpnm1 = new DecorateParameters(0,mfpnm1zWertheim);//0th factor

            p = decorate.apply(lnfXi, zWertheim, dp);
            lnfXi.clear();
            p = isoFree.apply(p, null);
            
           
            if (isInteractive) {
//                topSet.clear();
//                topSet.addAll(p);
//                System.out.println("\np(rho0)");
//                for (Graph g : topSet) {
//                    System.out.println(g);
//                }
                System.out.println("p(rho0) has "+p.size()+" graphs");
                //ClusterViewer.createView("p(rho0)", topSet);
            }

          
            for(Graph g:rhoA){
	                g.setNumFactors(1);
	                g.addFactors(new int[]{g.nodeCount()});
	               }
               for(Graph g:rhoB){
	                g.setNumFactors(1);
	                g.addFactors(new int[]{g.nodeCount()});
	               }
               for(Graph g:rhoC){
	                g.setNumFactors(1);
	                g.addFactors(new int[]{g.nodeCount()});
	               }
               for(Graph g:rhoAB){
	               g.setNumFactors(1);
	               g.addFactors(new int[]{g.nodeCount()});
	              }
               for(Graph g:rhoAC){
                   g.setNumFactors(1);
                   g.addFactors(new int[]{g.nodeCount()});
                  }
               for(Graph g:rhoBC){
                   g.setNumFactors(1);
                   g.addFactors(new int[]{g.nodeCount()});
                  }
               for(Graph g:rhoABC){
                  g.setNumFactors(1);
                  g.addFactors(new int[]{g.nodeCount()});
                 }
               for(Graph g:zWertheim){
                  g.setNumFactors(1);
                 }
                  
              rhoA = isoFree.apply(decorate.apply(rhoA, zWertheim, dpnm1),null);
              rhoB = isoFree.apply(decorate.apply(rhoB, zWertheim, dpnm1),null);
              rhoC = isoFree.apply(decorate.apply(rhoC, zWertheim, dpnm1),null);
              rhoAB = isoFree.apply(decorate.apply(rhoAB, zWertheim, dpnm1),null);
              rhoAC = isoFree.apply(decorate.apply(rhoAC, zWertheim, dpnm1),null);
              rhoBC = isoFree.apply(decorate.apply(rhoBC, zWertheim, dpnm1),null);
              rhoABC = isoFree.apply(decorate.apply(rhoABC, zWertheim, dpnm1),null);
              zWertheim.clear();
	          for(Graph g:rhoA){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoB){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoC){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoAB){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoAC){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoBC){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
	          for(Graph g:rhoABC){
	              g.setNumFactors(2);
	              g.addFactors(new int[]{g.nodeCount(),0});
	             }
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Factor factor = new Factor();

            DecorateWertheim3SiteRho decorateWertheim3Site = new DecorateWertheim3SiteRho();
            DecorateWertheimParameters3Site dpWertheim3Site = new DecorateWertheimParameters3Site(mfp, capFACBond, capFBCBond, capFCABond, capFCBBond,mCapFACBond, mCapFBCBond, mCapFCABond, mCapFCBBond,rhoA,rhoB,rhoC,rhoAB,rhoAC,rhoBC,rhoABC);//decorate rho point with rho0C1
            Set<Graph> pWertheim3Site = decorateWertheim3Site.apply(p, dpWertheim3Site);
            
            rhoA.clear();
            rhoB.clear();
            rhoC.clear();
            rhoAB.clear();
            rhoAC.clear();
            rhoBC.clear();
            rhoABC.clear();
 
            pWertheim3Site = isoFree.apply(pWertheim3Site, null);

            MaxIsomorphWertheim maxIso = new MaxIsomorphWertheim();
            List<Character> aBonds = new ArrayList<Character>();
            List<Character> bBonds = new ArrayList<Character>();
            aBonds.add(capFACBond);
            aBonds.add(mCapFACBond);
            bBonds.add(capFBCBond);
            bBonds.add(mCapFBCBond);
            
            
            BondConnectedSubstitution3Site bondConnected = new BondConnectedSubstitution3Site();
            BondConnectedSubstitutionParameters3Site bcp3Site = new BondConnectedSubstitutionParameters3Site(fRBond, eBond); 

            pWertheim3Site = bondConnected.apply(pWertheim3Site,bcp3Site);

            //combine graphs
            Set<Graph> newpWertheim3Site = new HashSet<Graph>();
            for (Graph g : pWertheim3Site) {
                boolean ap = hap.check(g);//graph with a.p.
                boolean con = hap.isConnected();//graph connected at least singly
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    Graph gf = factor.apply(g, mfp);//gf = disconnected diagram
                    newpWertheim3Site.add(gf);
                }
                else {
                    newpWertheim3Site.add(g.copy());
                }
            }
            pWertheim3Site = newpWertheim3Site;
            
            CombineABSiteParameters abParams = new CombineABSiteParameters(aBonds, bBonds); 
            pWertheim3Site = maxIso.apply(pWertheim3Site, new MaxIsomorphWertheimParameters(abParams));

            for (Graph g:pWertheim3Site){
                int[] factors = g.factors();
                int[] newFactors = new int[6];
                newFactors[0] = factors[0];
                newFactors[1] = factors[1]+factors[2];
                newFactors[2] = factors[3];
                newFactors[3] = factors[4];
                newFactors[4] = factors[5]+factors[6];
                newFactors[5] = factors[7];
                g.setNumFactors(6);
                g.addFactors(newFactors);
            }

            MetadataImpl.rootPointsSpecial = false;// every root point is indistinguishable
            pWertheim3Site = isoFree.apply(pWertheim3Site, null);

            if (isInteractive) {
                topSet.clear();
                topSet.addAll(pWertheim3Site);
                System.out.println("\npwertheim");
                for (Graph g : topSet) {
                    System.out.println(g);
                }
                //ClusterViewer.createView("pwertheim", topSet);
                System.out.println("pwertheim has "+pWertheim3Site.size()+" graphs");
            }
            writeSetToFile(pWertheim3Site, "pw"+n);

            // clear these out -- we don't need them and (in extreme cases) we might need the memory
            z.clear();
        }
    }

    protected static void writeSetToFile(Set<Graph> set, String fileName) {
        try {
            FileWriter fw = new FileWriter(fileName);
            Set<Graph> topSet = makeGraphList();
            topSet.addAll(set);
            for (Graph g : topSet) {
                fw.write(g.toString()+"\n");
            }
            fw.close();
        }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}