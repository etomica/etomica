/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNodeColors;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.Decorate;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.FactorOnce;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;

public class VirialDiagramsMix {

    protected final int n;
    protected char[] nodeColors;
    protected final boolean[] flex;
    protected final boolean[] multibody;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Set<Graph> multiP;
    protected Set<Graph> rhoA, rhoB;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Graph> cancelMap;
    protected boolean doShortcut;
    protected boolean doMinimalMulti;
    protected boolean doMinimalBC;
    protected boolean doKeepEBonds;
    protected final char nodeColor = Metadata.COLOR_CODE_0;
    protected char[] flexColors;
    public char fBond, eBond, excBond, mBond, MBond, efbcBond;
    
    public static void main(String[] args) {
        int n = 4;
        boolean[] multibody = new boolean[]{false,false};
        boolean[] flex = new boolean[]{true,true};
        boolean doKeepEBonds = false;
        boolean doReeHoover = true;
        VirialDiagramsMix virialDiagrams = new VirialDiagramsMix(n, multibody, flex, true);
        virialDiagrams.setDoReeHoover(doReeHoover);
        virialDiagrams.setDoKeepEBonds(doKeepEBonds);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.setDoMinimalMulti(true);
        virialDiagrams.setDoMinimalBC(true);
        virialDiagrams.makeVirialDiagrams();
    }
    
    public VirialDiagramsMix(int n, boolean[] multibody, boolean[] flex) {
        this(n, multibody, flex, false);
    }

    public VirialDiagramsMix(int n, boolean[] multibody, boolean[] flex, boolean interactive) {
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

    public void setDoMinimalMulti(boolean newDoMinimalMulti) {
        doMinimalMulti = newDoMinimalMulti;
    }

    public void setDoMinimalBC(boolean newDoMinimalBC) {
        doMinimalBC = newDoMinimalBC;
    }

    public Set<Graph> makeGraphList() {
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNodeColors(nodeColors));
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        return new GraphList(comp);
    }
    
    public void makeRhoDiagrams() {

        final char nodeA = Metadata.COLOR_CODE_0;
        final char nodeB = Metadata.COLOR_CODE_1;
        // we'll pretend that everything is flexible until the end
        // if we allow rigid multiplication to happen during intermediate
        // steps, we get confused because multiplications happen in an order
        // that makes things unhappy (root points in the "wrong" place).  We
        // could work around this by having multiplication move root points
        // around to an appropriate color, but that seems icky.
        nodeColors = new char[]{nodeA,nodeB};
        flexColors = new char[]{};
        if (flex[0]) {
            flexColors = new char[]{nodeA};
            if (flex[1]) {
                flexColors = new char[]{nodeA,nodeB};
            }
        }
        else if (flex[1]) {
            flexColors = new char[]{nodeB};
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

        Metadata.COLOR_MAP.put(eBond, "red");
        Metadata.COLOR_MAP.put(fBond, "green");
        Metadata.COLOR_MAP.put(mBond, "blue");
        Metadata.COLOR_MAP.put(MBond, "orange");
        Metadata.COLOR_MAP.put(efbcBond, "fuchsia");
        Metadata.COLOR_MAP.put(excBond, "red");
        Metadata.DASH_MAP.put(excBond, 3);

        
        Set<Graph> topSet = makeGraphList();
        
        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
        System.out.println("Xi");
        // factors: zA, zB, rhoA, rhoB
        for (byte i=1; i<n+1; i++) {
            for (byte j=0; j<i+1; j++) {
                // j points of color A, i-j points of color B
                Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                g.coefficient().setDenominator((int)(SpecialFunctions.factorial(i-j)*SpecialFunctions.factorial(j)));
                g.setNumFactors(4);
                g.addFactors(new int[]{j,i-j,0,0});
                for (byte k=j; k<i; k++) {
                    g.getNode(k).setColor(nodeB);
                }
                for (Node node1 : g.nodes()) {
                    for (Node node2 : g.nodes()) {
                        if (node2.getId() <= node1.getId()) continue;
                        g.getEdge(node1.getId(), node2.getId()).setColor(eBond);
                    }
                }
                eXi.add(g);
                
/*                if ((multiA && j > 2) || (multiB & i-j > 2) || (multiA && multiB && i > 2)) { 
                    g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                    g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                    g.setNumFactors(4);
                    g.addFactors(new int[]{j,i-j,0,0});
                    for (byte k=j; k<i; k++) {
                        g.getNode(k).setColor(nodeB);
                    }
                    for (Node node1 : g.nodes()) {
                        for (Node node2 : g.nodes()) {
                            if (node2.getId() <= node1.getId()) continue;
                            if (node1.getColor() == nodeA) {
                                if (node2.getColor() == nodeA) {
                                    g.getEdge(node1.getId(), node2.getId()).setColor(multiA ? mBond : eBondAA);
                                }
                                else {
                                    g.getEdge(node1.getId(), node2.getId()).setColor((multiA && multiB) ? mBond : eBondAB);
                                }
                            }
                            else if (node2.getColor() == nodeA) {
                                g.getEdge(node1.getId(), node2.getId()).setColor((multiA && multiB) ? mBond : eBondAB);
                            }
                            else {
                                g.getEdge(node1.getId(), node2.getId()).setColor(multiB ? mBond : eBondBB);
                            }
                        }
                    }

                    eXi.add(g);
                }*/
            }
        }
        topSet.addAll(eXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
//        ClusterViewer.createView("eXi", topSet);

        Split split = new Split();
        SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
        Set<Graph> setOfSubstituted = split.apply(eXi, bonds);

        DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
        DeleteEdge deleteEdge = new DeleteEdge();
        //set of full star diagrams with f bonds
        System.out.println("\nXi with f bonds");
        Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
        IsoFree isoFree = new IsoFree();
        fXi = isoFree.apply(fXi, null);
//        ClusterViewer.createView("fXi", fXi);
        topSet.clear();
        topSet.addAll(fXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
        
        MulFlexible mulFlex = new MulFlexible();
        
        MulFlexibleParameters mfpn = MulFlexibleParameters.makeParameters(nodeColors, (byte)n);
        lnfXi = new HashSet<Graph>();
        Set<Graph> fXipow = new HashSet<Graph>();
        fXipow.addAll(fXi);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();
        for (int i=1; i<n+1; i++) {

            lnfXi.addAll(fXipow);
            lnfXi = isoFree.apply(lnfXi, null);
            msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
            fXipow = isoFree.apply(mulScalar.apply(mulFlex.apply(fXipow, fXi, mfpn), msp), null);
        }
        topSet.clear();
        topSet.addAll(lnfXi);
        System.out.println("\nlnfXi");
        for (Graph g : topSet) {
            System.out.println(g);
        }

        DifByNode opzdlnXidzA = new DifByNode();
        DifParameters difParams = new DifParameters('A');
        rhoA = isoFree.apply(opzdlnXidzA.apply(lnfXi, difParams), null);
        
        topSet.clear();
        topSet.addAll(rhoA);
        ClusterViewer.createView("rhoA", topSet);
        System.out.println("\nrhoA");
        for (Graph g : topSet) {
            System.out.println(g);
        }

        DifByNode opzdlnXidzB = new DifByNode();
        difParams = new DifParameters('B');
        rhoB = opzdlnXidzB.apply(lnfXi, difParams);

        rhoB = isoFree.apply(rhoB, null);

        topSet.clear();
        topSet.addAll(rhoB);
        ClusterViewer.createView("rhoB", topSet);
        System.out.println("\nrhoB");
        for (Graph g : topSet) {
            System.out.println(g);
        }
    }
    
    public void makeVirialDiagrams() {
        if (rhoA == null) {
            makeRhoDiagrams();
        }

        MulScalarParameters msp = new MulScalarParameters(new CoefficientImpl(-1,1));
        MulScalar mulScalar = new MulScalar();
        Set<Graph> zAz = new HashSet<Graph>();
        Set<Graph> zBz = new HashSet<Graph>();
        Set<Graph> zA = new HashSet<Graph>();
        Set<Graph> zB = new HashSet<Graph>();
        for (Graph g : rhoA) {
            if (g.nodeCount() == 1) {
                // switch zA to rhoA
                g = g.copy();
                g.factors()[0] = 0;
                g.factors()[2] = 1;
                // zA = rhoA is our initial approximation for zA in terms of rho
                zA.add(g.copy());
            }
            else {
                // each of our terms in zAz is the term from rhoA, but subtracted
                g = mulScalar.apply(g, msp);
            }
            zAz.add(g);
        }
        for (Graph g : rhoB) {
            if (g.nodeCount() == 1) {
                g = g.copy();
                g.factors()[1] = 0;
                g.factors()[3] = 1;
                zB.add(g.copy());
            }
            else {
                g = mulScalar.apply(g, msp);
            }
            zBz.add(g);
        }
        // we have  zAz = rhoA - a20 zA^2 - a11 zAzB
                
        Decorate decorate = new Decorate();
        MulFlexibleParameters mfpnm1 = MulFlexibleParameters.makeParameters(nodeColors, (byte)(n-1));
        for (int i=2; i<n+1; i++) {
            // now decorate zAz with zA and zB
            // we actually only need zAz to ith order, but that's more work.  Decorate will truncate for us.
            Set<Graph> newZA = decorate.apply(zAz, zA, new DecorateParameters(0, mfpnm1));
            newZA = decorate.apply(newZA, zB, new DecorateParameters(1, mfpnm1));
            Set<Graph> newZB = decorate.apply(zBz, zA, new DecorateParameters(0, mfpnm1));
            newZB = decorate.apply(newZB, zB, new DecorateParameters(1, mfpnm1));
            
            zA = newZA;
            zB = newZB;
        }
        
        IsoFree isoFree = new IsoFree();
        zA = isoFree.apply(zA, null);
        zB = isoFree.apply(zB, null);

        Set<Graph> topSet = makeGraphList();
        topSet.clear();
        topSet.addAll(zA);
        System.out.println("\nzA");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("zA", topSet);

        topSet.clear();
        topSet.addAll(zB);
        System.out.println("\nzB");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("zB", topSet);

        MulFlexibleParameters mfpn = MulFlexibleParameters.makeParameters(nodeColors, (byte)n);
        p = decorate.apply(lnfXi, zA, new DecorateParameters(0, mfpn));
        p = decorate.apply(p, zB, new DecorateParameters(1, mfpn));
        p = isoFree.apply(p, null);

        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        if (flexColors.length < nodeColors.length) {
            Set<Graph> newP = new HashSet<Graph>();
            Factor factor = new Factor();
            MulFlexibleParameters factorParameters = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
                    boolean factorable = false;
                    for (byte nodeID : hap.getArticulationPoints()) {
                        char color = g.getNode(nodeID).getColor();
                        factorable = true;
                        for (int i=0; i<flexColors.length; i++) {
                            if (flexColors[i] == color) {
                                factorable = false;
                                break;
                            }
                        }
                        if (factorable) break;
                    }
                    if (factorable) {
                        Graph gf = factor.apply(g, factorParameters);
                        newP.add(gf);
                    }
                    else {
                        // graph had an articulation point, but it was flexible
                        newP.add(g.copy());
                    }
                }
                else {
                    newP.add(g.copy());
                }
            }
            p = isoFree.apply(newP, null);
        }
        HashMap<Graph,Set<Graph>> cancelSet = new HashMap<Graph,Set<Graph>>();
        if (flexColors.length > 0) {
            // pretend everything is fully flexible
            FactorOnce factor = new FactorOnce();
            Relabel relabel = new Relabel();
            FactorOnceParameters fopA = new FactorOnceParameters((byte)0, false);
            FactorOnceParameters fopB = new FactorOnceParameters((byte)1, false);
            Set<Graph> newP = new HashSet<Graph>();
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if (con && ap) {
                    FactorOnceParameters fop = null;
                    boolean factorable0 = g.getNode((byte)0).getColor() == nodeColors[0] && hap.getArticulationPoints().contains(0);
                    if (!factorable0) {
                        byte articulationId = -1;
                        for (byte nodeId : hap.getArticulationPoints()) {
                            if (g.getNode(nodeId).getColor() == nodeColors[0]) {
                                articulationId = nodeId;
                                break;
                            }
                        }
                        if (articulationId != -1) {
                            factorable0 = true;
                            byte[] permutations = new byte[g.nodeCount()];
                            for (int i=0; i<permutations.length; i++) {
                                permutations[i] = (byte)i;
                            }
                            permutations[0] = articulationId;
                            permutations[articulationId] = 0;
                            RelabelParameters rp = new RelabelParameters(permutations);
                            g = relabel.apply(g, rp);
                        }
                    }
                    boolean factorableB = g.getNode((byte)1).getColor() == nodeColors[1] && hap.getArticulationPoints().contains(1);
                    if (!factorable0 && !factorableB) {
                        byte articulationId = -1;
                        for (byte nodeId : hap.getArticulationPoints()) {
                            if (g.getNode(nodeId).getColor() == nodeColors[1]) {
                                articulationId = nodeId;
                                break;
                            }
                        }
                        if (articulationId != -1) {
                            factorableB = true;
                            byte[] permutations = new byte[g.nodeCount()];
                            for (int i=0; i<permutations.length; i++) {
                                permutations[i] = (byte)i;
                            }
                            permutations[1] = articulationId;
                            permutations[articulationId] = 1;
                            RelabelParameters rp = new RelabelParameters(permutations);
                            g = relabel.apply(g, rp);
                        }
                    }
                    g = g.copy();
                    newP.add(g);
                    if (factorable0) {
                        fop = fopA;
                    }
                    else if (factorableB) {
                        fop = fopB;
                    }
                    else {
                        throw new RuntimeException("oops");
                    }
                    Set<Graph> gf = factor.apply(g, fop);
                    cancelSet.put(g, gf);
                    // add a copy here so that it can be modified by IsoMorphism filter without
                    // adversely affecting our cancelSet graph
                    newP.addAll(gf);
                }
                else {
                    newP.add(g.copy());
                }
            }
            // we have to be extra careful because we made a graph-graph hasmap,
            // so we can't make any graph copies
            IteratorWrapper wrapper = new IteratorWrapper(newP.iterator());
            GraphIterator isomorphs = new IsomorphismFilter(wrapper);
            p = new HashSet<Graph>();
            while (isomorphs.hasNext()) {
                p.add(isomorphs.next());
            }
        }

        topSet.clear();
        topSet.addAll(p);
        System.out.println("\nP");
        for (Graph g : topSet) {
            System.out.println(g);
            Set<Graph> cancelGraph = cancelSet.get(g);
            if (cancelGraph != null) {
                for (Graph cg : cancelGraph) {
                    System.out.println(" -  "+cg);
                }
            }
        }
        ClusterViewer.createView("P", topSet);
    }
}
