package etomica.virial.cluster;

import java.util.HashMap;
import java.util.HashSet;
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

    public static void main(String[] args) {
        final int n = 3;
        final char nodeA = Metadata.COLOR_CODE_0;
        final char nodeB = Metadata.COLOR_CODE_1;
        // we'll pretend that everything is flexible until the end
        // if we allow rigid multiplication to happen during intermediate
        // steps, we get confused because multiplications happen in an order
        // that makes things unhappy (root points in the "wrong" place).  We
        // could work around this by having multiplication move root points
        // around to an appropriate color, but that seems icky.
        char[] allColors = new char[]{nodeA,nodeB};
        char[] flexColors = new char[]{nodeA, nodeB};
        boolean multiA = false;
        boolean multiB = false;
        
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNodeColors(new char[]{nodeA,nodeB}));
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList<Graph> topSet = new GraphList<Graph>(comp);
        
        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
        System.out.println("Xi");
        char eBondAA = Metadata.COLOR_CODE_0;
        char eBondAB = Metadata.COLOR_CODE_1;
        char eBondBB = Metadata.COLOR_CODE_2;
        char mBond = Metadata.COLOR_CODE_3;
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
                        if (node1.getColor() == nodeA) {
                            if (node2.getColor() == nodeB) {
                                g.getEdge(node1.getId(), node2.getId()).setColor(eBondAB);
                            }
                        }
                        else if (node2.getColor() == nodeA) {
                            g.getEdge(node1.getId(), node2.getId()).setColor(eBondAB);
                        }
                        else {
                            g.getEdge(node1.getId(), node2.getId()).setColor(eBondBB);
                        }
                    }
                }
                eXi.add(g);
                
                if ((multiA && j > 2) || (multiB & i-j > 2) || (multiA && multiB && i > 2)) { 
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
                }
            }
        }
        topSet.addAll(eXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
//        ClusterViewer.createView("eXi", topSet);

        char fBondAA = Metadata.COLOR_CODE_4;
        char fBondAB = Metadata.COLOR_CODE_5;
        char fBondBB = Metadata.COLOR_CODE_6;
        char oneBond = 'o';
        
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(eBondAA, fBondAA, oneBond);
        Set<Graph> setOfSubstituted = split.apply(eXi, bonds);
        bonds = new SplitParameters(eBondAB, fBondAB, oneBond);
        setOfSubstituted = split.apply(setOfSubstituted, bonds);
        bonds = new SplitParameters(eBondBB, fBondBB, oneBond);
        setOfSubstituted = split.apply(setOfSubstituted, bonds);

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
        
        MulFlexibleParameters mfpn = new MulFlexibleParameters(allColors, (byte)n);
        MulFlexibleParameters mfpnm1 = new MulFlexibleParameters(allColors, (byte)(n-1));
        Set<Graph> lnfXi = new HashSet<Graph>();
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
        Set<Graph> rhoA = isoFree.apply(opzdlnXidzA.apply(lnfXi, difParams), null);
        
        topSet.clear();
        topSet.addAll(rhoA);
        ClusterViewer.createView("rhoA", topSet);
        System.out.println("\nrhoA");
        for (Graph g : topSet) {
            System.out.println(g);
        }

        DifByNode opzdlnXidzB = new DifByNode();
        difParams = new DifParameters('B');
        Set<Graph> rhoB = opzdlnXidzB.apply(lnfXi, difParams);

        rhoB = isoFree.apply(rhoB, null);

        topSet.clear();
        topSet.addAll(rhoB);
        ClusterViewer.createView("rhoB", topSet);
        System.out.println("\nrhoB");
        for (Graph g : topSet) {
            System.out.println(g);
        }

        msp = new MulScalarParameters(new CoefficientImpl(-1,1));
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
        
        zA = isoFree.apply(zA, null);
        zB = isoFree.apply(zB, null);

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

        Set<Graph> p = decorate.apply(lnfXi, zA, new DecorateParameters(0, mfpn));
        p = decorate.apply(p, zB, new DecorateParameters(1, mfpn));
        p = isoFree.apply(p, null);

        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        if (flexColors.length < allColors.length) {
            Set<Graph> newP = new HashSet<Graph>();
            Factor factor = new Factor();
            MulFlexibleParameters factorParameters = new MulFlexibleParameters(flexColors, (byte)n);
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
            FactorOnceParameters fopA = new FactorOnceParameters((byte)0, new char[0], false);
            FactorOnceParameters fopB = new FactorOnceParameters((byte)1, new char[0], false);
            Set<Graph> newP = new HashSet<Graph>();
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if (con && ap) {
                    FactorOnceParameters fop = null;
                    boolean factorable0 = g.getNode((byte)0).getColor() == nodeA && hap.getArticulationPoints().contains(0);
                    if (!factorable0) {
                        byte articulationId = -1;
                        for (byte nodeId : hap.getArticulationPoints()) {
                            if (g.getNode(nodeId).getColor() == nodeA) {
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
                    boolean factorableB = g.getNode((byte)1).getColor() == nodeB && hap.getArticulationPoints().contains(1);
                    if (!factorable0 && !factorableB) {
                        byte articulationId = -1;
                        for (byte nodeId : hap.getArticulationPoints()) {
                            if (g.getNode(nodeId).getColor() == nodeB) {
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
