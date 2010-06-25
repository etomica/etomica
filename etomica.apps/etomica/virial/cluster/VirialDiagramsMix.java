package etomica.virial.cluster;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
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
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.SetPropertyFilter;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.FieldNodeCountMax;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;

public class VirialDiagramsMix {

    public static void main(String[] args) {
        final int n = 4;
        final char nodeA = Metadata.COLOR_CODE_0;
        final char nodeB = Metadata.COLOR_CODE_1;
        char[] flexColors = new char[]{nodeA,nodeB};
        
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
            }
        }
        topSet.addAll(eXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
//        dump(topSet, n, nodeA);
//        ClusterViewer.createView("eXi", eXi);
        char fBondAA = Metadata.COLOR_CODE_3;
        char fBondAB = Metadata.COLOR_CODE_4;
        char fBondBB = Metadata.COLOR_CODE_5;
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
        
        MulFlexibleParameters mfpn = new MulFlexibleParameters(flexColors, (byte)n);
        MulFlexibleParameters mfpnm1 = new MulFlexibleParameters(flexColors, (byte)(n-1));
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

        // rhoA in powers of zA and zB
        HashSet<Graph>[][] allRhoA = new HashSet[n+1][n+1];
        // rhoB in powers of zA and zB
        HashSet<Graph>[][] allRhoB = new HashSet[n+1][n+1];
        for (int i=0; i<n+1; i++) {
            for (int j=0; j<n+1; j++) {
                allRhoA[i][j] = new HashSet<Graph>();
                allRhoB[i][j] = new HashSet<Graph>();
            }
        }
        for (Graph g : rhoA) {
            int nodeACount = 0;
            for (Node node : g.nodes()) {
                if (node.getColor() == nodeA) {
                    nodeACount++;
                }
            }
            allRhoA[nodeACount][g.nodeCount()-nodeACount].add(g);
        }
        for (Graph g : rhoB) {
            int nodeBCount = 0;
            for (Node node : g.nodes()) {
                if (node.getColor() == nodeB) {
                    nodeBCount++;
                }
            }
            allRhoB[nodeBCount][g.nodeCount()-nodeBCount].add(g);
        }

        Set<Graph> zA = new HashSet<Graph>();
        zA.addAll(allRhoA[1][0]);
        for (Graph g : zA) {
            g.setNumFactors(4);
            g.addFactors(new int[]{0,0,1,0});
        }
        Set<Graph> zB = new HashSet<Graph>();
        zB.addAll(allRhoB[1][0]);
        for (Graph g : zB) {
            g.setNumFactors(4);
            g.addFactors(new int[]{0,0,0,1});
        }
        for (int i=2; i<n+1; i++) {
            Set<Graph>[] zAPow = new Set[n+1];
            zAPow[1] = new HashSet<Graph>();
            zAPow[1].addAll(zA);
            Set<Graph>[] zBPow = new Set[n+1];
            zBPow[1] = new HashSet<Graph>();
            zBPow[1].addAll(zB);
            for (int j=2; j<i+1; j++) {
                zAPow[j] = new HashSet<Graph>();
                zAPow[j] = isoFree.apply(mulFlex.apply(zAPow[j-1], zA, mfpnm1), null);
                zBPow[j] = new HashSet<Graph>();
                zBPow[j] = isoFree.apply(mulFlex.apply(zBPow[j-1], zB, mfpnm1), null);
            }
            zA = new HashSet<Graph>();
            zA.addAll(allRhoA[1][0]);
            zB = new HashSet<Graph>();
            zB.addAll(allRhoB[1][0]);
            msp = new MulScalarParameters(new CoefficientImpl(-1,1));
            for (int j=2; j<i+1; j++) {
                for (int k=1; k<j+1; k++) {
                    Set<Graph> term = mulFlex.apply(allRhoA[k][j-k], zAPow[k], mfpnm1);
                    if (j-k > 0) {
                        term = mulFlex.apply(term, zBPow[j-k], mfpnm1);
                    }
                    zA.addAll(mulScalar.apply(term, msp));

                    term = mulFlex.apply(allRhoB[k][j-k], zBPow[k], mfpnm1);
                    if (j-k > 0) {
                        term = mulFlex.apply(term, zAPow[j-k], mfpnm1);
                    }
                    zB.addAll(mulScalar.apply(term, msp));
                }
            }
            for (Graph g : zA) {
                g.factors()[0] = 0;
                g.factors()[1] = 0;
            }
            for (Graph g : zB) {
                g.factors()[0] = 0;
                g.factors()[1] = 0;
            }
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

        Decorate decorate = new Decorate();
        Set<Graph> p = decorate.apply(lnfXi, zA, new DecorateParameters(0, mfpn));
        p = decorate.apply(p, zB, new DecorateParameters(1, mfpn));
        p = isoFree.apply(p, null);

        topSet.clear();
        topSet.addAll(p);
        System.out.println("\nP");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("P", topSet);
    }
    
    public static class Term {
        public int coefficient;
        public Set<Graph>[] factors;
        public Term(int coefficient, Set<Graph>[] factors) {
            this.coefficient = coefficient;
            this.factors = factors;
        }
    }

    protected static Set<Graph> invertSeries(Set<Graph>[][][] allRho, int n, String comp0, String comp1, MulFlexibleParameters mfp) {
        MulScalar mulScalar = new MulScalar();
        MulFlexible mulFlex = new MulFlexible();
        IsoFree isoFree = new IsoFree();
        SetPropertyFilter truncater = new SetPropertyFilter(new FieldNodeCountMax(n-1));
        Set<Graph> z = new HashSet<Graph>();
        
        int[][][][] terms = new int[][][][]{
                {{{1},{0,1,0}}},
                
                {{{-1},{0,2,0}}},
                {{{-1},{0,1,1}}},
                
                {{{-1},{0,3,0}}, {{2},{0,2,0},{0,2,0}}},
                {{{-1},{0,2,1}}, {{3},{0,1,1},{0,2,0}}, {{1},{0,1,1},{1,1,1}}},
                {{{-1},{0,1,2}}, {{1},{0,1,1},{0,1,1}}, {{1},{0,1,1},{1,2,0}}},

                {{{-5}, {0,2,0}, {0,2,0}, {0,2,0}}, {{5}, {0,2,0}, {0,3,0}}, {{-1}, {0,4,0}}}, 
                {{{4}, {0,2,0}, {0,2,1}}, {{-1}, {0,3,1}}, {{1}, {0,2,1}, {1,1,1}}, {{-10}, {0,1,1}, {0,2,0}, {0,2,0}}, {{4}, {0,1,1}, {0,3,0}}, {{-4}, {0,1,1}, {0,2,0}, {1,1,1}}, {{-1}, {0,1,1}, {1,1,1}, {1,1,1}}, {{1}, {0,1,1}, {1,1,2}}}, 
                {{{-1}, {0,2,2}}, {{1}, {0,2,1}, {1,2,0}}, {{-6}, {0,1,1}, {0,1,1}, {0,2,0}}, {{-3}, {0,1,1}, {0,1,1}, {1,1,1}}, {{3}, {0,1,2}, {0,2,0}}, {{2}, {0,1,2}, {1,1,1}}, {{3}, {0,1,1}, {0,2,1}}, {{-3}, {0,1,1}, {0,2,0}, {1,2,0}}, {{-3}, {0,1,1}, {1,2,0}, {1,1,1}}, {{1}, {0,1,1}, {1,2,1}}}, 
                {{{-1}, {0,1,1}, {0,1,1}, {0,1,1}}, {{-1}, {0,1,3}}, {{-2}, {0,1,1}, {0,1,1}, {1,2,0}}, {{2}, {0,1,2}, {1,2,0}}, {{2}, {0,1,1}, {0,1,2}}, {{-2}, {0,1,1}, {1,2,0}, {1,2,0}}, {{1}, {0,1,1}, {1,3,0}}}
        };
              

        int nT = 1, n1 = -1; 
        for (int k=0; k<terms.length; k++ ) {
            n1++;
            if (n1 == nT) {
                nT++;
                n1 = 0;
            }
            int[] factors;
            if (comp0.equals("A")) {
                factors = new int[]{0,0,nT-n1,n1};
            }
            else {
                factors = new int[]{0,0,n1,nT-n1};
            }
                    
            int[][][] termSet = terms[k];
            for (int i=0; i<termSet.length; i++ ) {
                int[][] term = termSet[i];
                int coefficient = term[0][0];
                int order = 1-term.length;
                Set<Graph> product = new HashSet<Graph>();
                for (int j=1; j<term.length; j++) {
                    order += term[j][1]+term[j][2];
                    if (order > n || term[j][1]+term[j][2] > n) {
                        return z;
                    }
                    if (j == 1) {
                        product.addAll(allRho[term[j][0]][term[j][1]][term[j][2]]);
                    }
                    else {
                        product = mulFlex.apply(product, allRho[term[j][0]][term[j][1]][term[j][2]], mfp);
                        product = isoFree.apply(truncater.apply(product, null), null);
                    }
                }
                product = mulScalar.apply(product, new MulScalarParameters(coefficient,1));
                for (Graph g : product) {
                    g.setNumFactors(4);
                    g.addFactors(factors);
                }
                z.addAll(product);
                z = isoFree.apply(z, null);
            }
        }
        return z;
    }
}
