package etomica.virial.cluster;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.awt.Color;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphSuperSet;
import etomica.graph.model.GraphSuperSet.SubSetJudge;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.SetPropertyFilter;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.FieldNodeCountMax;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;

public class VirialDiagramsMix {

    public static void main(String[] args) {
        final int n = 4;
        final char nodeA = Metadata.COLOR_CODE_0;
        final char nodeB = Metadata.COLOR_CODE_1;
        
        SubSetJudge nFieldJudge = new SubSetJudge() {
            public int getSetId(Graph g) {
                int fieldCount = 0;
                for (Node node : g.nodes()) {
                    if (node.getType() == 'F') {
                        fieldCount++;
                    }
                }
                return fieldCount;
            }
        };
        SubSetJudge biConJudge = new SubSetJudge() {
            public int getSetId(Graph g) { return isBiCon.check(g) ? 0 : 1;}
            protected final IsBiconnected isBiCon = new IsBiconnected();
        };
        SubSetJudge nAPointsJudge = new SubSetJudge() {
            public int getSetId(Graph g) {
                int[] factors = g.factors();
                if (factors.length > 0) {
                    return factors[1];
                }
                int nodeACount = 0;
                for (Node node : g.nodes()) {
                    if (node.getColor() == nodeA) {
                        nodeACount++;
                    }
                }
                return g.nodeCount() - nodeACount;
            }
        };
        SubSetJudge nEdgesJudge = new SubSetJudge() {
            public int getSetId(Graph g) {
                return n*(n-1)/2 - g.edgeCount();
            }
        };
        SubSetJudge nRootJudge = new SubSetJudge() {
            public int getSetId(Graph g) {
                int fieldCount = 0;
                for (Node node : g.nodes()) {
                    if (node.getType() == 'F') {
                        fieldCount++;
                    }
                }
                return g.nodeCount() - fieldCount;
            }
        };

        // array of sets containing different # of field points
        GraphSuperSet<Graph>[] fieldPointsSets = new GraphSuperSet[n+1];
        for (int l=0; l<n+1; l++) {
            // array of sets containing biconnected and not-biconnected sets
            GraphSuperSet<Graph>[] biConSets = new GraphSuperSet[2];
            for (int k=0; k<2; k++) {
                // array of sets containing different # of A points
                GraphSuperSet<Graph>[] aPointsSets = new GraphSuperSet[n+1];
                for (int i=0; i<n+1; i++) {
                    // array of sets containing different # of edges
                    GraphSuperSet<Graph>[] nEdgesSets = new GraphSuperSet[n*(n+1)/2+1];
                    for (int j=0; j<n*(n+1)/2+1; j++) {
                        // array of sets containing different # of root points
                        Set<Graph>[] rootSets = new Set[n+1];  // low-level sets
                        for (int m=0; m<n+1; m++) {
                            rootSets[m] = new HashSet<Graph>();
                        }
                        nEdgesSets[j] = new GraphSuperSet<Graph>(nRootJudge, rootSets);
                    }
                    aPointsSets[i] = new GraphSuperSet<Graph>(nEdgesJudge, nEdgesSets);
                }
                biConSets[k] = new GraphSuperSet<Graph>(nAPointsJudge, aPointsSets);
            }
            fieldPointsSets[l] = new GraphSuperSet<Graph>(biConJudge, biConSets);
        }
        GraphSuperSet<Graph> topSet = new GraphSuperSet<Graph>(nFieldJudge, fieldPointsSets);
        
        
        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
        System.out.println("Xi");
        char eBondAA = Metadata.COLOR_CODE_0;
        char eBondAB = Metadata.COLOR_CODE_1;
        char eBondBB = Metadata.COLOR_CODE_2;
        for (byte i=1; i<n+1; i++) {
            for (byte j=0; j<i+1; j++) {
                // j points of color A, i-j points of color B
                Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                g.coefficient().setDenominator((int)(SpecialFunctions.factorial(i-j)*SpecialFunctions.factorial(j)));
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
        
        GraphIterator iterator = new PropertyFilter(new IteratorWrapper(fXi.iterator()), new IsConnected());
        Set<Graph> lnfXi = new HashSet<Graph>();
        while (iterator.hasNext()) {
            Graph g = iterator.next();
            lnfXi.add(g);
        }

        DifByNode opzdlnXidzA = new DifByNode();
        DifParameters difParams = new DifParameters('A');
        Set<Graph> rhoA = opzdlnXidzA.apply(lnfXi, difParams);
        
        iterator = new IsomorphismFilter(new IteratorWrapper(rhoA.iterator()));
        rhoA = new HashSet<Graph>();
        while (iterator.hasNext()) {
            rhoA.add(iterator.next());
        }
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
//        ClusterViewer.createView("rhoB", rhoB);
//        System.out.println("\nrhoB");
//        dump(rhoB, n-1, nodeA);

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
        MulFlexible mulFlex = new MulFlexible();
        
        Set<Graph> zA = invertSeries(new Set[][][]{allRhoA, allRhoB}, n, "A", "B");
        Set<Graph> zB = invertSeries(new Set[][][]{allRhoB, allRhoA}, n, "B", "A");
        
        SetPropertyFilter truncater = new SetPropertyFilter(new FieldNodeCountMax(n-1));
        zA = isoFree.apply(truncater.apply(zA, null), null);
        zB = isoFree.apply(truncater.apply(zB, null), null);
        topSet.clear();
        topSet.addAll(zA);
        System.out.println("\nzA");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("zA", topSet);

        if (false) {
            zA = invertSeries(allRhoA, allRhoB, n);
            zA = isoFree.apply(truncater.apply(zA, null), null);
            System.out.println("\nold zA");
            dump(zA, n-1, nodeA);
            zB = invertSeries(allRhoB, allRhoA, n);
            zB = isoFree.apply(truncater.apply(zB, null), null);
            System.out.println("\nold zB");
            dump(zB, n-1, nodeA);
            System.exit(1);
        }


        // now we have a series expansion for rhoA in powers of zB, with rhoA in the coefficients

        // now we have a series expansion for rhoB in powers of zB, with rhoA in the coefficients
        // we can invert that series and get zB in terms of rhoA and rhoB
        
        HashSet<Graph> lnfXiOverV = new HashSet<Graph>();
        HashSet<Graph>[][] alllnfXiOverV = new HashSet[n+1][n+1];
        for (int i=0; i<n+1; i++) {
            for (int j=0; j<n+1; j++) {
                alllnfXiOverV[i][j] = new HashSet<Graph>();
            }
        }
        for (Graph g : lnfXi) {
            g = g.copy();
            g.setNumFactors(2);
            g.getNode((byte)0).setType(TYPE_NODE_ROOT);
            lnfXiOverV.add(g);
            int nodeACount = 0;
            for (Node node : g.nodes()) {
                if (node.getColor() == nodeA) {
                    nodeACount++;
                }
            }
            alllnfXiOverV[nodeACount][g.nodeCount()-nodeACount].add(g);
        }
        
        Set<Graph> p = new HashSet<Graph>();

        Set<Graph> zApow = new HashSet<Graph>();
        zApow.addAll(zA);
        Set<Graph> zBpow = new HashSet<Graph>();
        zBpow.addAll(zB);
        for (int j=1; j<n+1; j++) {
            p.addAll(truncater.apply(mulFlex.apply(alllnfXiOverV[0][j], zBpow, null), null));
            p = isoFree.apply(p, null);
            zBpow = isoFree.apply(truncater.apply(mulFlex.apply(zBpow, zB, null), null), null);
        }
        for (int i=1; i<n+1; i++) {
            Set<Graph> zABpow = new HashSet<Graph>();
            zABpow.addAll(zApow);
            for (int j=0; j<n+1-i; j++) {
                p.addAll(truncater.apply(mulFlex.apply(alllnfXiOverV[i][j], zABpow, null), null));
                p = isoFree.apply(p, null);
                zABpow = isoFree.apply(truncater.apply(mulFlex.apply(zABpow, zB, null), null), null);
            }
            zApow = isoFree.apply(truncater.apply(mulFlex.apply(zApow, zA, null), null), null);
        }
        
        topSet.clear();
        topSet.addAll(p);
        System.out.println("\nP");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("P", topSet);

//        iterator = new FieldNodeCount(new IteratorWrapper(p.iterator()), n-1);
//        System.out.println("\nP (in terms of rho)");
//        while (iterator.hasNext()) {
//            System.out.println(iterator.next());
//        }

    }
    
    public static class Term {
        public int coefficient;
        public Set<Graph>[] factors;
        public Term(int coefficient, Set<Graph>[] factors) {
            this.coefficient = coefficient;
            this.factors = factors;
        }
    }

    protected static Set<Graph> invertSeries(Set<Graph>[][][] allRho, int n, String comp0, String comp1) {
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
                factors = new int[]{nT-n1,n1};
            }
            else {
                factors = new int[]{n1,nT-n1};
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
                        product = mulFlex.apply(product, allRho[term[j][0]][term[j][1]][term[j][2]], null);
                        product = isoFree.apply(truncater.apply(product, null), null);
                    }
                }
                product = mulScalar.apply(product, new MulScalarParameters(coefficient,1));
                for (Graph g : product) {
                    g.setNumFactors(2);
                    g.addFactors(factors);
                }
                z.addAll(product);
                z = isoFree.apply(z, null);
            }
        }
        return z;
    }

    
    protected static Set<Graph> invertSeries(Set<Graph>[][] allRho1, Set<Graph>[][] allRho2, int n) {
        MulScalar mulScalar = new MulScalar();
        MulScalarParameters minus = new MulScalarParameters(-1, 1);
        MulFlexible mulFlex = new MulFlexible();
        Set<Graph> z1 = new HashSet<Graph>();
//        ArrayList<Term> terms = new ArrayList<Term>();
//        
//        terms.add(new Term(1, new Set[]{allRho1[1][0]}));
//        int[] coefs = new int[]{1, -1, -1, -1, 2, -1, 3};
//        int[][][] terms = new int[][][]{{{1,1,0}}, {{1,2,0}}, {{1,1,1}}, {{1,3,0}}, {{1,2,0},{2,0}}, {{2,1}}, {{1,1},{2,0}}, };
        z1.addAll(allRho1[1][0]);

        if (n>1) {
            z1.addAll(mulScalar.apply(allRho1[2][0], minus));
            z1.addAll(mulScalar.apply(allRho1[1][1], minus));

            if (n>2) {
                // msp = -1
                z1.addAll(mulScalar.apply(allRho1[3][0], minus));
                z1.addAll(mulScalar.apply(allRho1[2][1], minus));
                z1.addAll(mulScalar.apply(allRho1[1][2], minus));

                // from C30
                Set<Graph> C20_2 = mulFlex.apply(allRho1[2][0], allRho1[2][0], null);
                z1.addAll(mulScalar.apply(C20_2, new MulScalarParameters(2,1)));

                // from C21
                Set<Graph> C20C11 = mulFlex.apply(allRho1[2][0], allRho1[1][1], null);
                z1.addAll(mulScalar.apply(C20C11, new MulScalarParameters(3,1)));
                Set<Graph> C11D11 = mulFlex.apply(allRho1[1][1], allRho2[1][1], null);
                z1.addAll(C11D11);

                // from C12
                Set<Graph> C11_2 = mulFlex.apply(allRho1[1][1], allRho1[1][1], null);
                z1.addAll(C11_2);
                Set<Graph> C11D20 = mulFlex.apply(allRho1[1][1], allRho2[2][0], null);
                z1.addAll(C11D20);
                
//                if (n>3) {
//                    Set<Graph> b3 = mulFlex.apply(b2, allRho1[2], null);
//                    msp = new MulScalarParameters(new CoefficientImpl(-5,1));
//                    zA.addAll(new MulScalar().apply(b3, msp));
//                    Set<Graph> bc = mulFlex.apply(allRho1[2], allRho1[3], null);
//                    msp = new MulScalarParameters(new CoefficientImpl(5,1));
//                    zA.addAll(new MulScalar().apply(bc, msp));
//                    msp = new MulScalarParameters(new CoefficientImpl(-1,1));
//                    zA.addAll(new MulScalar().apply(allRho1[4], msp));
//                    
//                    if (n>4) {
//                        Set<Graph> b4 = mulFlex.apply(b3, allRho1[2], null);
//                        msp = new MulScalarParameters(new CoefficientImpl(14,1));
//                        zA.addAll(new MulScalar().apply(b4, msp));
//                        Set<Graph> b2c = mulFlex.apply(b2, allRho1[3], null);
//                        msp = new MulScalarParameters(new CoefficientImpl(-21,1));
//                        zA.addAll(new MulScalar().apply(b2c, msp));
//                        Set<Graph> c2 = mulFlex.apply(allRho1[3], allRho1[3], null);
//                        msp = new MulScalarParameters(new CoefficientImpl(3,1));
//                        zA.addAll(new MulScalar().apply(c2, msp));
//                        Set<Graph> bd = mulFlex.apply(allRho1[2], allRho1[4], null);
//                        msp = new MulScalarParameters(new CoefficientImpl(6,1));
//                        zA.addAll(new MulScalar().apply(bd, msp));
//                        msp = new MulScalarParameters(new CoefficientImpl(-1,1));
//                        zA.addAll(new MulScalar().apply(allRho1[5], msp));
//                    }
//                }
            }
        }
        
        return z1;
    }
    
    protected static void dump(Set<Graph> set, int max, char colorA) {
        dump(set, max, false, colorA);
    }
    
    protected static void dump(Set<Graph> set, int max, boolean verbose, char colorA) {
        IsBiconnected isBi = new IsBiconnected();
        for (int i=0; i<max+1; i++) {  // # of points
            System.out.println("** "+i+" **");
            for (int d=1; d>-1; d--) {  // isBiconnected
                if (verbose && i>2) System.out.println(d==1 ? "* biconnected *" : "* not biconnected *");
                for (int j=i*(i+1)/2; j>-1; j--) { // # of edges
                    for (int l=2*i+2; l>-1; l--) {
                        for (int k=0; k<2*i+2; k++) {  // # of root points
                            for (Graph g : set) {
                                if (g.edgeCount() != j || (isBi.check(g) != (d==1))) continue;
                                int[] factors = g.factors();
                                if (factors.length > 0) {
                                    if (factors[0] != l) continue;
                                }
                                int fieldCount = 0;
                                int nodeACount = 0;
                                for (Node node : g.nodes()) {
                                    if (node.getType() == 'F') {
                                        fieldCount++;
                                    }
                                    if (node.getColor() == colorA) {
                                        nodeACount++;
                                    }
                                }
                                if (factors.length == 0 && nodeACount != l) continue;
                                if (fieldCount == i && (g.nodeCount()-fieldCount == k)) {
                                    System.out.println(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
