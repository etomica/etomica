package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.FieldNodeCount;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.Node;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.property.IsBiconnected;
import etomica.graph.viewer.ClusterViewer;

public class VirialDiagrams {

    public static void main(String[] args) {
        final int n = 5;
        boolean multibody = false;
        boolean flex = false;
        
        char[] flexColors = new char[0];
        if (flex) {
           flexColors = new char[]{COLOR_CODE_0};
        }

        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList<Graph> topSet = new GraphList<Graph>(comp);

        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
        System.out.println("Xi");
        char eBond = 'A';//color of edge
        char fBond = 'f';
        char oneBond = 'o';
        char mBond = 'm';  // multi-body
        Set<Graph>[] multiXi = new Set[n+1];
        for (byte i=1; i<n+1; i++) {
            Graph g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
            g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
            eXi.add(g);
            
            multiXi[i] = new HashSet<Graph>();

            if (multibody && i>2) {
                g = GraphFactory.createGraph(i, BitmapFactory.createBitmap(i,true));
                g.coefficient().setDenominator((int)etomica.math.SpecialFunctions.factorial(i));
                for (Edge e : g.edges()) {
                    e.setColor(mBond);
                }
                eXi.add(g);
                multiXi[i].add(g);
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
        topSet.clear();
        topSet.addAll(fXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
//        ClusterViewer.createView("fXi", topSet);
        
        Set<Graph> lnfXi = new HashSet<Graph>();
        Set<Graph> fXipow = new HashSet<Graph>();
        MulFlexible mulFlex = new MulFlexible();
        MulFlexibleParameters mfp = new MulFlexibleParameters(flexColors);
        IsoFree isoFree = new IsoFree();
        fXipow.addAll(fXi);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();
        for (int i=1; i<n+1; i++) {

            lnfXi.addAll(fXipow);
            lnfXi = isoFree.apply(lnfXi, null);
            msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
            fXipow = isoFree.apply(mulScalar.apply(mulFlex.apply(fXipow, fXi, mfp), msp), null);
            FieldNodeCount truncater = new FieldNodeCount(new IteratorWrapper(fXipow.iterator()), n);
            Set<Graph> truncatedfXipow = new HashSet<Graph>();
            while (truncater.hasNext()) {
                truncatedfXipow.add(truncater.next());
            }
            fXipow = truncatedfXipow;

        }
        topSet.clear();
        topSet.addAll(lnfXi);
        System.out.println("\nlnfXi");
        for (Graph g : topSet) {
            System.out.println(g);
        }
//        ClusterViewer.createView("lnfXi", topSet);

        DifByNode opzdlnXidz = new DifByNode();
        DifParameters difParams = new DifParameters('A');
        Set<Graph> rho = opzdlnXidz.apply(lnfXi, difParams);
        
        GraphIterator iterator = new IsomorphismFilter(new IteratorWrapper(rho.iterator()));
        rho = new HashSet<Graph>();
        while (iterator.hasNext()) {
            rho.add(iterator.next());
        }
        System.out.println("\nrho");
        topSet.clear();
        topSet.addAll(rho);
        for (Graph g : topSet) {
            System.out.println(g);
        }
        
        HashSet<Graph>[] allRho = new HashSet[n+1];
        for (int i=0; i<n+1; i++) {
            allRho[i] = new HashSet<Graph>();
        }
        iterator = new IteratorWrapper(rho.iterator());
        while (iterator.hasNext()) {
            Graph g = iterator.next();
            allRho[g.nodeCount()].add(g);
        }
        
        Set<Graph> z = new HashSet<Graph>();
        
        // r = z + b*z^2 + c*z^3 + d*z^4 + e*z^5
        // z = r - b*r^2 + (2 b^2 - c) r^3 + (-5 b^3 + 5*b*c - d) r^4
        //     + (14 b^4 - 21 b^2*c + 3 c^2 + 6 b*d - e) r^5
        z.addAll(allRho[1]);
        if (n>1) {
            msp = new MulScalarParameters(new CoefficientImpl(-1,1));
            z.addAll(new MulScalar().apply(allRho[2], msp));
            
            if (n>2) {
                msp = new MulScalarParameters(new CoefficientImpl(2,1));
                Set<Graph> b2 = mulFlex.apply(allRho[2], allRho[2], mfp);
                z.addAll(new MulScalar().apply(b2, msp));
                msp = new MulScalarParameters(new CoefficientImpl(-1,1));
                z.addAll(new MulScalar().apply(allRho[3], msp));
                
                if (n>3) {
                    Set<Graph> b3 = mulFlex.apply(b2, allRho[2], mfp);
                    msp = new MulScalarParameters(new CoefficientImpl(-5,1));
                    z.addAll(new MulScalar().apply(b3, msp));
                    Set<Graph> bc = mulFlex.apply(allRho[2], allRho[3], mfp);
                    msp = new MulScalarParameters(new CoefficientImpl(5,1));
                    z.addAll(new MulScalar().apply(bc, msp));
                    msp = new MulScalarParameters(new CoefficientImpl(-1,1));
                    z.addAll(new MulScalar().apply(allRho[4], msp));
                    
                    if (n>4) {
                        Set<Graph> b4 = mulFlex.apply(b3, allRho[2], mfp);
                        msp = new MulScalarParameters(new CoefficientImpl(14,1));
                        z.addAll(new MulScalar().apply(b4, msp));
                        Set<Graph> b2c = mulFlex.apply(b2, allRho[3], mfp);
                        msp = new MulScalarParameters(new CoefficientImpl(-21,1));
                        z.addAll(new MulScalar().apply(b2c, msp));
                        Set<Graph> c2 = mulFlex.apply(allRho[3], allRho[3], mfp);
                        msp = new MulScalarParameters(new CoefficientImpl(3,1));
                        z.addAll(new MulScalar().apply(c2, msp));
                        Set<Graph> bd = mulFlex.apply(allRho[2], allRho[4], mfp);
                        msp = new MulScalarParameters(new CoefficientImpl(6,1));
                        z.addAll(new MulScalar().apply(bd, msp));
                        msp = new MulScalarParameters(new CoefficientImpl(-1,1));
                        z.addAll(new MulScalar().apply(allRho[5], msp));
                    }
                }
            }
        }

        System.out.println("\nz");
        z = isoFree.apply(z, null);
        topSet.clear();
        topSet.addAll(z);
        for (Graph g : topSet) {
            System.out.println(g);
        }
        
        HashSet<Graph> lnfXiOverV = new HashSet<Graph>();
        iterator = new IteratorWrapper(lnfXi.iterator());
        HashSet<Graph>[] alllnfXiOverV = new HashSet[n+1];
        for (int i=0; i<n+1; i++) {
            alllnfXiOverV[i] = new HashSet<Graph>();
        }
        while (iterator.hasNext()) {
            Graph g = iterator.next().copy();
//            g.getNode((byte)0).setType(TYPE_NODE_ROOT);
            lnfXiOverV.add(g);
            alllnfXiOverV[g.nodeCount()].add(g);
        }
        
        Set<Graph> p = new HashSet<Graph>();
        p.addAll(mulFlex.apply(alllnfXiOverV[1], z, mfp));
        p = isoFree.apply(p, null);


        if (n>1) {
            Set<Graph> z2 = new HashSet<Graph>();
            iterator = new FieldNodeCount(new IteratorWrapper(mulFlex.apply(z, z, mfp).iterator()), n-2);
            while (iterator.hasNext()) {
                z2.add(iterator.next());
            }

            p.addAll(mulFlex.apply(alllnfXiOverV[2], z2, mfp));

            p = isoFree.apply(p, null);

            if (n>2) {
                Set<Graph> z3 = new HashSet<Graph>();
                iterator = new FieldNodeCount(new IteratorWrapper(mulFlex.apply(z2, z, mfp).iterator()), n-3);
                while (iterator.hasNext()) {
                    z3.add(iterator.next());
                }
                p.addAll(mulFlex.apply(alllnfXiOverV[3], z3, mfp));
                p = isoFree.apply(p, null);
                
                if (n>3) {
                    Set<Graph> z4 = new HashSet<Graph>();
                    iterator = new FieldNodeCount(new IteratorWrapper(mulFlex.apply(z3, z, mfp).iterator()), n-4);
                    while (iterator.hasNext()) {
                        z4.add(iterator.next());
                    }
                    p.addAll(mulFlex.apply(alllnfXiOverV[4], z4, mfp));
                    p = isoFree.apply(p, null);

                    if (n>4) {
                        Set<Graph> z5 = new HashSet<Graph>();
                        iterator = new FieldNodeCount(new IteratorWrapper(mulFlex.apply(z4, z, mfp).iterator()), n-5);
                        while (iterator.hasNext()) {
                            z5.add(iterator.next());
                        }
                        p.addAll(mulFlex.apply(alllnfXiOverV[5], z5, mfp));
                        p = isoFree.apply(p, null);

                    }

                }

            }
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
    
    protected static void dump(Set<Graph> set, int max) {
        dump(set, max, false);
    }
    
    protected static void dump(Set<Graph> set, int max, boolean verbose) {
        IsBiconnected isBi = new IsBiconnected();
        for (int i=0; i<max+1; i++) {  // # of points
            System.out.println("** "+i+" **");
            for (int d=1; d>-1; d--) {  // isBiconnected
                if (verbose && i>2) System.out.println(d==1 ? "* biconnected *" : "* not biconnected *");
                for (int j=i*(i+1)/2; j>-1; j--) { // # of edges 
                    for (int k=0; k<2*i+2; k++) {  // # of root points
                        Iterator<Graph> iterator = set.iterator();
                        HashSet<Graph> printSet = new HashSet<Graph>();
                        while (iterator.hasNext()) {
                            Graph g = iterator.next();
                            if (g.edgeCount() != j || (isBi.check(g) != (d==1))) continue;
                            int fieldCount = 0;
                            for (Node node : g.nodes()) {
                                if (node.getType() == 'F') {
                                    fieldCount++;
                                }
                            }
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
