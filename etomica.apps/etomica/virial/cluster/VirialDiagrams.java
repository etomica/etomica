package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
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
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitOne;
import etomica.graph.operations.SplitOne.SplitOneParameters;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
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
        MulFlexibleParameters mfp = new MulFlexibleParameters(flexColors, (byte)n);
        IsoFree isoFree = new IsoFree();
        fXipow.addAll(fXi);
        MulScalarParameters msp = null;
        MulScalar mulScalar = new MulScalar();
        for (int i=1; i<n+1; i++) {

            lnfXi.addAll(fXipow);
            lnfXi = isoFree.apply(lnfXi, null);
            msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
            fXipow = isoFree.apply(mulScalar.apply(mulFlex.apply(fXipow, fXi, mfp), msp), null);
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
        Set<Graph> rho = isoFree.apply(opzdlnXidz.apply(lnfXi, difParams), null);
        
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

        z.addAll(allRho[1]);
        for (int i=2; i<n+1; i++) {
            Set<Graph>[] zPow = new Set[n+1];
            zPow[1] = new HashSet<Graph>();
            zPow[1].addAll(z);
            for (int j=2; j<i+1; j++) {
                zPow[j] = new HashSet<Graph>();
                zPow[j] = isoFree.apply(mulFlex.apply(zPow[j-1], z, mfp), null);
            }
            z = new HashSet<Graph>();
            z.addAll(allRho[1]);
            msp = new MulScalarParameters(new CoefficientImpl(-1,1));
            for (int j=2; j<i+1; j++) {
                z.addAll(mulScalar.apply(mulFlex.apply(allRho[j], zPow[j], mfp), msp));
            }
        }

        System.out.println("\nz");
        z = isoFree.apply(z, null);
        topSet.clear();
        topSet.addAll(z);
        for (Graph g : topSet) {
            System.out.println(g);
        }
        
        Decorate decorate = new Decorate();
        DecorateParameters dp = new DecorateParameters(COLOR_CODE_0, mfp);
        
        Set<Graph> p = decorate.apply(lnfXi, z, dp);
        p = isoFree.apply(p, null);
        
        Set<Graph> newP = new HashSet<Graph>();

        // attempt to factor any graphs with an articulation point
        HashMap<Graph,Graph> cancelSet = new HashMap<Graph,Graph>();
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
            p = isoFree.apply(newP, null);

            // perform Ree-Hoover substitution (brute-force)
            IsBiconnected isBi = new IsBiconnected();
            char nfBond = 'F';
            SplitOneParameters splitOneParameters = new SplitOneParameters(eBond, nfBond);
            SplitOne splitOne = new SplitOne();
            msp = new MulScalarParameters(-1, 1);
            newP.clear();
            for (Graph g : p) {
                if (isBi.check(g)) {
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
                else {
                    newP.add(g);
                }
            }
            p = isoFree.apply(newP, null);
        }
        else {
            MulFlexibleParameters mfp2 = new MulFlexibleParameters(new char[0], (byte)n);
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            FactorOnce factor = new FactorOnce();
            newP.clear();
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                newP.add(g.copy());
                if (con && ap) {
//                    System.out.println("\nfactoring\n"+g);
                    Graph gf = factor.apply(g, mfp2);
                    newP.add(gf);
                }
            }
            p = isoFree.apply(newP, null);
            
            for (Graph g : p) {
                boolean ap = hap.check(g);
                boolean con = hap.isConnected();
                if (con && ap) {
                    Graph gf = factor.apply(g, mfp2);
                    cancelSet.put(g,gf);
                }
            }
        }

        topSet.clear();
        topSet.addAll(p);
        System.out.println("\nP");
        for (Graph g : topSet) {
            System.out.println(g);
            Graph cancelGraph = cancelSet.get(g);
            if (cancelGraph != null) {
                System.out.println(" -  "+cancelGraph);
            }
        }
        ClusterViewer.createView("P", topSet);
    }
}
