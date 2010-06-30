package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.IsomorphismFilter;
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
import etomica.graph.operations.PCopy;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitOne;
import etomica.graph.operations.SplitOne.SplitOneParameters;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.viewer.ClusterViewer;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterSumShell;
import etomica.virial.MayerFunction;

public class VirialDiagrams {

    protected final int n;
    protected final boolean flex;
    protected final boolean multibody;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Map<Graph,Graph> cancelMap;
    
    public static void main(String[] args) {
        final int n = 4;
        boolean multibody = false;
        boolean flex = true;
        new VirialDiagrams(n, multibody, flex, true);
    }
    
    public VirialDiagrams(int n, boolean multibody, boolean flex) {
        this(n, multibody, flex, false);
    }
    
    public VirialDiagrams(int n, boolean multibody, boolean flex, boolean interactive) {
        this.multibody = multibody;
        this.flex = flex;
        this.n = n;
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        makeVirialDiagrams(interactive);
    }
    
    public Set<Graph> getVirialGraphs() {
        return p;
    }

    public ClusterSum makeVirialCluster(MayerFunction f, MayerFunction e) {
        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        GraphList<Graph> allP = new GraphList<Graph>();
        allP.addAll(p);
        allP.addAll(cancelP);
        for (Graph g : allP) {
            int fieldCount = 0;
            for (Node node : g.nodes()) {
              if (node.getType() == 'F') {
                fieldCount++;
              }
            }
            if (fieldCount == n) {
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
    
    public void populateEFBonds(Graph g, ArrayList<ClusterBonds> allBonds, boolean swap) {
        ArrayList<int[]> fbonds = new ArrayList<int[]>();
        ArrayList<int[]> ebonds = new ArrayList<int[]>();
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
                    if (edgeColor == 'f') {
                        fbonds.add(new int[]{n1,n2});
                    }
                    else if (edgeColor == 'A') {
                        ebonds.add(new int[]{n1,n2});
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

    public ClusterSumShell[] makeSingleVirialClusters(ClusterSum coreCluster, MayerFunction e, MayerFunction f) {
        ArrayList<ClusterSumShell> allClusters = new ArrayList<ClusterSumShell>();
        double[] w = new double[]{1};
        if (flex) {
            w = new double[]{0.5,0.5};
        }
        for (Graph g : p) {
            int fieldCount = 0;
            for (Node node : g.nodes()) {
              if (node.getType() == 'F') {
                fieldCount++;
              }
            }
            if (fieldCount == n) {
                System.out.println(allClusters.size()+" "+g);
                double[] thisW = w;
                ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
                populateEFBonds(g, allBonds, false);
                populateEFBonds(g, allBonds, true);
                if (flex && cancelMap.get(g) != null) {
                    Graph cg = cancelMap.get(g);
                    populateEFBonds(cg, allBonds, false);
                    populateEFBonds(cg, allBonds, true);
                    thisW = new double[]{0.5,0.5,-0.5,-0.5};
                }
                if (n > 3 && !flex && !multibody) {
                    allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{e,f}));
                }
                else if (!multibody) {
                    allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, new MayerFunction[]{f}));
                }
            }
        }
        return allClusters.toArray(new ClusterSumShell[0]);
    }
    
    public Set<Graph> getExtraDisconnectedVirialGraphs() {
        return disconnectedP;
    }
    
    public void makeVirialDiagrams(boolean interactive) {
        
        char[] flexColors = new char[0];
        if (flex) {
           flexColors = new char[]{COLOR_CODE_0};
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

        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        GraphList<Graph> topSet = new GraphList<Graph>(comp);

        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
        char eBond = 'A';//color of edge
        char fBond = 'f';
        char oneBond = 'o';
        char mBond = 'm';  // multi-body
        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(mBond, 1);
        colorOrderMap.put(eBond, 2);
        colorOrderMap.put(fBond, 3);
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

        if (interactive) {
            System.out.println("Xi");
            topSet.addAll(eXi);
            for (Graph g : topSet) {
                System.out.println(g);
            }
//            ClusterViewer.createView("eXi", topSet);
        }

        
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
        Set<Graph> setOfSubstituted = split.apply(eXi, bonds);

        DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
        DeleteEdge deleteEdge = new DeleteEdge();
        //set of full star diagrams with f bonds
        Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
        if (interactive) {
            System.out.println("\nXi with f bonds");
            topSet.clear();
            topSet.addAll(fXi);
            for (Graph g : topSet) {
                System.out.println(g);
            }
//            ClusterViewer.createView("fXi", topSet);
        }
        
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
        if (interactive) {
            topSet.clear();
            topSet.addAll(lnfXi);
            System.out.println("\nlnfXi");
            for (Graph g : topSet) {
                System.out.println(g);
            }
//            ClusterViewer.createView("lnfXi", topSet);
        }

        DifByNode opzdlnXidz = new DifByNode();
        DifParameters difParams = new DifParameters('A');
        Set<Graph> rho = isoFree.apply(opzdlnXidz.apply(lnfXi, difParams), null);

        if (interactive) {
            System.out.println("\nrho");
            topSet.clear();
            topSet.addAll(rho);
            for (Graph g : topSet) {
                System.out.println(g);
            }
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

        MulFlexibleParameters mfpnm1 = new MulFlexibleParameters(flexColors, (byte)(n-1));
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
        if (interactive) {
            System.out.println("\nz");
            topSet.clear();
            topSet.addAll(z);
            for (Graph g : topSet) {
                System.out.println(g);
            }
        }
        
        Decorate decorate = new Decorate();
        DecorateParameters dp = new DecorateParameters(COLOR_CODE_0, mfp);
        
        p = decorate.apply(lnfXi, z, dp);
        p = isoFree.apply(p, null);
        
        Set<Graph> newP = new HashSet<Graph>();

        // attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Graph>();
        cancelP = new GraphList<Graph>();
        disconnectedP = new HashSet<Graph>();
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
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Relabel relabel = new Relabel();
            FactorOnce factor = new FactorOnce();
            FactorOnceParameters fop = new FactorOnceParameters((byte)0, new char[0]);
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
                    Graph gf = factor.apply(g, fop);
                    disconnectedP.add(gf);
                    gf = mulScalar.apply(gf, msp);
                    cancelP.add(gf);
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
            // we don't need to re-isofree p, we know that's still good.
            // some of our new disconnected diagrams might condense with the old ones
            disconnectedP = isoFree.apply(disconnectedP, null);

            // we want to condense cancelP (in case multiple diagrams were factored into the
            // same one), but want to be careful not to permute bonds.
            PCopy pcopy = new PCopy();
            IteratorWrapper wrapper = new IteratorWrapper(pcopy.apply(cancelP, null).iterator());
            GraphIterator isomorphs = new IsomorphismFilter(wrapper);
            cancelP = new GraphList<Graph>();
            while (isomorphs.hasNext()) {
                cancelP.add(isomorphs.next());
            }
        }
        p = newP;
        
        if (interactive) {
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

        GraphList<Graph> pFinal = new GraphList<Graph>();
        pFinal.addAll(p);
        p = pFinal;
        GraphList<Graph> disconnectedPFinal = new GraphList<Graph>();
        disconnectedPFinal.addAll(disconnectedP);
        disconnectedP = disconnectedPFinal;

    }
}
