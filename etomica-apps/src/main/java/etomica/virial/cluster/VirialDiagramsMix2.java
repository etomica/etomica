/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Edge;
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
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.Decorate;
import etomica.graph.operations.Decorate.DecorateParameters;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.DifByNode;
import etomica.graph.operations.DifParameters;
import etomica.graph.operations.Factor;
import etomica.graph.operations.FactorOnce;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.FactorOnce.FactorOnceParameters;
import etomica.graph.operations.GraphOpMaxRoot;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.Property;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumShell;
import etomica.virial.MayerFunction;

public class VirialDiagramsMix2 {

    protected final int n;
    protected char[] nodeColors;
    protected final boolean[] flex;
    protected char[] flexColors;
    protected final boolean isInteractive;
    protected boolean doReeHoover;
    protected Set<Graph> p, cancelP, disconnectedP;
    protected Set<Graph> rhoA, rhoB;
    protected Set<Graph> lnfXi;
    protected Map<Graph,Graph> cancelMap;
    protected boolean doShortcut;
    protected boolean doKeepEBonds;
    protected boolean doDisconnectedMatching = true;
    public char fBond, eBond, excBond, efbcBond;
    protected boolean allPermutations = false;
    
    
    public static void main(String[] args) {
        int n = 3;
        int[] num = new int[]{2,1};
        boolean[] flex = new boolean[]{false, false};// 1st species:black; 2nd species:red
        boolean doKeepEBonds = false;
        boolean doReeHoover = false;
        VirialDiagramsMix2 virialDiagrams = new VirialDiagramsMix2(n, flex, true);
        virialDiagrams.setDoReeHoover(doReeHoover);
        virialDiagrams.setDoKeepEBonds(doKeepEBonds);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.makeVirialDiagrams();
        virialDiagrams.getMSMCGraphs(true, -1, num);
    }


	public VirialDiagramsMix2(int n , boolean[] flex) {
        this( n, flex, false);
    }

    public VirialDiagramsMix2(int n, boolean[] flex, boolean interactive) {
        this.flex = flex;
        this.n = n;
        this.isInteractive = interactive;
        doReeHoover = false;
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

    public void setDoDisconnectedMatching(boolean newDoDisconnectedMatching) {
        doDisconnectedMatching = newDoDisconnectedMatching;
    }
    
    public Set<Graph> getExtraDisconnectedVirialGraphs(int[] numPoints) {
        if (p == null) {
            makeVirialDiagrams();
        }
        GraphList dpn = makeGraphList();
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        Property mic = new MaxIsomorphCriteriaMixture();
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), mic);
        for (Graph g : disconnectedP) {
        	boolean correctComposition = true;
        	for ( int i = 0;correctComposition&& i<numPoints.length ;i++){
        		correctComposition = numPoints[i]==g.factors()[nodeColors.length+i];
        	}
            if (correctComposition) {
                g = maxIsomorph.apply(g, mip);
	
                dpn.add(g);
            }
        }
       
        return dpn;
    }

    public Set<Graph> getSplitDisconnectedVirialGraphs(Graph g) {
        SplitGraph splitGraph = new SplitGraph();
        MaxIsomorph maxIsomorph = new MaxIsomorph();
//        Property happyArticulation = new ArticulatedAt0(false, '0');
//        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), happyArticulation);
        // apply MaxIsomorphCriteriaMixture so that in simulation, the diagram label apprears in a desired way!
        Property mic = new MaxIsomorphCriteriaMixture();
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), mic);
        // we want gSplit unsorted
        Set<Graph> gSplit = new GraphList(null);
        Set<Graph> gSplit1 = splitGraph.apply(g);
        for (Graph gs : gSplit1) {
            // the graph we get from splitting might not be in our preferred bonding arrangement
            Graph gsmax = maxIsomorph.apply(gs, mip);
            gSplit.add(gsmax);
        }
        return gSplit;
    }

    public GraphList makeGraphList() {
    	   	
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNodeColors(nodeColors));
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        return new GraphList(comp);
    }
    
    public static boolean graphHasEdgeColor(Graph g, char color) {
        for (Edge edge : g.edges()) {
            if (edge.getColor() == color) {
                return true;
            }
        }
        return false;
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
        efbcBond = 'b';

        colorOrderMap.put(oneBond, 0);
        colorOrderMap.put(efbcBond, 2);
        colorOrderMap.put(fBond, 3);
        colorOrderMap.put(eBond, 4);

        Metadata.COLOR_MAP.put(eBond, "red");
        Metadata.COLOR_MAP.put(fBond, "green");
        Metadata.COLOR_MAP.put(efbcBond, "fuchsia");
        Metadata.COLOR_MAP.put(excBond, "red");
        Metadata.DASH_MAP.put(excBond, 3);
        
        Set<Graph> topSet = makeGraphList();
        // ==================================================== eXi ======================================================= //
        Set<Graph> eXi = new HashSet<Graph>();//set of full star diagrams with e bonds
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
            }
        }
        if (isInteractive) {
            System.out.println("eXi");
        	topSet.addAll(eXi);
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        	ClusterViewer.createView("eXi", topSet);
        }
        // ==================================================== fXi ======================================================= //
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
        Set<Graph> setOfSubstituted = split.apply(eXi, bonds);
        DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
        DeleteEdge deleteEdge = new DeleteEdge();
        
        Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
        IsoFree isoFree = new IsoFree();
        fXi = isoFree.apply(fXi, null);
        if (isInteractive) {
        	System.out.println("\nXi with f bonds");//set of full star diagrams with f bonds
        	ClusterViewer.createView("fXi", fXi);
        	topSet.clear();
        	topSet.addAll(fXi);
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        }
        // ==================================================== lnfXi ======================================================= //
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
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(lnfXi);
        	System.out.println("\nlnfXi");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        	ClusterViewer.createView("lnfXi", lnfXi);
        }
        // put lnfXi ordered
        if (isInteractive) {
        	System.out.println("\nlnfXi ordered");
            topSet.clear();
            topSet.addAll(lnfXi);
            for (Graph g : topSet) {
            	System.out.println(g);
            }
            ClusterViewer.createView("lnfXi, ordered", topSet);
        }
        
                
        // ----------------- rhoA and rhoB  ------------------ //
        DifByNode opzdlnXidzA = new DifByNode();
        DifParameters difParams = new DifParameters('A');
        rhoA = isoFree.apply(opzdlnXidzA.apply(lnfXi, difParams), null);
        if(isInteractive){
        	topSet.clear();
        	topSet.addAll(rhoA);
        	ClusterViewer.createView("rhoA", topSet);
        	System.out.println("\nrhoA");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        }
        DifByNode opzdlnXidzB = new DifByNode();
        difParams = new DifParameters('B');
        rhoB = isoFree.apply(opzdlnXidzB.apply(lnfXi, difParams), null);
        if(isInteractive){
        	topSet.clear();
        	topSet.addAll(rhoB);
        	ClusterViewer.createView("rhoB", topSet);
        	System.out.println("\nrhoB");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        }
    } //end makeRhoDiagrams method
    
    public void makeVirialDiagrams() {
        if (rhoA == null) {
            makeRhoDiagrams();
        }
        // ======================= inverse series, zA and zB in power series of rhoA and rhoB=================== //
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
                g.factors()[0] = 0;/////// ??????????????????????????
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
        if(isInteractive){
        	topSet.clear();
        	topSet.addAll(zA);
        	System.out.println("\nzA");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        	ClusterViewer.createView("zA", topSet);
        }
        if(isInteractive){
        	topSet.clear();
        	topSet.addAll(zB);
        	System.out.println("\nzB");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        	ClusterViewer.createView("zB", topSet);
        }
        // refer to VirialDiagrams2 Line 1452
        MulFlexibleParameters mfpn = MulFlexibleParameters.makeParameters(nodeColors, (byte)n);///????????????????????
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        Property mic = new MaxIsomorphCriteriaMixture();
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), mic);
        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        
        //----------------------- decorate lnfXi by using zA and zB  to get p -----------------------------//
        p = decorate.apply(lnfXi, zA, new DecorateParameters(0, mfpn));
        p = decorate.apply(p, zB, new DecorateParameters(1, mfpn));
        p = maxIsomorph.apply(isoFree.apply(p, null), mip);

   //     p = isoFree.apply(maxIsomorph.apply(p, mip), null);// based on maxisomorphCriteriaMixture
        if(isInteractive){
        	topSet.clear();
        	topSet.addAll(p);
        	System.out.println("\nfirst P");
        	for (Graph g : topSet) {
        		System.out.println(g);
        	}
        	ClusterViewer.createView("first P", topSet);
        }
        Set<Graph> newP = new HashSet<Graph>();// refer to line 1455
        //----------------------- factor diagrams with rigid articulation points-----------------------------------------//
        if (flexColors.length < nodeColors.length) {//at least 1 species is rigid, refer to line 1459, !flex
        	Factor factor = new Factor();
        	MulFlexibleParameters factorParameters = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        	for (Graph g : p) {
        		boolean ap = hap.check(g);
        		boolean con = hap.isConnected();
        		if (ap) {//?????????????????
        //		if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {
        			//factor the diagrams if the a.p. is rigid
        			boolean factorable = false;// can be decomposed to 2 or more diagrams
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
        				// graph has an articulation point, but it was flexible
        				newP.add(g.copy());
        			}
        		}
        		else {
        			newP.add(g.copy());
        		}
        	}
            if (isInteractive) {
            	topSet.clear();
            	topSet.addAll(newP);
            	System.out.println("\nnewP");
            	ClusterViewer.createView("newP after factoring away the rigid ap diagrams", topSet);
            }
        	p = isoFree.apply(newP, null);// combine all the isomorphs in newP to get P
        	newP.clear();
        	newP.addAll(maxIsomorph.apply(p, mip));
        	p = newP;
        	newP = makeGraphList();
        	
        	// or do this(2 ways after newP):
        	// newP = isoFree.apply(maxIsomorph.apply(newP,mip),null);
        	// (1)p = newP;
        	// (1)newP = makeGraphList();
        	// (2)p.clear();
        	// (2)p.addAll(newP); then p and newP wont link together
            if(isInteractive){
            	topSet.clear();
            	topSet.addAll(p);
            	System.out.println("\nP");
            	ClusterViewer.createView("p after factoring away the rigid ap diagrams", topSet);
            }    	
        }
        
        //----------------------- factor diagrams with flex articulation points-----------------------------------------//
        disconnectedP = new HashSet<Graph>();// attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Graph>();
        if (flexColors.length > 0) {  // refer to VirialDiagram2 Line 1489
            // pretend everything is fully flexible
            FactorOnce factor = new FactorOnce();
            FactorOnceParameters fop = null;
            // match up singly-connected (in p) with disconnected diagrams.
            // we have to do this last so that our cancelMap remains valid.
            newP.clear();///??????????????????
            msp = new MulScalarParameters(-1, 1);

            for (Graph g : p) {
            	boolean ap = hap.check(g);
            	boolean con = hap.isConnected();
            	if (con && ap ) {// search singly-connected graph
            		fop = null;
            		//loop over all articulation points in g to find where to factor the diagram
            		for( Node node : g.nodes() ) {
            			byte nodeID = node.getId();
            			boolean isArticulationPt = hap.getArticulationPoints().contains(nodeID);//check whether the node is an a.p.
            			if (isArticulationPt){
            				fop = new FactorOnceParameters(nodeID, false);
            				break;// find 1 a.p. and we are good!
            			}
        				
            		}
            		// newP will contain connected diagrams
            		g = g.copy();//??????????????????????refer to ViralDiagram2 Line 1509
            		newP.add(g);// newP now has singly-connected diagrams ONLY
            		Set<Graph> gfSet = factor.apply(g, fop);// factored diagrams set of g
            		Graph gf = gfSet.iterator().next();// only 1 iterate
            		disconnectedP.add(gf);
            		gf = mulScalar.apply(gf, msp);
            		cancelMap.put(g, gf);
            	}
            	else if (con) {// this is a biconnected diagram;
            		newP.add(g.copy());
            	}
            	else {// this is a disconnected diagram;
            		disconnectedP.add(g.copy());
            	}
            }
            p = makeGraphList();
            p.addAll(newP);
            if(isInteractive){
            	topSet.clear();
            	topSet.addAll(disconnectedP);
            	System.out.println("\ndisconnectedP");
            	ClusterViewer.createView("disconnectedP,before isofree", topSet);
            }
            // we don't need to re-isofree p, we know that's still good.
            // some of our new disconnected diagrams might condense with the old ones
    //        disconnectedP = isoFree.apply(maxIsomorph.apply(disconnectedP, mip), null);
            disconnectedP = maxIsomorph.apply(isoFree.apply(disconnectedP, null), mip);
            
            if(isInteractive){
            	topSet.clear();
            	topSet.addAll(newP);
            	System.out.println("\nnewP");
            	ClusterViewer.createView("newP in flexible ap branch", topSet);
            }
            if(isInteractive){
            	topSet.clear();
            	topSet.addAll(p);
            	System.out.println("\nP");
            	ClusterViewer.createView("p in flexible ap branch", topSet);
            }
            if(isInteractive){
            	topSet.clear();
            	topSet.addAll(disconnectedP);
            	System.out.println("\ndisconnectedP");
            	ClusterViewer.createView("disconnectedP,after isofree", topSet);
            }
             
            // all connected diagrams have been treated, now handle disconnected diagrams, refer to VirialDiagrams2 line 1551
            // disconnected graphs with i-1 components
            @SuppressWarnings("unchecked")
            Set<Graph>[] newDisconnectedP = new HashSet[n+1];
            for (int i=0; i<n+1; i++) {
            	newDisconnectedP[i] = new HashSet<Graph>();
            }
            SplitGraph graphSplitter = new SplitGraph();
            for (Graph g : disconnectedP) {
            	Set<Graph> gSplit = graphSplitter.apply(g);
            	newDisconnectedP[gSplit.size()].add(g);
            }
       
            disconnectedP.clear();
            for (int i = 0; i<n; i++) {
            	// looking for graphs with i components
            	for (Graph g : newDisconnectedP[i]) {
            		boolean ap = hap.check(g);
            		if (ap) {
            			// this is a disconnected diagram with a singly-connected component.
            			// we need to match it up with a more disconnected diagram
            			fop = null;
            			//loop over all articulation points in g to find where to factor the diagrams
            			for( Node node : g.nodes() ) {
            				byte nodeID = node.getId();
            				boolean isArticulationPt = hap.getArticulationPoints().contains(nodeID);//check whether the node is an a.p.
            				if (isArticulationPt){
            					fop = new FactorOnceParameters(nodeID, false);
            					break;// find 1 a.p. and we are good!
            				}
            				
            			}
            			Set<Graph> gfSet = factor.apply(g, fop);
            			Graph gf = gfSet.iterator().next();//only 1 iterate
            			gf = mulScalar.apply(gf, msp);
            			cancelMap.put(g,gf);
            			newDisconnectedP[i+1].addAll(gfSet);
            		}
            	}
            	
            	newDisconnectedP[i+1] = maxIsomorph.apply(isoFree.apply(newDisconnectedP[i+1], null), mip);
            	//newDisconnectedP[i+1] = isoFree.apply(maxIsomorph.apply(newDisconnectedP[i+1], mip), null);
            	disconnectedP.addAll(newDisconnectedP[i]);
            }
            
        }
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(disconnectedP);
        	System.out.println("\ndisconnectedP");
    	    ClusterViewer.createView("disconnectedP,after process C", topSet);
        }
        Set<Graph> cancel = new HashSet<Graph>();
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(p);
        	topSet.addAll(disconnectedP);
        	System.out.println("\nP");
        	for (Graph g : topSet) {
        		System.out.println(g);
        		Graph cancelGraph = cancelMap.get(g);
        		if (cancelGraph != null) {
        			System.out.println("   "+cancelGraph);
        			cancel.add(cancelGraph);
        		}
        	}
        	ClusterViewer.createView("P diagrams", topSet);
        }
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(cancel);
        	ClusterViewer.createView("cancelMap diagrams", topSet);
        }
    }// end makeVirialDiagrams method
    
    // add property class
    public static final class MaxIsomorphCriteriaMixture implements Property {
    	protected final HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    	public MaxIsomorphCriteriaMixture() {
    	
    	}
    	public boolean check(Graph graph) {
    		char prevColor = graph.getNode((byte)0).getColor();
    		for (byte i=1; i<graph.nodeCount(); i++) {
    			char color = graph.getNode(i).getColor();
    			if (color < prevColor) return false;
    			prevColor = color;
    		}
    		boolean ap = hap.check(graph);
    		if (!ap) return true;
    		prevColor = '0';
    		if (ap){//either A or B is a.p.
    			for (byte i=0; i<graph.nodeCount(); i++) {
    				char color = graph.getNode(i).getColor();
    				if (hap.getArticulationPoints().contains(i)) {
    					return color != prevColor ;
    				}
    				// this is equivalent with the if loop above
    				//		boolean isArticulationPt = hap.getArticulationPoints().contains(i);
    		//		if ( (color == prevColor) && isArticulationPt ) return false; // if still the 1st species, but it is a a.p., return false
    		//		if ( (color!= prevColor) && isArticulationPt) return true;// if the 2nd species has an a.p. and appears right after species A, accept this graph

    				prevColor = color;
				
    			}
    		}
    		return false;
    	}
    }// end of inner class "MaxIsomorphCriteria"
    

    
    public Set<Graph> getMSMCGraphs(boolean connectedOnly,int flexID, int[] numPoints) {
    	if (p == null) {
    		makeVirialDiagrams();
    	}
    	GraphList allP = makeGraphList();
        Set<Graph> topSet = makeGraphList();

    	for (Graph g : p) {
        	int[] numSpecies = new int[nodeColors.length];
        	boolean correctComposition = true;
    		// add graph to allP if the composition is correct
            for (Node node : g.nodes()) {//loop over all nodes
            	for (int t=0; t<numSpecies.length;t++){// loop over all species
            		if(node.getColor()==nodeColors[t]){
            			numSpecies[t]++;
            			break;
            		}
            	}
            }
            
            for (int t=0; t<numSpecies.length;t++){// loop over all species
            	if (numSpecies[t]!=numPoints[t]){
            		correctComposition = false;
            		break;
            	}
            }
            if(correctComposition){
            	allP.add(g);
            }
    	}
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(allP);
        	ClusterViewer.createView("allP,getMSMC,B"+numPoints[0]+numPoints[1]+" "+connectedOnly + " "+flexID, topSet);
        }
        // ==============================================================================================================//
        GraphList pn = makeGraphList();// add graph to pn based on flexID
        GraphList pnTemp = makeGraphList();
        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        for (Graph g : allP) {
        	boolean ap = hap.check(g);
        	if ( flexID==-1){
        		if (!ap){//rigid graph is added only==> biconnected graph, no corresponding graph in cancelMap
        			pn.add(g);
        		}
        		continue;
        	}
          
        	if(ap){// no biconnected graph
        		pnTemp.add(g);
        		//Flex A graph only, the ap must be species A
        	    //con or ! con are both fine ? Doesn't matter?????????????????????????????????????????????????
        		List<Byte> bytelist = hap.getArticulationPoints();
        		Collections.sort(bytelist);
        		for (byte nodeID : bytelist) {
        			char color = g.getNode(nodeID).getColor();
        			//System.out.println("node ID:"+ nodeID);
        			if (color == nodeColors[flexID]){//???????????????????????
        			//	System.out.println(" color == nodeColors[flexID] :"+ color );
        				pn.add(g);
        				if (!connectedOnly) {
        					Graph c = cancelMap.get(g);
        					pn.add(c);
        				}
        			}
        			else {
        		//		System.out.println("color != nodeColors[flexID]");
        			}
        			break;//don't need to check the next ap,if there is more than one ap in the graph
        		}
        	}
        		 	
        }
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(pnTemp);
        	ClusterViewer.createView("pnTemp, getMSMC, B"+numPoints[0]+numPoints[1]+" "+connectedOnly + " "+flexID, topSet);
        }
        if (isInteractive) {
        	topSet.clear();
        	topSet.addAll(pn);
        	ClusterViewer.createView("getMSMC, B"+numPoints[0]+numPoints[1]+" "+connectedOnly + " "+flexID, topSet);
        }
        return pn;
    }
    
    public Map<Graph,Graph> getCancelMap() {
        if (p == null) {
            makeVirialDiagrams();
        }
        return cancelMap;
    }
    
    // get graphs from getMSMC graphs
    public ClusterSum makeVirialCluster(MayerFunction[][] f, int flexID, int[] numPoints) {
        if (p == null) {
            makeVirialDiagrams();
        }
//        System.out.println("makeVirialCluster is executed!");
        Set<Graph> pn = getMSMCGraphs(false,flexID,numPoints);
        return makeVirialCluster(pn, f, flexID, numPoints);
    }
    
    public ClusterSum makeVirialCluster(Set<Graph> graphs, MayerFunction[][] f, int flexID, int[] numPoints) {
    	ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
        ArrayList<Double> weights = new ArrayList<Double>();
        // get firstPt and lastPt from flexID
        // flexID = -1 ==> rigid
        // flexID = 0  ==> flex A
        // flexID = 1  ==> flex B
        int firstPoint = 0;
        int lastPoint  = 0;
        for ( int i=0; i<flexID;i++){
        	firstPoint += numPoints[i];
        }
        if (flexID > -1){
        	lastPoint = firstPoint+numPoints[flexID]-1;
        }
        
        for (Graph g : graphs) {
        	boolean connected = (g.nodeCount() == n);
        	int nDiagrams = populateEFBonds(g, allBonds, weights, firstPoint,lastPoint,false, connected, flexID);
            if ( (nDiagrams>0)  && (flexColors.length>0) && (flexID!=-1) ) {
            	populateEFBonds(g, allBonds, weights, firstPoint, lastPoint, true, connected,flexID);
            }
        }

        double[] w = new double[weights.size()];
        for (int i=0; i<w.length; i++) {
            w[i] = weights.get(i);
        }
        // pour f in 1D array
        //add more information here
        int newfLength = (f.length + 1 ) * f.length / 2; 
        MayerFunction[] newf = new MayerFunction[newfLength];
        int k = 0; 
        for ( int i = 0 ; i < f.length; i++) {
        	for ( int j = i ; j < f.length ; j++){
        		newf[k] = f[i][j];
        		k++;
        	}
        }
        return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, newf);
       
    }// end makeVirialCluster
    
    protected byte swap0n(byte i, int firstPoint, int lastPoint, boolean connectedGraph) {
        if (connectedGraph){
        	if (i==firstPoint) return (byte)(lastPoint+1);
        	if (i>lastPoint) return (byte)(i+1);
        }
        else {
        	if (i==firstPoint) return (byte)(lastPoint+1);
            if (i==(lastPoint+1)) return (byte)firstPoint;
            
        }
        return i;
    }
    public int populateEFBonds(Graph g, List<ClusterBonds> allBonds, List<Double> weights, int firstPoint, int lastPoint, 
    		                   boolean swapControl, boolean connectedGraph, int flexID) {
        int rv = 0;
        int numSpecies = nodeColors.length;//number of species
        ArrayList<int[]>[] ebonds = new ArrayList[numSpecies*(numSpecies+1)/2];
        ArrayList<int[]>[] fbonds = new ArrayList[numSpecies*(numSpecies+1)/2];
        for ( int t=0;t<ebonds.length;t++){
        	ebonds[t]=new ArrayList<int[]>();
        	fbonds[t]=new ArrayList<int[]>();
        }
        rv = 1;
        for (Node node1 : g.nodes()) {
        	for (Node node2 : g.nodes()) {
        		if (node1.getId() >= node2.getId()) continue;
        		if (g.hasEdge(node1.getId(), node2.getId())) {
        			byte n1 = node1.getId();
        			byte n2 = node2.getId();
        			// get t1 and t2, the species number of each end
        			int t1=-1,t2=-1;
     	   			for ( int j=0;j<numSpecies; j++){
     	   		        if (node1.getColor() == nodeColors[j]){
     	   		        	t1 = j;
     	   		        }
     	   		        if (node2.getColor() == nodeColors[j]){
     	   		        	t2 = j;
     	   		        }
     	   			}
    			//	System.out.println("t1 and t2: "+t1+" "+t2);
    			//	System.out.println("node1 and node2: "+node1+" "+node2+" "+nodeColors[0]+" "+nodeColors[1]);

     	   			if ( (t1<0)||(t2<0) ){
     	   				throw new RuntimeException("wrong t1 or t2 value! t1:"+t1+",t2:"+t2);
     	   			}
     	   			
     	   			char edgeColor = g.getEdge(n1, n2).getColor();
        			if (swapControl) n1 = swap0n(n1,firstPoint, lastPoint,connectedGraph);
        			if (swapControl) n2 = swap0n(n2,firstPoint, lastPoint, connectedGraph);
//        			if ( (n1==3) || (n2==3) ){
//        				throw new RuntimeException("n1 or n2 = 3!");
//        			}
        			int index = t1*(2*numSpecies-t1+1)/2 +t2-t1;
        	        if ( (index<0) || (index>(fbonds.length-1)) ){
        	        	throw new RuntimeException("wrong index!, index is:"+index);
        	        }
        			if (edgeColor == fBond) {
//        				System.out.println("n1:"+n1+",n2:"+n2+",t1:"+t1+",t2:"+t2 +",index is:"+index);
//        				System.out.println("fbonds["+index+"] is: "+fbonds[index]);
        				fbonds[index].add(new int[]{n1,n2});
//        				System.out.println("fbonds["+index+"] is: "+fbonds[index]);
        			}
        			
        			else if (edgeColor == eBond) {
//        				System.out.println("n1:"+n1+",n2:"+n2+",t1:"+t1+",t2:"+t2 +",index is:"+index);
//        				System.out.println("ebonds["+index+"] is: "+fbonds[index]);
        				ebonds[index].add(new int[]{n1,n2});
//        				System.out.println("ebonds["+index+"] is: "+fbonds[index]);

        			}
        			else {
        				throw new RuntimeException("oops, unknown bond,edgeColor is: "+edgeColor);
        			}
        		}
        	}
        }
        
		int[][][] efBonds = new int[( (n>3&&doReeHoover)?2:1) * fbonds.length][0][0];//???????????????????
        for (int t=0;t<fbonds.length;t++){
    		efBonds[t]= fbonds[t].toArray(new int[0][0]);
    		if (ebonds[t].size()>0){
    			efBonds[t+fbonds.length] = ebonds[t].toArray(new int[0][0]);
    		}
        }
        allBonds.add(new ClusterBonds(  ( (flexColors.length>0) && (flexID!=-1) ) ? n+1 : n,efBonds) );

        double w = g.coefficient().getValue();
        if (flexColors.length> 0){
        	w *= 0.5;
        }
        weights.add(w);
            
        return rv;
    } //end of populateEFBonds

    public ClusterSumShell[] makeSingleVirialClusters(ClusterSum coreCluster, MayerFunction[][] f, int flexID, int[] numPoints) {
        if (p == null) {
            makeVirialDiagrams();
        }
        ArrayList<ClusterSumShell> allClusters = new ArrayList<ClusterSumShell>();
        Set<Graph> pn = getMSMCGraphs(true ,flexID, numPoints);
//        System.out.println("makeSingleVirialClusters is executed!");

        int firstPoint = 0;
        int lastPoint  = 0;
        for ( int i=0; i<flexID;i++){
        	firstPoint += numPoints[i];
        }
        if (flexID > -1){
        	lastPoint = firstPoint+numPoints[flexID]-1;
        }
        for (Graph g : pn) {
        	List<Double> weights = new ArrayList<Double>();
            ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            populateEFBonds(g, allBonds, weights, firstPoint, lastPoint, false, true, flexID);//this will be always executed 
            
            if ( (flexColors.length>0)  && (flexID !=-1)) {
                populateEFBonds(g, allBonds, weights, firstPoint, lastPoint,true, true, flexID);
            
                if (cancelMap.get(g) != null) {
                	Graph cg = cancelMap.get(g);
                	populateEFBonds(cg, allBonds, weights, firstPoint, lastPoint, false, false, flexID);
                	populateEFBonds(cg, allBonds, weights, firstPoint, lastPoint, true, false, flexID);
                }
            }
            int newfLength = (f.length + 1 ) * f.length / 2; 
            MayerFunction[] newf = new MayerFunction[newfLength];
            int k = 0; 
            for ( int i = 0 ; i < f.length; i++) {
            	for ( int j = i ; j < f.length ; j++){
            		newf[k] = f[i][j];
            		k++;
            	}
            }
            double[] thisW = new double[weights.size()];
            double gCoef = g.coefficient().getValue();
            for (int i=0; i<thisW.length; i++) {
                thisW[i] = weights.get(i) / gCoef;
            }
            allClusters.add(new ClusterSumShell(coreCluster, allBonds.toArray(new ClusterBonds[0]), thisW, newf));// newf
        }
        return allClusters.toArray(new ClusterSumShell[0]);
    }// end of makeSingleVirialClusters
    
}
