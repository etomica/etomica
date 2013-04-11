package etomica.virial.cluster;

import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsConnected;
import etomica.graph.property.Property;
import etomica.graph.traversal.CVisitor;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.virial.cluster.VirialDiagrams.ArticulatedAt0;

public class VirialDiagramsMix2 {

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
    protected boolean doDisconnectedMatching = true;
    protected final char nodeColor = Metadata.COLOR_CODE_0;
    protected char[] flexColors;
    public char fBond, eBond, excBond, mBond, MBond, efbcBond;
    
    public static void main(String[] args) {
        int n = 4;
        boolean[] multibody = new boolean[]{false,false};
        boolean[] flex = new boolean[]{true, false};// 1st species:black; 2nd species:red
        boolean doKeepEBonds = false;
        boolean doReeHoover = false;
        VirialDiagramsMix2 virialDiagrams = new VirialDiagramsMix2(n, multibody, flex, true);
        virialDiagrams.setDoReeHoover(doReeHoover);
        virialDiagrams.setDoKeepEBonds(doKeepEBonds);
        virialDiagrams.setDoShortcut(false);
        virialDiagrams.setDoMinimalMulti(false);
        virialDiagrams.setDoMinimalBC(false);
        virialDiagrams.makeVirialDiagrams();
    }
    
    public VirialDiagramsMix2(int n, boolean[] multibody, boolean[] flex) {
        this(n, multibody, flex, false);
    }

    public VirialDiagramsMix2(int n, boolean[] multibody, boolean[] flex, boolean interactive) {
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
    public void setDoDisconnectedMatching(boolean newDoDisconnectedMatching) {
        doDisconnectedMatching = newDoDisconnectedMatching;
    }
    public Set<Graph> makeGraphList() {
        ComparatorChain comp = new ComparatorChain();
        comp.addComparator(new ComparatorNumFieldNodes());
        comp.addComparator(new ComparatorBiConnected());
        comp.addComparator(new ComparatorNodeColors(nodeColors));
        comp.addComparator(new ComparatorNumEdges());
        comp.addComparator(new ComparatorNumNodes());
        return new GraphList<Graph>(comp);
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
        // ==================================================== eXi ======================================================= //
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
            }
        }
        topSet.addAll(eXi);
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("eXi", topSet);
        // ==================================================== fXi ======================================================= //
        Split split = new Split();
        SplitParameters bonds = new SplitParameters(eBond, fBond, oneBond);
        Set<Graph> setOfSubstituted = split.apply(eXi, bonds);
        DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(oneBond);
        DeleteEdge deleteEdge = new DeleteEdge();
        System.out.println("\nXi with f bonds");//set of full star diagrams with f bonds
        Set<Graph> fXi = deleteEdge.apply(setOfSubstituted, deleteEdgeParameters);
        IsoFree isoFree = new IsoFree();
        fXi = isoFree.apply(fXi, null);
        ClusterViewer.createView("fXi", fXi);
        topSet.clear();
        topSet.addAll(fXi);
        for (Graph g : topSet) {
            System.out.println(g);
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
        topSet.clear();
        topSet.addAll(lnfXi);
        System.out.println("\nlnfXi");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("lnfXi", lnfXi);
        
        
         // put lnfXi ordered
        if (false) {
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
        
        topSet.clear();
        topSet.addAll(rhoA);
        ClusterViewer.createView("rhoA", topSet);
        System.out.println("\nrhoA");
        for (Graph g : topSet) {
            System.out.println(g);
        }

        DifByNode opzdlnXidzB = new DifByNode();
        difParams = new DifParameters('B');
        rhoB = isoFree.apply(opzdlnXidzB.apply(lnfXi, difParams), null);

        topSet.clear();
        topSet.addAll(rhoB);
        ClusterViewer.createView("rhoB", topSet);
        System.out.println("\nrhoB");
        for (Graph g : topSet) {
            System.out.println(g);
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
     

        MulFlexibleParameters mfpn = MulFlexibleParameters.makeParameters(nodeColors, (byte)n);///????????????????????
        MaxIsomorph maxIsomorph = new MaxIsomorph();
        // refer to VirialDiagram2 Line 1457
        //compare Property happyArticulation = new ArticulatedAt0(doExchange, multibody ? mmBond : '0');
 //       Property happyArticulation = new ArticulatedAt0(false, '0');// articulation point is at 0
//        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), happyArticulation);
 //       HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();//implements property
        Property mic = new MaxIsomorphCriteria();
        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOpMaxRoot(), mic);
        HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
        
        //----------------------- decorate P using zA and zB -----------------------------------------//
        p = decorate.apply(lnfXi, zA, new DecorateParameters(0, mfpn));
        p = decorate.apply(p, zB, new DecorateParameters(1, mfpn));
        p = isoFree.apply(maxIsomorph.apply(p, mip), null);// based on maxisomorphCriteria
        
        topSet.clear();
        topSet.addAll(p);
        System.out.println("\nfirst P");
        for (Graph g : topSet) {
            System.out.println(g);
        }
        ClusterViewer.createView("first P", topSet);
     
        Set<Graph> newP = new HashSet<Graph>();// refer to line 1455
        disconnectedP = new HashSet<Graph>();// attempt to factor any graphs with an articulation point
        cancelMap = new HashMap<Graph,Graph>();
        
        //----------------------- factor diagrams with rigid articulation points-----------------------------------------//
        if (flexColors.length < nodeColors.length) {//at least 1 species is rigid, refer to line 1459, !flex
        	Factor factor = new Factor();
        	MulFlexibleParameters factorParameters = MulFlexibleParameters.makeParameters(flexColors, (byte)n);
        	for (Graph g : p) {
        		boolean ap = hap.check(g);
        		boolean con = hap.isConnected();
        		if ((con && ap) || (!con && hap.getArticulationPoints().size() > 0)) {//?????????????????
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
        }
       //----------------------- factor diagrams with flex articulation points-----------------------------------------//
        if (flexColors.length > 0) {  // refer to VirialDiagram2 Line 1489
            // pretend everything is fully flexible
            FactorOnce factor = new FactorOnce();
            FactorOnceParameters fop = null;

            // match up singly-connected (in p) with disconnected diagrams.
            // we have to do this last so that our cancelMap remains valid.
            newP.clear();
            msp = new MulScalarParameters(-1, 1);//????????????

            for (Graph g : p) {
            	boolean ap = hap.check(g);
            	boolean con = hap.isConnected();
            	if (con && ap ) {
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
            		
            		// newP will contain connected diagrams
            		g = g.copy();//??????????????????????refer to ViralDiagram2 Line 1509
            		newP.add(g);// newP now has singly-connected diagrams ONLY
            		Set<Graph> gfSet = factor.apply(g, fop);// factored diagrams set of g
            		Graph gf = gfSet.iterator().next();// only 1 iterate
            		disconnectedP.add(gf);// disconnectedP has factored diagrams///?????????????????/ 
            		gf = mulScalar.apply(gf, msp);
            		cancelMap.put(g, gf);
            	}
            	else if (con) {////??????????????????
            		// this is a biconnected diagram;
            		newP.add(g.copy());
            	}
            	else {
            		// this is a disconnected diagram;
            		disconnectedP.add(g.copy());
            	}
            }
            p = makeGraphList();
            p.addAll(newP);
            
            // we don't need to re-isofree p, we know that's still good.
            // some of our new disconnected diagrams might condense with the old ones
            disconnectedP = isoFree.apply(maxIsomorph.apply(disconnectedP, mip), null);
            	   
            topSet.clear();
            topSet.addAll(disconnectedP);
            System.out.println("\ndisconnectedP");
            ClusterViewer.createView("disconnectedP", topSet);

             
            // refer to VirialDiagram line 1541
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
            			// all of the articulation points we're interested in will be at 0
            			// (because of happyArticulation.  Other articulation points will be
            			// from exchange groups, but we can't split those.
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
            			Graph gf = gfSet.iterator().next();
            			gf = mulScalar.apply(gf, msp);
            			cancelMap.put(g,gf);
            			newDisconnectedP[i+1].addAll(gfSet);
            		}
            	}
            	
            	newDisconnectedP[i+1] = isoFree.apply(maxIsomorph.apply(newDisconnectedP[i+1], mip), null);

            	disconnectedP.addAll(newDisconnectedP[i]);
            }
            

        }
        
        topSet.clear();
        topSet.addAll(p);
        topSet.addAll(disconnectedP);
        System.out.println("\nP");
        for (Graph g : topSet) {
            System.out.println(g);
            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                System.out.println("   "+cancelGraph);
                }
            }
        
        ClusterViewer.createView("P, all connected diagrams", topSet);
    }// end makeVirialDiagrams method
    
    // add property class
    public static final class MaxIsomorphCriteria implements Property {
    	protected final HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    	public MaxIsomorphCriteria() {
    	
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
    
}
