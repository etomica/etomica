/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.GraphOp;
import etomica.graph.operations.IsoFree;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.viewer.ClusterViewer;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumExternalField;
import etomica.virial.MayerFunction;

public class ExternalVirialDiagrams extends VirialDiagrams{

	protected boolean doSurfacetension;
	protected boolean doExpDensityprofile;
	public static void main(String[] args) {
    	MetadataImpl.rootPointsSpecial = true;
        final int n = 4;
        ExternalVirialDiagrams virialDiagrams = new ExternalVirialDiagrams(n, true, false, true);
        virialDiagrams.setDoShortcut(true);
        virialDiagrams.makeRhoDiagrams();
    }
	public ExternalVirialDiagrams(int n, boolean doSurfacetension, boolean doExpDensityprofile, boolean interactive) {
		super(n, false, false, interactive);
		this.doSurfacetension = doSurfacetension;
		this.doExpDensityprofile = doExpDensityprofile;
		
	}
	 public ClusterSum makeRhoCluster(MayerFunction f, boolean useExternalField) {
	       if (rho == null) {
	           makeRhoDiagrams();
	       }
	        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
	        ArrayList<Double> weights = new ArrayList<Double>();
	        Set<Graph> rhon = getMSMCGraphsEX(false);
	        //System.out.println(rhon.size());
	        for (Graph g : rhon) {
	            populateEFBonds(g, allBonds, weights, false);
	        }
	        double[] w = new double[weights.size()];
	        for (int i=0; i<w.length; i++) {
	            w[i] = weights.get(i);
	        }
	        if (useExternalField) {
		        return new ClusterSumExternalField(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f});
		    }
		    else {
		        return new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f});
		    }
	        
	    }
	 public Set<Graph> getMSMCGraphsEX(boolean connectedOnly) {
	        if (rho == null) {
	            makeRhoDiagrams();
	        }
	        GraphList rhon = makeGraphList();
	       
	        for (Graph g : rho) {
	            
	            
	            if (g.nodeCount()==n) {
	                rhon.add(g);
	            }
	        }
	        return rhon;
	    }
	@Override
	public void makeRhoDiagrams() {
		// TODO Auto-generated method stub
		super.makeRhoDiagrams();
		//System.out.println(rho.size());
//        Set<Graph> fXipow = new HashSet<Graph>();
//        Set<Graph> excessRho = new HashSet<Graph>();
//        for(Graph g : rho){
//        	if(g.nodeCount()>1){
//        		excessRho.add(g);
//        	}
//        }
//        fXipow.addAll(excessRho);
//        IsoFree isoFree = new IsoFree();
//        MulScalarParameters msp = null;
//        MulScalar mulScalar = new MulScalar();
//
//        MulFlexible mulFlex = new MulFlexible();
//        MulFlexibleParameters mfp = MulFlexibleParameters.makeParameters(flexColors, (byte)(n-1));
//        Set<Graph> logRho = new HashSet<Graph>();
//        for (int i=1; i<n+1; i++) {
//
//        	logRho.addAll(fXipow);
//        	logRho = isoFree.apply(logRho, null);
//            msp = new MulScalarParameters(new CoefficientImpl(-i,(i+1)));
//            fXipow = isoFree.apply(mulScalar.apply(mulFlex.apply(fXipow, excessRho, mfp), msp), null);
//        }

//        if (isInteractive) {
//        	GraphList<Graph> topSet = makeGraphList();
//            System.out.println("\nlogRho");
//            topSet.clear();
//            topSet.addAll(logRho);
//            for (Graph g : topSet) {
//                System.out.println(g);
//            }
//            ClusterViewer.createView("logRho", topSet);
//        }
        if (doSurfacetension) {
        	Set<Graph> rhoBi = new HashSet<Graph>();
        	Set<Graph> rhoNotBi = new HashSet<Graph>();
        	IsBiconnected isBi = new IsBiconnected();
        	for (Graph g :rho){        		 
        		for (Node node : g.nodes()){
        			if (node.getType()==Metadata.TYPE_NODE_ROOT){
        				node.setType(Metadata.TYPE_NODE_FIELD);
        				break;
        			}
        		}
        	      if (isBi.check(g)){
        	    	  rhoBi.add(g);
        	      }
        	      else{
        	    	  rhoNotBi.add(g);
        	    	  }
        	}
        	
	        MaxIsomorph maxIsomorph = new MaxIsomorph();
	        MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), new ArticulatedAt0(false));
	        IsoFree isoFree = new IsoFree();
	        
	        rho = isoFree.apply(maxIsomorph.apply(rhoNotBi, mip), null);
	        rho.addAll(isoFree.apply(maxIsomorph.apply(rhoBi, MaxIsomorph.PARAM_ALL), null));
	        for ( Graph g : rho){
	        	g.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
	        	
	        }
	        if (isInteractive) {
	        	GraphList topSet = makeGraphList();
	            System.out.println("\nsurface tension");
	            topSet.clear();
	            topSet.addAll(rho);
	            for (Graph g : topSet) {
	                System.out.println(g);
	            }
	            ClusterViewer.createView("surface tension", topSet);
	        }
        }
        else if (doExpDensityprofile){
        	HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            Set<Graph> rhoNew = new HashSet<Graph>();
            for (Graph g : rho) {
                if (!hap.check(g) || !hap.getArticulationPoints().contains((byte)0)) {
                    rhoNew.add(g);
                    continue;
                }              
                
                
              
            }
	        rho = rhoNew;
	        if (isInteractive) {
	        	GraphList topSet = makeGraphList();
	            System.out.println("\nExpDensityprofile");
	            topSet.clear();
	            topSet.addAll(rho);
	            for (Graph g : topSet) {
	                System.out.println(g);
	            }
	            ClusterViewer.createView("ExpDensityprofile", topSet);
	        }
        }
        
       
	}

}      
