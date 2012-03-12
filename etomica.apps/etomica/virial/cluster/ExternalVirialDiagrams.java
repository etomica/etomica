package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumExternalField;
import etomica.virial.MayerFunction;

public class ExternalVirialDiagrams extends VirialDiagrams{

	public ExternalVirialDiagrams(int n, boolean multibody, boolean flex) {
		super(n, multibody, flex);
		// TODO Auto-generated constructor stub
	}
	 public ClusterSum makeRhoCluster(MayerFunction f, boolean useExternalField) {
	       if (rho == null) {
	           makeRhoDiagrams();
	       }
	        ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
	        ArrayList<Double> weights = new ArrayList<Double>();
	        Set<Graph> rhon = getMSMCGraphsEX(false);
	        for (Graph g : rhon) {
	            populateEFBonds(g, allBonds, false);
	            
	            double w = ((double)g.coefficient().getNumerator())/g.coefficient().getDenominator();
	           
	            weights.add(w);
	            
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
	        GraphList<Graph> rhon = makeGraphList();
	       
	        for (Graph g : rho) {
	            
	            
	            if (g.nodeCount()==n) {
	                rhon.add(g);
	            }
	        }
	        return rhon;
	    }

}      
