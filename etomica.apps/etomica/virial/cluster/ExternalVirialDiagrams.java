package etomica.virial.cluster;

import static etomica.graph.model.Metadata.COLOR_CODE_0;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IdenticalGraphFilter;
import etomica.graph.iterators.filters.PropertyFilter;
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
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.MulScalar;
import etomica.graph.operations.MulScalarParameters;
import etomica.graph.operations.PCopy;
import etomica.graph.operations.ReeHoover;
import etomica.graph.operations.ReeHoover.ReeHooverParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitGraph;
import etomica.graph.operations.SplitOne;
import etomica.graph.operations.SplitOne.SplitOneParameters;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;
import etomica.graph.viewer.ClusterViewer;
import etomica.math.SpecialFunctions;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumEF;
import etomica.virial.ClusterSumExternalField;
import etomica.virial.ClusterSumShell;
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
