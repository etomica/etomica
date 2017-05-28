/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.operations.Exclude;
import etomica.graph.operations.ExcludeParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.property.HasNoRootEdge;
import etomica.graph.property.IsBiconnected;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterBonds;
import etomica.virial.ClusterSum;
import etomica.virial.MayerFunction;

/**
 * @author Hye Min Kim
 *
 * Class that provides some standard pair sets  used in specification of clusters for 
 * having association bonds.
 */
public final class ClusterAssociation {

    public static ClusterAbstract virialCluster(int whichGraph, MayerFunction fR, MayerFunction fAA, 
             byte numPoint) {
    	
        	Set<Graph> setOfBiConnected = new HashSet<Graph>();//set of doubly connected diagrams with f bonds
          Iterator<Graph> iterator = new IsomorphismFilter(new PropertyFilter(new PropertyFilter(new DefaultIterator(numPoint),
              new HasNoRootEdge()), new IsBiconnected()));
          while (iterator.hasNext()) {
                setOfBiConnected.add(iterator.next());
            }
        char fTotal = 'A';//color of edge
        char fRef = 'R';
        char fAssociation = 'F';
    	  
    	  ExcludeParameters excludeParameters = new ExcludeParameters();
    	  
    	  excludeParameters.bondMap.put(fAssociation, new char[]{'A','A'});
          
      Split split = new Split();
      SplitParameters bonds = new SplitParameters(fTotal,fRef,fAssociation);
      Set<Graph> setOfSubstituted = split.apply(setOfBiConnected, bonds);
      Exclude exclude = new Exclude();
      Set<Graph> result = exclude.apply(setOfSubstituted, excludeParameters);
      Graph graph = (Graph)result.toArray()[whichGraph];

    	int nBody = graph.nodeCount();//number of nodes in a diagram
        int numBonds = graph.edgeCount();//number of edges in a diagram
        List<Edge> edge = graph.edges();
        int numfRBonds = 0;
        for (int j = 0; j < edge.size(); j++){
        if (edge.get(j).getColor() == fRef){
            numfRBonds++;	
        	}
        }

            int[][][] bondList = new int[2][][];//fR = 0, fAA = 1
            int numfAABonds = numBonds - numfRBonds;
            bondList[0] = new int[numfRBonds][2];
            bondList[1] = new int[numfAABonds][2];
            int k = 0;
            int l = 0;
            for (Node node1 : graph.nodes()) {
    			for (Node node2 : graph.nodes()){
    				if (node1 == node2 || !graph.hasEdge(node1.getId(), node2.getId())){
    					continue;
    				}
    				Edge edge12 = graph.getEdge(node1.getId(), node2.getId());
    				if (edge12.getColor()== fRef){
    					
    					bondList[0][k][0] = node1.getId();
    					bondList[0][k][1] = node2.getId();
    					k++;
    				} else {
    					bondList[1][l][0] = node1.getId();
    					bondList[1][l][1] = node2.getId();
    					l++;
    				}
    			}
            }

        ClusterBonds[] clusters = new ClusterBonds[]{new ClusterBonds(nBody, bondList, false)};

        return new ClusterSum(clusters,new double[]{graph.coefficient().getNumerator()/(double)graph.coefficient().getDenominator()},new MayerFunction[]{fR, fAA});
    }
}
