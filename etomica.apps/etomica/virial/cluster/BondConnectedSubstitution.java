/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graph.operations.Parameters;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.operations.Split;
import etomica.graph.operations.SplitParameters;
import etomica.graph.operations.Unary;
import etomica.graph.operations.DecorateWertheim2Site.DecorateWertheimParameters2Site;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * If Wertheim diagrams have bond connected network, fR bond is split into eR and -one bond.
 * This is for 2 site model.
 * @author Hye Min Kim
 */
public class BondConnectedSubstitution implements Unary {
	
	public Set<Graph> apply(Set<Graph> argument, Parameters params) {


	    assert (params instanceof BondConnectedSubstitutionParameters2Site);
	    Set<Graph> result = new HashSet<Graph>();
	    for (Graph g : argument) {
	    	Set<Graph> newGraph = apply(g, (BondConnectedSubstitutionParameters2Site) params);
	      result.addAll(newGraph);
	    }
	    return result;
	  }

	  public Set<Graph> apply(Graph g, BondConnectedSubstitutionParameters2Site params) {
	      char fRBond = params.fRBond;
	      char eRBond = params.eRBond;
	      char capFABBond = MetadataImpl.edgeColorPairs.get(0).get(0);
	      char capFBABond = MetadataImpl.edgeColorPairs.get(0).get(1);
	      Set<Graph> result = new HashSet<Graph>();

      	Graph gOnlyF=g.copy();
          DeleteEdgeParameters deleteEdgeParameters = new DeleteEdgeParameters(fRBond);
          DeleteEdge deleteEdge = new DeleteEdge();
          gOnlyF = deleteEdge.apply(gOnlyF, deleteEdgeParameters);
          List<List<Byte>> components = new ArrayList<List<Byte>>();
          Map<Edge,Boolean> edgeMap = new HashMap<Edge,Boolean>();
          ArrayList<Edge> pool = new ArrayList<Edge>();
          Set<Edge> edgeTraveled = new HashSet<Edge>();
          Set<Node> nodeTraveledA = new HashSet<Node>();
          Set<Node> nodeTraveledB = new HashSet<Node>();
          Set<Node> nodeTraveledACurrent = new HashSet<Node>();
          Set<Node> nodeTraveledBCurrent = new HashSet<Node>();
          for (byte j = 0;j<(gOnlyF.nodeCount()-1);j++){
            	Node node = gOnlyF.getNode(j);
            	for (char site = 'A';site <= 'B';site++){
            		boolean startA = site=='A';
            		boolean startB = site=='B';
                 	if ((startA && !nodeTraveledA.contains(node))||(startB && !nodeTraveledB.contains(node))){
              		for (byte i=(byte)(j+1); i <gOnlyF.nodeCount();i++){
              			if(!gOnlyF.hasEdge(j, i))continue;
              			Edge edge = gOnlyF.getEdge(j, i); 
              			if ((edge.getColor()== capFABBond)&& startA){
          					nodeTraveledA.add(node);
          					nodeTraveledACurrent.add(node);
              			}
          				else if ((edge.getColor()== capFBABond)&& startB) {
          					nodeTraveledB.add(node);
          					nodeTraveledBCurrent.add(node);
          				} else continue;
          				edgeMap.put(edge,true);//this edge goes forward
          				pool.add(edge);
              		}
                 		while(!pool.isEmpty()){
              			Edge edge = pool.remove(pool.size()-1);
              			edgeTraveled.add(edge);
              			boolean forward = edgeMap.get(edge);
              			byte newNode = forward?gOnlyF.getToNode(edge.getId()):gOnlyF.getFromNode(edge.getId());
              			boolean siteA, siteB;
              			if (forward){
              				siteA = (edge.getColor() == capFBABond);
              				siteB = (edge.getColor() == capFABBond);
              			} 
              			else {
              				siteA = (edge.getColor() == capFABBond);
              				siteB = (edge.getColor() == capFBABond);
              			}
              			if(siteA) {
          					nodeTraveledA.add(gOnlyF.getNode(newNode));
          					nodeTraveledACurrent.add(gOnlyF.getNode(newNode));
          				}
              			else if(siteB) {
          					nodeTraveledB.add(gOnlyF.getNode(newNode));
          					nodeTraveledBCurrent.add(gOnlyF.getNode(newNode));
          				}
                  		for (byte i=0; i <gOnlyF.nodeCount();i++){
                  			if (i == newNode)continue;
                  			if(!gOnlyF.hasEdge(newNode, i))continue;
              				Edge newEdge = i>newNode?gOnlyF.getEdge(newNode, i):gOnlyF.getEdge(i,newNode);
              				if(edgeTraveled.contains(newEdge)){
              					continue;
              				}
              				boolean newForward = i>newNode;
              				boolean newSiteA, newSiteB;
              				if (newForward){
                  				newSiteA = (newEdge.getColor() == capFABBond);
                  				newSiteB = (newEdge.getColor() == capFBABond);

                  			} 
                  			else {
                  				newSiteA = (newEdge.getColor() == capFBABond);
                  				newSiteB = (newEdge.getColor() == capFABBond);
                  			}
              				if ((siteA && newSiteA)||(siteB && newSiteB)){
              					if(pool.contains(newEdge)){
                  					pool.remove(newEdge);
                  					continue;
                  				}

                  				if(i>newNode){
                  					edgeMap.put(newEdge, true);
                  				} else edgeMap.put(newEdge, false);
                  				pool.add(newEdge);
              				}
                  		}
                  		
              		}
                	}
                    List<Byte> allNodes = new ArrayList<Byte>();	
                    for (Node nodeTraveled:nodeTraveledACurrent){
                    	allNodes.add(nodeTraveled.getId());
                    }
                    for (Node nodeTraveled:nodeTraveledBCurrent){
                    	allNodes.add(nodeTraveled.getId());
                    }
                    if (!allNodes.isEmpty()){
                        components.add(allNodes);
                        nodeTraveledACurrent.clear();
                        nodeTraveledBCurrent.clear();
                    }
            
            	}
            }
          boolean noFR = true;

          for (int iComp = 0; iComp<components.size(); iComp++) {
              List<Byte> nodes = components.get(iComp);
              if (nodes.size() < 3)continue;
              for (byte node1:nodes){
              	for (byte node2:nodes){
              		if(node2<= node1)continue;
              		if(g.hasEdge(node1, node2)){
              			if (g.getEdge(node1, node2).getColor()== fRBond){
              				g.getEdge(node1, node2).setColor('Z');
              				noFR = false;
              			}
              		}
              	}
              }
          }
          if(!noFR){
        	  char oneBond = '1';
                Split split = new Split();
                SplitParameters bonds = new SplitParameters('Z', eRBond, oneBond);
                Set<Graph> setOfSubstituted = split.apply(g, bonds);
                for (Graph gSet:setOfSubstituted){
                	int count = 0;
                	for (Edge edge:gSet.edges()){
                		if (edge.getColor() == oneBond)count++;
                	}
                	if (count%2==1)gSet.coefficient().multiply(new CoefficientImpl(-1));
                }
                DeleteEdgeParameters deleteOneEdgeParameters = new DeleteEdgeParameters(oneBond);
                setOfSubstituted = deleteEdge.apply(setOfSubstituted, deleteOneEdgeParameters);
                result.addAll(setOfSubstituted);
          } else result.add(g);

		return result;
	}
	  public static class BondConnectedSubstitutionParameters2Site implements Parameters {
	    public final char fRBond, eRBond;
	    public BondConnectedSubstitutionParameters2Site (char fRBond, char eRBond){
	    	this.fRBond = fRBond;
	    	this.eRBond = eRBond;
	    }
	  }
}
