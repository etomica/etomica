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
 *
 * @author Hye Min Kim
 */
public class BondConnectedSubstitution3Site implements Unary {
	
	public Set<Graph> apply(Set<Graph> argument, Parameters params) {


	    assert (params instanceof BondConnectedSubstitutionParameters3Site);
	    Set<Graph> result = new HashSet<Graph>();
	    for (Graph g : argument) {
	    	Set<Graph> newGraph = apply(g, (BondConnectedSubstitutionParameters3Site) params);
	      result.addAll(newGraph);
	    }
	    return result;
	  }

	  public Set<Graph> apply(Graph g, BondConnectedSubstitutionParameters3Site params) {
	      char fRBond = params.fRBond;
	      char eRBond = params.eRBond;
	      char capFACBond = MetadataImpl.edgeColorPairs.get(0).get(0);//1st pair, 1st bond from 1st pair
	      char capFCABond = MetadataImpl.edgeColorPairs.get(0).get(1);//1st pair, 1st bond from 1st pair
	      char capFBCBond = MetadataImpl.edgeColorPairs.get(1).get(0);//1st pair, 1st bond from 1st pair
	      char capFCBBond = MetadataImpl.edgeColorPairs.get(1).get(1);//1st pair, 1st bond from 1st pair
	      
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
          Set<Node> nodeTraveledC = new HashSet<Node>();
          Set<Node> nodeTraveledACurrent = new HashSet<Node>();
          Set<Node> nodeTraveledBCurrent = new HashSet<Node>();
          Set<Node> nodeTraveledCCurrent = new HashSet<Node>();
          for (byte j = 0;j<(gOnlyF.nodeCount()-1);j++){//looping over nodes
          	Node node = gOnlyF.getNode(j);
          	for (char site = 'A';site <= 'C';site++){
          		boolean startA = site=='A';
        		boolean startB = site=='B';
        		boolean startC = site=='C';
               	if ((startA && !nodeTraveledA.contains(node))||(startB && !nodeTraveledB.contains(node))||(startC && !nodeTraveledC.contains(node))){
            		for (byte i=(byte)(j+1); i <gOnlyF.nodeCount();i++){
            			if(!gOnlyF.hasEdge(j, i))continue;
            			Edge edge = gOnlyF.getEdge(j, i); 
            			if ((edge.getColor()== capFACBond)&& startA){
        					nodeTraveledA.add(node);
        					nodeTraveledACurrent.add(node);
            			}
        				else if ((edge.getColor()== capFBCBond)&& startB) {
        					nodeTraveledB.add(node);
        					nodeTraveledBCurrent.add(node);
        				}
        				else if ((edge.getColor()== capFCABond||edge.getColor()== capFCBBond)&& startC) {
        					nodeTraveledC.add(node);
        					nodeTraveledCCurrent.add(node);
        				} else continue;
        				edgeMap.put(edge,true);//this edge goes forward
        				pool.add(edge);
            		}
               		while(!pool.isEmpty()){
            			Edge edge = pool.remove(pool.size()-1);
            			edgeTraveled.add(edge);
            			boolean forward = edgeMap.get(edge);
            			byte newNode = forward?gOnlyF.getToNode(edge.getId()):gOnlyF.getFromNode(edge.getId());
            			boolean siteA, siteB, siteC;
            			if (forward){
            				siteA = (edge.getColor() == capFCABond);
            				siteB = (edge.getColor() == capFCBBond);
            				siteC = (edge.getColor() == capFACBond||edge.getColor() == capFBCBond);
            			} 
            			else {
            				siteA = (edge.getColor() == capFACBond);
            				siteB = (edge.getColor() == capFBCBond);
            				siteC = (edge.getColor() == capFCABond||edge.getColor() == capFCBBond);
            			}
            			if(siteA) {
        					nodeTraveledA.add(gOnlyF.getNode(newNode));
        					nodeTraveledACurrent.add(gOnlyF.getNode(newNode));
        				}
            			else if(siteB) {
        					nodeTraveledB.add(gOnlyF.getNode(newNode));
        					nodeTraveledBCurrent.add(gOnlyF.getNode(newNode));
        				}
            			else if(siteC) {
        					nodeTraveledC.add(gOnlyF.getNode(newNode));
        					nodeTraveledCCurrent.add(gOnlyF.getNode(newNode));
        				}
                		for (byte i=0; i <gOnlyF.nodeCount();i++){
                			if (i == newNode)continue;
                			if(!gOnlyF.hasEdge(newNode, i))continue;
            				Edge newEdge = i>newNode?gOnlyF.getEdge(newNode, i):gOnlyF.getEdge(i,newNode);
            				if(edgeTraveled.contains(newEdge)){
            					continue;
            				}
            				boolean newForward = i>newNode;
            				boolean newSiteA, newSiteB, newSiteC;
            				if (newForward){
                				newSiteA = (newEdge.getColor() == capFACBond);
                				newSiteB = (newEdge.getColor() == capFBCBond);
                				newSiteC = (newEdge.getColor() == capFCABond||newEdge.getColor() == capFCBBond);

                			} 
                			else {
                				newSiteA = (newEdge.getColor() == capFCABond);
                				newSiteB = (newEdge.getColor() == capFCBBond);
                				newSiteC = (newEdge.getColor() == capFACBond||newEdge.getColor() == capFBCBond);
                			}
            				if ((siteA && newSiteA)||(siteB && newSiteB)||(siteC && newSiteC)){
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
                  for (Node nodeTraveled:nodeTraveledCCurrent){
                  	allNodes.add(nodeTraveled.getId());
                  }
                  if (!allNodes.isEmpty()){
                      components.add(allNodes);
                      nodeTraveledACurrent.clear();
                      nodeTraveledBCurrent.clear();
                      nodeTraveledCCurrent.clear();
                  }
          
          	}
          }
          boolean noFR = true;
          //System.out.println("gOnlyF "+gOnlyF+" components "+components);
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
	  public static class BondConnectedSubstitutionParameters3Site implements Parameters {
	    public final char fRBond, eRBond;
	    public BondConnectedSubstitutionParameters3Site (char fRBond, char eRBond){
	    	this.fRBond = fRBond;
	    	this.eRBond = eRBond;
	    }
	  }
}
