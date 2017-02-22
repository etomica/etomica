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
 * revised from Hye Min's "BondConnectedSubstitution"
 * If Wertheim diagrams have bond connected network, fR bond is split into eR and -one bond.
 * Naphthalene 2 site model.
 * @author shu
 * Date:Sep 10 2011
 * 
 */
public class BondConnectedSubstitutionNaphthalene implements Unary {
	
	public Set<Graph> apply(Set<Graph> argument, Parameters params) {


	    assert (params instanceof BondConnectedSubstitutionParameters2SiteNa);
	    Set<Graph> result = new HashSet<Graph>();
	    for (Graph g : argument) {
	    	Set<Graph> newGraph = apply(g, (BondConnectedSubstitutionParameters2SiteNa) params);
	      result.addAll(newGraph);
	    }
	    return result;
	  }

	public Set<Graph> apply(Graph g, BondConnectedSubstitutionParameters2SiteNa params) {
		  char fRBond = params.fRBond;
	      char eRBond = params.eRBond;
	      List<Character> characterStartA = params.characterStartA;
	      List<Character> characterEndA = params.characterEndA;
	      List<Character> characterStartB = params.characterStartB;
	      List<Character> characterEndB = params.characterEndB;

	                    
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
          Set<Node> nodeTraveledCurrent = new HashSet<Node>();
          for (byte j = 0;j<(gOnlyF.nodeCount()-1);j++){
        	  Node node = gOnlyF.getNode(j);
        	  for (char site = 'A';site <= 'B';site++){
        		  boolean startA = site =='A';
        		  boolean startB = site =='B';
        		  if ((startA && !nodeTraveledA.contains(node))||(startB && !nodeTraveledB.contains(node))){
        			  for (byte i=(byte)(j+1); i <gOnlyF.nodeCount();i++){
        				  if(!gOnlyF.hasEdge(j, i))continue;
        				  Edge edge = gOnlyF.getEdge(j, i); 
        				  if (startA  && characterStartA.contains(edge.getColor() )){
        					  nodeTraveledA.add(node);
        					  nodeTraveledCurrent.add(node);
        				  }
        				  else if (startB  && characterStartB.contains(edge.getColor() )) {
        					  nodeTraveledB.add(node);
        					  nodeTraveledCurrent.add(node);
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
        					  siteA =  characterEndA.contains(edge.getColor());
        					  siteB =  characterEndB.contains(edge.getColor());
        				  } 
        				  else {
        					  siteA = characterStartA.contains(edge.getColor());
        					  siteB = characterStartB.contains(edge.getColor());
        				  }
        				  if(siteA) {
        					  nodeTraveledA.add(gOnlyF.getNode(newNode));
        					  nodeTraveledCurrent.add(gOnlyF.getNode(newNode));
        				  }
        				  else if(siteB) {
        					  nodeTraveledB.add(gOnlyF.getNode(newNode));
        					  nodeTraveledCurrent.add(gOnlyF.getNode(newNode));
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
        						  newSiteA = characterStartA.contains(newEdge.getColor());
        						  newSiteB = characterStartB.contains(newEdge.getColor());

        					  } 
        					  else {
        						  newSiteA = characterEndA.contains(newEdge.getColor());
        						  newSiteB = characterEndB.contains(newEdge.getColor());
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
        		  for (Node nodeTraveled:nodeTraveledCurrent){
        			  allNodes.add(nodeTraveled.getId());
        		  }

        		  if (!allNodes.isEmpty()){
        			  components.add(allNodes);
        			  nodeTraveledCurrent.clear();
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
          if(!noFR){// if this is not a fR bond, then ???
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
	//inner class
	public static class BondConnectedSubstitutionParameters2SiteNa implements Parameters {//info
		public final char fRBond, eRBond;
		public final List<Character> characterStartA, characterStartB,characterEndA,characterEndB;
		//constructor 
        public BondConnectedSubstitutionParameters2SiteNa (char fRBond, char eRBond, List<Character> characterStartA, 
        		List<Character> characterStartB,List<Character> characterEndA,List<Character> characterEndB){
        	this.fRBond = fRBond;
        	this.eRBond = eRBond;
	    	this.characterEndA = characterEndA;
	    	this.characterEndB =characterEndB ;
	    	this.characterStartA=characterStartA;
	    	this.characterStartB=characterStartB;
	    	
	    	
	    }
	  }
}
