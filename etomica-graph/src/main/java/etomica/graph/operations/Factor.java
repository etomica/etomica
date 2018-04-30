/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Metadata;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.traversal.BCVisitor;

/**
 * Factor graphs.  Each graph with an articulation point is factored into as
 * many biconnected components as possible.  Biconnected components connected
 * by flexible molecule nodes are not separated.
 * 
 * @author Andrew Schultz
 */
public class Factor implements Unary {

  protected HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert(params instanceof MulFlexibleParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.add(apply(g, (MulFlexibleParameters)params));
    }
    return result;
  }

  @SuppressWarnings("unchecked")
  public Graph apply(Graph g, MulFlexibleParameters mfp) {
    if (!hap.check(g)) {
      return g.copy();
    }
    List<List<Byte>> biComponents = BCVisitor.getBiComponents(g);
    List<List<Byte>> newRootNodes = new ArrayList<List<Byte>>();
    for (int i=0; i<biComponents.size(); i++) {
      newRootNodes.add(new ArrayList<Byte>());
    }

    for (int i=0; i<biComponents.size()-1; i++) {
      List<Byte> iBiComponent = biComponents.get(i);
      for (int k=0; k<iBiComponent.size(); k++) {
        byte iNodeID = iBiComponent.get(k);
        for (int j=i+1; j<biComponents.size(); j++) {
          for (byte jNodeID : biComponents.get(j)) {
            if (iNodeID == jNodeID) {
              boolean flexColor = false;
              for (int ic=0; ic<mfp.flexColors.length; ic++) {
                if (mfp.flexColors[ic] == g.getNode(iNodeID).getColor()) {
                  flexColor = true;
                }
              }
              if (flexColor) {
                // flexible-molecule color, can't break it up.
                for (byte jNodeID2 : biComponents.get(j)) {
                  if (jNodeID2 != iNodeID) {
                    biComponents.get(i).add(jNodeID2);
                  }
                }
                newRootNodes.get(i).addAll(newRootNodes.get(j));
                newRootNodes.remove(newRootNodes.get(j));
                biComponents.remove(biComponents.get(j));
                j--;
                break;
              }
              if (g.getNode(iNodeID).getType() != TYPE_NODE_ROOT) {
                newRootNodes.get(j).add(iNodeID); 
              }
            }
          }
        }
      }
    }
    
    if (biComponents.size() == 1) return g.copy();

    byte newNodeCount = 0;
    HashMap<Byte,Byte>[] byteMaps = new HashMap[biComponents.size()];
    byte myNode = -1;
    for (int i=0; i<biComponents.size(); i++) {
      newNodeCount += biComponents.get(i).size();
      byteMaps[i] = new HashMap<Byte,Byte>();
      for (byte id : biComponents.get(i)) {
        myNode++;
        byteMaps[i].put(id,myNode);
      }
    }
    Graph result = GraphFactory.createGraph(newNodeCount);
    for (int iComp = 0; iComp<biComponents.size(); iComp++) {
      for (byte nodeID : biComponents.get(iComp)) {
        byte newNodeId = byteMaps[iComp].get(nodeID);
        result.getNode(newNodeId).setColor(g.getNode(nodeID).getColor());
        result.getNode(newNodeId).setType(g.getNode(nodeID).getType());
        if (newRootNodes.get(iComp).contains(nodeID)) {
          result.getNode(newNodeId).setType(TYPE_NODE_ROOT);
        }
        for (byte nodeID2 : biComponents.get(iComp)) {
          if (nodeID2 <= nodeID || !g.hasEdge(nodeID, nodeID2)) continue;
          byte newNodeId2 = byteMaps[iComp].get(nodeID2);
          result.putEdge(newNodeId, newNodeId2);
          if (newNodeId > newNodeId2) {
            // node order is reversed, but we can't assign color for a reverse edge (if that is in effect).
            // retrieve the new forward edge and make it equal to the color of the reverse edge from the original graph.
            result.getEdge(newNodeId2, newNodeId).setColor(g.getEdge(nodeID2, nodeID).getColor());//calling setColor on Edge
          }
          else {
            result.getEdge(newNodeId, newNodeId2).setColor(g.getEdge(nodeID, nodeID2).getColor());//calling setColor on Edge
          }
        }
      }
    } 

    if (!MetadataImpl.rootPointsSpecial) {
      // we might have started with a graph with a single point and then created a root point (through factorization)
      // try now to condense such a graph by superimposing single point on the root node. 
      while (true) {
        byte singlePoint = -1;
        byte rootNode = -1;
        for (byte i=0; i<result.nodeCount(); i++) {
          if (result.getOutDegree(i) == 0) {
            for (byte j=0; j<result.nodeCount(); j++) {
              if (result.getNode(j).getType() == Metadata.TYPE_NODE_ROOT && result.getNode(j).getColor() == result.getNode(i).getColor()) {
                singlePoint = i;
                rootNode = j;
                break;
              }
            }
            if (rootNode > -1) break;
          }
        }
        if (singlePoint == -1) {
          break;
        }
        Graph newResult = GraphFactory.createGraph((byte)(result.nodeCount()-1));
        byte j = 0;
        for (byte i=0; i<result.nodeCount(); i++) {
          if (i==singlePoint) continue;
          newResult.getNode(j).setColor(result.getNode(i).getColor());
          if (i==rootNode) {
            newResult.getNode(j).setType(result.getNode(singlePoint).getType());
          }
          else {
            newResult.getNode(j).setType(result.getNode(i).getType());
          }
          
          for (byte k=(byte)(i+1); k<result.nodeCount(); k++) {
            if (k == singlePoint || !result.hasEdge(i,k)) continue;
            byte kNew = k;
            if (k>singlePoint) kNew--;
            newResult.putEdge(j,kNew);
            newResult.getEdge(j,kNew).setColor(result.getEdge(i,k).getColor());
          }
          j++;
        }
        result = newResult;
      }
    }

    result.coefficient().multiply(g.coefficient());

    result.setNumFactors(g.factors().length);
    result.addFactors(g.factors());
    result.createReverseEdges();

    return result;
  }
}
