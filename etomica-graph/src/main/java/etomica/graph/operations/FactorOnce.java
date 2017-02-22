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
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.traversal.BCVisitor;

/**
 * Factor graphs once.  Each graph with an articulation point is factored into
 * two components.  Articulation points that correspond to flexible molecules
 * are ignored.
 *
 * This operation throws an exception if a graph cannot be factored.
 * 
 * @author Andrew Schultz
 */
public class FactorOnce implements Unary {

  protected HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    
  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
    assert(params instanceof FactorOnceParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      result.addAll(apply(g, (FactorOnceParameters)params));
    }
    return result;
  }
    
  public Set<Graph> apply(Graph g, FactorOnceParameters fop) {
    if (!hap.check(g)) {
      throw new RuntimeException("unfactorable");
    }
    List<List<Byte>> biComponents = BCVisitor.getBiComponents(g);

    List<Integer> myComponents = new ArrayList<Integer>();
    for (int i=0; i<biComponents.size(); i++) {
      if (biComponents.get(i).contains(fop.nodeId)) {
        myComponents.add(i);
      }
    }
    Set<Graph> resultSet = new HashSet<Graph>();
    if (myComponents.size() < 2) {
      resultSet.add(g);
      return resultSet;
    }
    
    List<Integer> iComps = new ArrayList<Integer>();
    List<Integer> jComps = new ArrayList<Integer>();
    for (int k=0; k<(1<<myComponents.size()-1); k++) {
      iComps.clear();
      jComps.clear();
      iComps.add(myComponents.get(0));
      for (int l=1; l<myComponents.size(); l++) {
        if ((k & (1<<(l-1))) > 0) {
          iComps.add(myComponents.get(l));
        }
        else {
          jComps.add(myComponents.get(l));
        }
      }
      if (jComps.isEmpty()) continue;
      resultSet.add(apply(g, iComps, jComps));
      if (!fop.allPermutations) return resultSet;
    }
    if (resultSet.size() > 1) {
      MulScalar mulScalar = new MulScalar();
      MulScalarParameters mspFactor = new MulScalarParameters(1, resultSet.size());
      resultSet = mulScalar.apply(resultSet, mspFactor);
    }
    return resultSet;
  }

  // components in iComps will be combined together
  // components not in iComps will be combined together
  public Graph apply(Graph g, List<Integer> iComps, List<Integer> jComps) {
    if (!hap.check(g)) {
      throw new RuntimeException("unfactorable or disconnected");
    }
    if (iComps.size() == 0 || jComps.size() == 0) {
      throw new RuntimeException("you gotta give me both i and j components");
    }
    List<List<Byte>> biComponents = BCVisitor.getBiComponents(g);

    int lastSize = biComponents.size()+1;
    boolean lastDitch = false;
    
    while (biComponents.size() > 2 && (biComponents.size() < lastSize || !lastDitch)) {
      if (!lastDitch && biComponents.size() == lastSize) {
        lastDitch = true;
      }
      lastSize = biComponents.size();
      for (int i=0; i<biComponents.size()-1; i++) {
        for (int in=0; in<biComponents.get(i).size(); in++) {
          byte iNodeID = biComponents.get(i).get(in);
          for (int j=i+1; j<biComponents.size(); j++) {
            for (byte jNodeID : biComponents.get(j)) {
              if (iNodeID == jNodeID) {
                
                // i and j in iComps => true
                // i and j in jComps => true
                // i iComps or jComps, j not in iComps or jComps => true
                // j iComps or jComps, i not in iComps or jComps => true
                // lastDitch (we've already condensed to iComps and jComps as much as possible), i and j both not in iComps or jComps
                if ((iComps.contains(i) && iComps.contains(j)) || (jComps.contains(i) && jComps.contains(j)) ||
                    ((!iComps.contains(i) && !jComps.contains(i)) && (iComps.contains(j) || jComps.contains(j))) || 
                    ((!iComps.contains(j) && !jComps.contains(j)) && (iComps.contains(i) || jComps.contains(i))) ||
                    (lastDitch && !iComps.contains(i) && !jComps.contains(i) && !iComps.contains(j) && !jComps.contains(j))) {
                  // we want to keep these components together
                  for (byte jNodeID2 : biComponents.get(j)) {
                    if (jNodeID2 != iNodeID) {
                      biComponents.get(i).add(jNodeID2);
                    }
                  }
                  biComponents.remove(biComponents.get(j));
                  for (int ic = 0; ic<iComps.size(); ic++) {
                    if (iComps.get(ic) == j) {
                      if (iComps.contains(i)) {
                        iComps.remove(ic);
                        ic--;
                      }
                      else {
                        iComps.set(ic, i);
                      }
                    }
                    else if (iComps.get(ic) > j) {
                      iComps.set(ic, iComps.get(ic)-1);
                    }
                  }
                  for (int jc = 0; jc<jComps.size(); jc++) {
                    if (jComps.get(jc) == j) {
                      if (jComps.contains(i)) {
                        jComps.remove(jc);
                        jc--;
                      }
                      else {
                        jComps.set(jc, i);
                      }
                    }
                    else if (jComps.get(jc) > j) {
                      jComps.set(jc, jComps.get(jc)-1);
                    }
                  }
                  j--;
                }
                break;
              }
            }
          }
        }
      }
    }
    List<Set<Byte>> newRootNodes = new ArrayList<Set<Byte>>();
    for (int i=0; i<biComponents.size(); i++) {
      newRootNodes.add(new HashSet<Byte>());
    }
    for (int i=0; i<biComponents.size()-1; i++) {
      for (int in=0; in<biComponents.get(i).size(); in++) {
        byte iNodeID = biComponents.get(i).get(in);
        for (int j=i+1; j<biComponents.size(); j++) {
          for (byte jNodeID : biComponents.get(j)) {
            if (iNodeID == jNodeID) {
              // we'll separate these components
              newRootNodes.get(j).add(iNodeID);
              break;
            }
          }
        }
      }
    }

    byte newNodeCount = 0;
    HashMap<Byte,Byte>[] byteMaps = new HashMap[biComponents.size()];
    byte newRootNodeId = g.nodeCount();
    for (int i=0; i<biComponents.size(); i++) {
      newNodeCount += biComponents.get(i).size();
      byteMaps[i] = new HashMap<Byte,Byte>();
      for (byte id : biComponents.get(i)) {
        if (newRootNodes.get(i).contains(id)) {
          byteMaps[i].put(id,newRootNodeId);
          newRootNodeId++;
        }
        else {
          byteMaps[i].put(id,id);
        }
      }
    }
    Graph result = GraphFactory.createGraph(newNodeCount);
    for (int i= 0; i<biComponents.size(); i++) {
      for (byte nodeID : biComponents.get(i)) {
        byte newNodeId = byteMaps[i].get(nodeID);
        result.getNode(newNodeId).setColor(g.getNode(nodeID).getColor());
        result.getNode(newNodeId).setType(g.getNode(nodeID).getType());
        if (newRootNodes.get(i).contains(nodeID)) {
          result.getNode(newNodeId).setType(TYPE_NODE_ROOT);
        }
        for (byte nodeID2 : biComponents.get(i)) {
          if (nodeID2 <= nodeID || !g.hasEdge(nodeID, nodeID2)) continue;
          byte newNodeId2 = byteMaps[i].get(nodeID2);
          result.putEdge(newNodeId, newNodeId2);
          result.getEdge(newNodeId, newNodeId2).setColor(g.getEdge(nodeID, nodeID2).getColor());
        }
      }
    } 

    result.coefficient().multiply(g.coefficient());

    result.setNumFactors(g.factors().length);
    result.addFactors(g.factors());
    result.createReverseEdges();

    return result;
  }

  public static class FactorOnceParameters implements Parameters {
    public final byte nodeId;
    public final boolean allPermutations;
    public FactorOnceParameters(byte nodeId, boolean allPermutations) {
      this.nodeId = nodeId;
      this.allPermutations = allPermutations;
    }
  }
}
