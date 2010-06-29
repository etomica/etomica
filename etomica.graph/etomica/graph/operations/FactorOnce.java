package etomica.graph.operations;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;
import static etomica.graph.traversal.Traversal.STATUS_ARTICULATION_POINT;
import static etomica.graph.traversal.Traversal.STATUS_START_BICOMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_BICOMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_NODE;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.traversal.Biconnected;
import etomica.graph.traversal.TraversalVisitor;

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
      result.add(apply(g, (FactorOnceParameters)params));
    }
    return result;
  }
    
  public Graph apply(Graph g, FactorOnceParameters mfp) {
    if (!hap.check(g)) {
      throw new RuntimeException("unfactorable");
    }
    List<List<Byte>> biComponents = new ArrayList<List<Byte>>();
    BCVisitor v = new BCVisitor(biComponents);
    new Biconnected().traverseAll(g, v);
    List<Set<Byte>> newRootNodes = new ArrayList<Set<Byte>>();
    for (int i=0; i<biComponents.size(); i++) {
      newRootNodes.add(new HashSet<Byte>());
    }

    for (int i=0; i<biComponents.size()-1; i++) {
      for (int in=0; in<biComponents.get(i).size(); in++) {
        byte iNodeID = biComponents.get(i).get(in);
        if (iNodeID != mfp.nodeId) continue;
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
              else {
                if (g.getNode(iNodeID).getType() != TYPE_NODE_ROOT) {
                  newRootNodes.get(j).add(iNodeID); 
                }
              }
            }
          }
        }
      }
    }

    if (biComponents.size() == 1) throw new RuntimeException("unfactorable");
    if (biComponents.size() > 2) {
      for (int i=0; i<biComponents.size()-1; i++) {
        for (int in=0; in<biComponents.get(i).size(); in++) {
          byte iNodeID = biComponents.get(i).get(in);
          for (int j=i+1; biComponents.size() > 2 && j<biComponents.size(); j++) {
            for (byte jNodeID : biComponents.get(j)) {
              if (iNodeID == jNodeID) {
                // we only want 2 components, so recombine i and j, which share iNodeID
                for (byte jNodeID2 : biComponents.get(j)) {
                  if (jNodeID2 != iNodeID) {
                    biComponents.get(i).add(jNodeID2);
                  }
                }
                newRootNodes.get(j).remove(jNodeID);
                newRootNodes.get(i).addAll(newRootNodes.get(j));
                newRootNodes.remove(newRootNodes.get(j));
                biComponents.remove(biComponents.get(j));
                j--;
                if (biComponents.size() == 2) {
                  break;
                }
              }
            }
          }
        }
      }
    }
    
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
          result.getEdge(newNodeId, newNodeId2).setColor(g.getEdge(nodeID, nodeID2).getColor());
        }
      }
    } 

    result.coefficient().multiply(g.coefficient());

    result.setNumFactors(g.factors().length);
    result.addFactors(g.factors());

    return result;
  }

  public static class BCVisitor implements TraversalVisitor {

    private List<List<Byte>> biComponents;
    private List<Byte> biComponent;
    private boolean isArticulation = false;

    public BCVisitor(List<List<Byte>> biComponents) {
      this.biComponents = biComponents;
    }

    public boolean visit(byte nodeID, byte status) {

      // the next node is an articulation point and should not be processed
      if (status == STATUS_START_BICOMPONENT) {
        biComponent = new ArrayList<Byte>();
        biComponents.add(biComponent);
      }
      else if (status == STATUS_VISITED_BICOMPONENT) {
        biComponent = null;
      }
      else if  (status == STATUS_ARTICULATION_POINT) {
        isArticulation = true;
      }
      // visiting a node in the current biconnected component
      else if (status == STATUS_VISITED_NODE) {
        // if it is an articulation point, ignore it
        if (isArticulation) {
          isArticulation = false;
        }
        else {
          biComponent.add(nodeID);
        }
      }
      return true;
    }
  }

  public static class FactorOnceParameters implements Parameters {
    public final byte nodeId;
    public final char[] flexColors;
    public FactorOnceParameters(byte nodeId, char[] flexColors) {
      this.nodeId = nodeId;
      this.flexColors = flexColors;
    }
  }
}
