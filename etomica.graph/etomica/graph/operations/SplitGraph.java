package etomica.graph.operations;

import static etomica.graph.traversal.Traversal.STATUS_ARTICULATION_POINT;
import static etomica.graph.traversal.Traversal.STATUS_START_COMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_COMPONENT;
import static etomica.graph.traversal.Traversal.STATUS_VISITED_NODE;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.property.IsConnected;
import etomica.graph.traversal.DepthFirst;
import etomica.graph.traversal.TraversalVisitor;

/**
 * Split a disconnected graph into multiple graphs.
 * 
 * @author Andrew Schultz
 */
public class SplitGraph {

  protected IsConnected isCon = new IsConnected();
    
  public Set<Graph> apply(Graph g) {
    Set<Graph> result = new HashSet<Graph>();
    if (isCon.check(g)) {
      result.add(g.copy());
      return result;
    }
    List<List<Byte>> components = new ArrayList<List<Byte>>();
    CVisitor v = new CVisitor(components);
    new DepthFirst().traverseAll(g, v);

    for (int iComp = 0; iComp<components.size(); iComp++) {
      List<Byte> nodes = components.get(iComp);
      Graph newG = GraphFactory.createGraph((byte)nodes.size());
      result.add(newG);
      for (byte id = 0; id<nodes.size(); id++) {
        byte nodeID = nodes.get(id);
        newG.getNode(id).setColor(g.getNode(nodeID).getColor());
        newG.getNode(id).setType(g.getNode(nodeID).getType());
        for (byte id2 = (byte)(id+1); id2<nodes.size(); id2++) {
          byte nodeID2 = nodes.get(id2);
          if (g.hasEdge(nodeID, nodeID2)) {
            newG.putEdge(id, id2);
            newG.getEdge(id, id2).setColor(g.getEdge(nodeID, nodeID2).getColor());
          }
        }
      }
      if (iComp == 0) {
        newG.coefficient().multiply(g.coefficient());
      }
    } 
    return result;
  }

  public static class CVisitor implements TraversalVisitor {

    private List<List<Byte>> components;
    private List<Byte> component;
    private boolean isArticulation = false;

    public CVisitor(List<List<Byte>> biComponents) {
      this.components = biComponents;
    }

    public boolean visit(byte nodeID, byte status) {

      // the next node is an articulation point and should not be processed
      if (status == STATUS_START_COMPONENT) {
        component = new ArrayList<Byte>();
        components.add(component);
      }
      else if (status == STATUS_VISITED_COMPONENT) {
        component = null;
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
          component.add(nodeID);
        }
      }
      return true;
    }
  }
}
