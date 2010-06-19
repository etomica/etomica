package etomica.graph.operations;

import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Node;

/**
 * Perform graph multiplication for flexible molecules.  Where possible, the
 * result is constructed by superimposing a from one diagram on top
 * of a node from another diagram.  The superimposing is possible when the
 * following criteria are met:
 * 1. the nodes being superimposed must be the same color.
 * 2. At least one of the nodes must be a root node.
 * 3a. The node color must correspond to a rigid molecule.
 * or
 * 3b. One of the nodes must have no edges.
 * 
 * When not possible to superimpose, the result of multiplication
 * for each pair of diagrams is a diagram composed by both original diagrams
 * without any connection between them (the new diagram is disconnected).
 *
 * @author Andrew Schultz
 */
public class MulFlexible implements Binary {

  public Set<Graph> apply(Set<Graph> argument, Set<Graph> argument2, Parameters params) {

    assert (params instanceof MulFlexibleParameters);
    Set<Graph> result = new HashSet<Graph>();
    for (Graph g : argument) {
      for (Graph g2 : argument2) {
        Graph newGraph = apply(g, g2, (MulFlexibleParameters)params);
        result.add(newGraph);
        apply(g,g2, (MulFlexibleParameters)params);
      }
    }
    return result;
  }

  public Graph apply(Graph g1, Graph g2, MulFlexibleParameters params) {
    Graph result;
    Node myNode1 = null;
    Node myNode2 = null;
    for (Node node1 : g1.nodes()) {
      boolean success = true;
      boolean flexColor = false;
      for (int i=0; i<params.flexColors.length; i++) {
        if (params.flexColors[i] == node1.getColor()) {
          flexColor = true;
        }
      }
      if (flexColor) {
        for (Node node2 : g1.nodes()) {
          if (node1 != node2 && g1.hasEdge(node1.getId(), node2.getId())) {
            success = false;
            break;
          }
        }
      }
      if (success) {
        // check if the other graph has a node of the same color
        for (Node node2 : g2.nodes()) {
          if (node2.getColor() == node1.getColor() && (node1.getType() == TYPE_NODE_ROOT || node2.getType() == TYPE_NODE_ROOT)) {
            myNode1 = node1;
            myNode2 = node2;
            break;
          }
        }
        if (myNode1 != null) break;
      }
    }
    if (myNode1 == null) {
      for (Node node1 : g2.nodes()) {
        if (node1.getType() != TYPE_NODE_ROOT) continue;
        boolean success = true;
        boolean flexColor = false;
        for (int i=0; i<params.flexColors.length; i++) {
          if (params.flexColors[i] == node1.getColor()) {
            flexColor = true;
          }
        }
        if (flexColor) {
          for (Node node2 : g2.nodes()) {
            if (node1 != node2 && g2.hasEdge(node1.getId(), node2.getId())) {
              success = false;
              break;
            }
          }
        }
        if (success) {
          // check if the other graph has a node of the same color
          for (Node node2 : g1.nodes()) {
            if (node2.getColor() == node1.getColor() && (node1.getType() == TYPE_NODE_ROOT || node2.getType() == TYPE_NODE_ROOT)) {
              Graph tmp = g1;
              g1 = g2;
              g2 = tmp;
              myNode1 = node1;
              myNode2 = node2;
              break;
            }
          }
          if (myNode1 != null) break;
        }
      }
    }

    byte nodesOffset = g1.nodeCount();
    byte newNodeCount = (byte)(nodesOffset + g2.nodeCount());
    if (myNode1 != null) {
      newNodeCount--;
    }
    result = GraphFactory.createGraph(newNodeCount);
    // add edges from g1
    for (Node node1 : g1.nodes()) {
      if (node1 != myNode1) {
        result.getNode(node1.getId()).setType(node1.getType());
      }
      else if (myNode1.getType() == TYPE_NODE_ROOT && myNode2.getType() == TYPE_NODE_ROOT) {
        result.getNode(node1.getId()).setType(TYPE_NODE_ROOT);
      }

      result.getNode(node1.getId()).setColor(node1.getColor());
      for (Node node2 : g1.nodes()) {
        if (node2.getId() <= node1.getId() || !g1.hasEdge(node1.getId(), node2.getId())) continue;
        result.putEdge(node1.getId(), node2.getId());
        result.getEdge(node1.getId(), node2.getId()).setColor(g1.getEdge(node1.getId(), node2.getId()).getColor());
      }
    }
    // now add edges from g2
    for (Node node1 : g2.nodes()) {
      byte newNodeId = (byte)(node1.getId()+nodesOffset);
      if (node1 != myNode2) {
        if (myNode1 != null && node1.getId() > myNode2.getId()) {
          newNodeId--;
        }
        result.getNode(newNodeId).setType(node1.getType());
        result.getNode(newNodeId).setColor(node1.getColor());
      }
      else if (myNode1 != null) {
        // don't need to set type or color
        newNodeId = myNode1.getId();
      }
      for (Node node2 : g2.nodes()) {
        if (node2.getId() <= node1.getId() || !g2.hasEdge(node1.getId(), node2.getId())) continue;
        byte newNode2Id = (byte)(node2.getId()+nodesOffset);
        if (node2 == myNode2) {
          newNode2Id = myNode1.getId();
        }
        else if (myNode1 != null && node2.getId() > myNode2.getId()) {
          newNode2Id--;
        }
        result.putEdge(newNodeId, newNode2Id);
        result.getEdge(newNodeId, newNode2Id).setColor(g2.getEdge(node1.getId(), node2.getId()).getColor());
      }
    }
    result.coefficient().multiply(g1.coefficient());
    result.coefficient().multiply(g2.coefficient());

    result.setNumFactors(g1.factors().length);
    result.addFactors(g1.factors());
    result.addFactors(g2.factors());
    return result;
  }

  public static class MulFlexibleParameters implements Parameters {
    public final char[] flexColors;
    public MulFlexibleParameters(char[] flexColors) {
      this.flexColors = flexColors;
    }
  }
}
