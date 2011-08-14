package etomica.graph.operations;

import static etomica.graph.model.Metadata.TYPE_NODE_FIELD;
import static etomica.graph.model.Metadata.TYPE_NODE_ROOT;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.Mul.MulParameters;

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
        if (newGraph != null) {
          result.add(newGraph);
        }
      }
    }
    return result;
  }

  public Graph apply(Graph g1, Graph g2, MulFlexibleParameters params) {

    int numNodes = 0;
    for (Node node : g1.nodes()) {
      if (node.getType() == TYPE_NODE_FIELD) {
        numNodes++;
      }
    }
    for (Node node : g2.nodes()) {
      if (node.getType() == TYPE_NODE_FIELD) {
        numNodes++;
      }
    }
    if (numNodes > params.nFieldPoints) {
      return null;
    }

    
    Graph result;
    Node myNode1 = null;
    Node myNode2 = null;
    // look for root nodes.  we might find superimposable nodes, neither of
    // which are root nodes.  if possible, just swap out the root node.
    byte g1RootNode = -1;
    byte g2RootNode = -1;
    boolean doubleRoot = false;
    for (Node node1 : g1.nodes()) {
      if (node1.getType() == Metadata.TYPE_NODE_ROOT) {
        g1RootNode = node1.getId();
        break;
      }
    }
    if (g1RootNode == -1) {
      for (Node node2 : g2.nodes()) {
        if (node2.getType() == Metadata.TYPE_NODE_ROOT) {
          g2RootNode = node2.getId();
          break;
        }
      }
    }
    for (Node node1 : g1.nodes()) {
      if (params.skipSuperimpose) break;
      if (params.node1ID > -1) {
        // only check the specified node
        node1 = g1.getNode(params.node1ID);
      }
      boolean flex1Unhappy = false;
      boolean flexColor = false;
      for (int i=0; i<params.flexColors.length; i++) {
        if (params.flexColors[i] == node1.getColor()) {
          flexColor = true;
        }
      }
      if (flexColor) {
        // we want node1 to be unbonded
        flex1Unhappy = g1.getOutDegree(node1.getId()) > 0;
      }
      // check if the other graph has a node of the same color
      for (Node node2 : g2.nodes()) {
        if (params.node2ID > -1) {
          // only check the specified node
          node2 = g2.getNode(params.node2ID);
        }
        if (node1.getColor() == node2.getColor()) {
          boolean success = true;
          if (flex1Unhappy) {
            // node1 was bonded.  node2 must unbonded
            success = g2.getOutDegree(node2.getId()) == 0;
          }
          if (success) {
            if (node1.getType() == TYPE_NODE_ROOT || node2.getType() == TYPE_NODE_ROOT) {
              // node1 and node2 are suitable for superimposing
              myNode1 = node1;
              myNode2 = node2;
              g1RootNode = -1;
              g2RootNode = -1;
              doubleRoot = node1.getType() == node2.getType();
              if (doubleRoot) break;
              // if both nodes are not root, then keep looking
            }
            else if (!MetadataImpl.rootPointsSpecial && (g1RootNode > -1 || g2RootNode > -1) && myNode1 == null) {
              // node1 and node2 are suitable for superimposing
              // but rootiness needs to come from somewhere else
              // this actually short-circuits the brute force search for superimposable nodes
              //   which is fine
              myNode1 = node1;
              myNode2 = node2;
              break;
            }
          }
        }
        if (params.node2ID > -1) {
          // unable to superimpose specified node2, so bail
          break;
        }
      }
      if (params.node1ID > -1) {
        // unable to superimpose specified node1, so bail
        break;
      }
      if (myNode1 != null && doubleRoot) break;
    }

    byte nodesOffset = g1.nodeCount();
    byte newNodeCount = (byte)(nodesOffset + g2.nodeCount());
    if (myNode1 != null) {
      newNodeCount--;
    }
    else {
      // we didn't find color-happy nodes so don't go moving root nodes around
      g1RootNode = -1;
      g2RootNode = -1;
    }
    result = GraphFactory.createGraph(newNodeCount);
    // add edges from g1
    for (Node node1 : g1.nodes()) {
      if (node1 != myNode1 && node1.getId() != g1RootNode) {
        result.getNode(node1.getId()).setType(node1.getType());
      }
      else if (myNode1 != null && (myNode1.getType() == TYPE_NODE_ROOT && myNode2.getType() == TYPE_NODE_ROOT)) {
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
        if (node1.getId() != g2RootNode) {
          result.getNode(newNodeId).setType(node1.getType());
        }
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

  public static class MulFlexibleParameters extends MulParameters {
    public final char[] flexColors;
    public final byte node1ID, node2ID;
    public final boolean skipSuperimpose;
    
    /**
     * Parameters that specify that the resulting product should include graphs
     * containing up to nFieldNodes field nodes and that nodes should never be
     * superimposed.
     */
    public MulFlexibleParameters(byte nFieldNodes) {
      super(nFieldNodes);
      skipSuperimpose = true;
      flexColors = new char[0];
      node1ID = node2ID = (byte)-1;
    }

    /**
     * Parameters that specify that the resulting product should include graphs
     * containing up to nFieldNodes field nodes and that nodes of a color
     * included in flexColors never be superimposed unless one of the nodes is
     * unbonded.
     */
    public MulFlexibleParameters(char[] flexColors, byte nFieldNodes) {
      this(flexColors, nFieldNodes, (byte)-1, (byte)-1);
    }

    /**
     * Parameters that specify that the resulting product should include graphs
     * containing up to nFieldNodes field nodes and that node1ID from the first
     * graph and node2ID from the second graph should be superimposed (if
     * possible).
     */
    public MulFlexibleParameters(char[] flexColors, byte nFieldNodes, byte node1ID, byte node2ID) {
      super(nFieldNodes);
      skipSuperimpose = false;
      this.flexColors = flexColors;
      this.node1ID = node1ID;
      this.node2ID = node2ID;
    }
  }
}
