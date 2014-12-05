/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.HashSet;
import java.util.Set;

import etomica.virial.cluster2.graph.impl.BCTraversal;
import etomica.virial.cluster2.graph.impl.BFTraversal;
import etomica.virial.cluster2.graph.impl.DFTraversal;
import etomica.virial.cluster2.graph.impl.NPTraversal;

/**
 * This class supports a number of graph algorithms through a simple delegation mechanism.
 * Improved algorithm implementations can be swapped by assigning algorithm instances to
 * the appropriate algorithm static fields. Most of the default implementations of these
 * algorithms are provided in this file as private classes.
 *
 * @author Demian Lessa
 */
public class Algorithms {

  public static GraphProperty algoArticulationPair = new CheckArticulatonPair();

  public static GraphProperty algoArticulationPoint = new CheckArticulationPoint();

  public static GraphProperty algoBiconnected = new CheckBiconnected();

  public static GraphProperty algoConnected = new CheckConnected();

  public static GraphProperty algoNodalPoint = new CheckNodalPoint();

  public static GraphTraversal algoBFT = new BFTraversal();

  public static GraphTraversal algoDFT = new DFTraversal();

  public static GraphTraversal algoNPT = new NPTraversal();

  public static GraphTraversal algoBCT = new BCTraversal();

  public static boolean hasArticulationPair(final Nodes nodes, final Edges edges) {

    return algoArticulationPair.check(nodes, edges);
  }

  public static boolean hasArticulationPoint(final Nodes nodes, final Edges edges) {

    return algoArticulationPoint.check(nodes, edges);
  }

  public static boolean hasNodalPoint(final Nodes nodes, final Edges edges) {

    return algoNodalPoint.check(nodes, edges);
  }

  public static boolean isBiconnected(final Nodes nodes, final Edges edges) {

    return algoBiconnected.check(nodes, edges);
  }

  public static boolean isConnected(final Nodes nodes, final Edges edges) {

    return algoConnected.check(nodes, edges);
  }

  public static boolean traverseBF(int nodeID, final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    return algoBFT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseBF(final Nodes nodes, final Edges edges, final NodesVisitor visitor) {

    algoBFT.traverseAll(nodes, edges, visitor);
  }

  public static boolean traverseDF(int nodeID, final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    return algoDFT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseDF(final Nodes nodes, final Edges edges, final NodesVisitor visitor) {

    algoDFT.traverseAll(nodes, edges, visitor);
  }

  public static boolean traverseBC(int nodeID, final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    return algoBCT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseBC(final Nodes nodes, final Edges edges, final NodesVisitor visitor) {

    algoBCT.traverseAll(nodes, edges, visitor);
  }
}

/**
 * Articulation Point algorithm.
 *
 * @author Demian Lessa
 */
class CheckArticulatonPair implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    // TODO Auto-generated method stub
    return false;
  }
}

/**
 * DEF 1. (K-Connected Graph) A graph G is said to be K-connected if there does not exist
 * a set of K-1 nodes whose removal disconnects. Null and singleton graphs are K-connected
 * for all K >= 1.
 *
 * Ref: Wolfram MathWorld, http://mathworld.wolfram.com/k-ConnectedGraph.html
 *
 *
 * PROP 2. (Virial Nodal Point) A graph G has a virial nodal point if: a) G is connected,
 * b) there exists a field node that is an articulation point, and c) all root nodes of G
 * «are not» in the same nodal bicomponent.
 *
 * DEF 2.1 (Nodal Bicomponent) A nodal bicomponent of a graph G is: a) a bicomponent BC of
 * G if all articulation points in BC are field nodes, or b) a bicomponent BC of G merged
 * with all nodal components of G that share a root articulation point with BC.
 *
 *
 * PROP 3. (Virial Articulation Point) A graph G has a virial articulation point if: a) G
 * is connected, b) there exists articulation points in G, and c) some leaf bicomponent of
 * G contains no root node, or its root node is also its articulation point.
 *
 * DEF 3.1 (Leaf Bicomponent) A bicomponent of a grap G is a leaf bicomponent if it only
 * contains one articulation point.
 *
 *
 * @author Demian Lessa
 */
class CheckConnected implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    // by definition, a null graph and a singleton graph are connected
    if (nodes.count() <= 1) {
      return true;
    }
    // invariant: a connected graph has at least N-1 edges
    if (edges.count() < (nodes.count() - 1)) {
      return false;
    }
    // invoke a DF traversal starting at the first node (nodeID = 0) and
    // return true IFF all nodes in the graph are traversed
    return Algorithms.traverseDF(0, nodes, edges, null);
  }
}

/**
 * A graph is biconnected if, upon removal of any node in the graph, the resulting graph
 * is connected. Based on this definition, null and singleton graphs are both biconnected.
 *
 * @author Demian Lessa
 */
class CheckBiconnected implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    // by definition, a null graph is not biconnected
    if (nodes.count() == 0) {
      return false;
    }
    // by definition, a singleton graph and a two node
    // graph are biconnected
    if (nodes.count() == 1) {
      return true;
    }
    // invariant: a biconnected graph has at least N edges
    if (edges.count() < (nodes.count() - 1)) {
      return false;
    }
    // invoke a BC traversal starting at the first node (nodeID = 0) and
    // return true IFF all nodes in the graph are traversed
    return Algorithms.traverseBC(0, nodes, edges, null);
  }
}

/**
 * Nodal Point algorithm.
 *
 * @author Demian Lessa
 */
class CheckNodalPoint implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    Set<Integer> rootNodes = new HashSet<Integer>();
    for (int nodeID = 0; nodeID < nodes.count(); nodeID++) {
      if (nodes.getAttributes(nodeID).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        rootNodes.add(nodeID);
      }
    }
    // a graph having no field node has no nodal point
    if (nodes.count() == rootNodes.size()) {
      return false;
    }
    // a graph having exactly one root node has no nodal point
    if (rootNodes.size() == 1) {
      return false;
    }
    // a graph having a single field node has a nodal point only
    // if the field node is not connected to all root nodes
    if (nodes.count() - rootNodes.size() == 1) {
      return edges.getOutDegree(nodes.count() - 1) != rootNodes.size();
    }
    // invoke a NP traversal starting and return true IFF the visitor
    // detected that the graph has a nodal point
    NPVisitor v = new NPVisitor(nodes, rootNodes);
    Algorithms.algoNPT.traverseAll(nodes, edges, v);
    return v.hasNodalPoint();
  }
}

/**
 * Visits all biconnected components of a graph whose traversals start with a field node.
 * A graph is nodal if all its root nodes are in the same biconnected component.
 *
 */
class NPVisitor implements NodesVisitor {

  private Nodes nodes;
  private int rootNodesVisited = 0;
  private boolean notNodal = false;
  private boolean isArticulation = false;
  private Set<Integer> rootNodes = new HashSet<Integer>();

  public NPVisitor(Nodes nodes, Set<Integer> rootNodes) {

    this.nodes = nodes;
    this.rootNodes = rootNodes;
  }

  public boolean visit(int nodeID) {

    if (nodeID == GraphTraversal.VISITED_BICOMPONENT) {
      notNodal = notNodal || (rootNodesVisited == rootNodes.size());
      rootNodesVisited = 0;
    }
    // the next node is an articulation point and should not be processed
    else if (nodeID == GraphTraversal.ARTICULATION_POINT) {
      isArticulation = true;
    }
    // visiting a node in the current biconnected component
    else if (nodeID > GraphTraversal.VISITED_NONE) {
      // if it is an articulation point, ignore it
      if (isArticulation) {
        isArticulation = !isArticulation;
      }
      else if (nodes.getAttributes(nodeID).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        rootNodesVisited++;
      }
    }
    return true;
  }

  public boolean hasNodalPoint() {

    return (!notNodal);
  }
}

/**
 * Articulation Point algorithm.
 *
 * @author Demian Lessa
 */
class CheckArticulationPoint implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    Set<Integer> rootNodes = new HashSet<Integer>();
    for (int nodeID = 0; nodeID < nodes.count(); nodeID++) {
      if (nodes.getAttributes(nodeID).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        rootNodes.add(nodeID);
      }
    }
    // invoke a BC traversal starting and return true IFF the visitor
    // detected that the graph has an articulation point
    APVisitor v = new APVisitor(nodes, rootNodes);
    System.out.println(edges);
    Algorithms.algoBCT.traverseAll(nodes, edges, v);
    return v.hasArticulationPoint();
  }
}

class APVisitor implements NodesVisitor {

  private Nodes nodes;
  private int rootNodesVisited = 0;
  private boolean isArticulated = false;
  private boolean isArticulation = false;
  private Set<Integer> rootNodes = new HashSet<Integer>();
  private boolean isConnected = false;
  private Set<Integer> visitedNodes = new HashSet<Integer>();

  public APVisitor(Nodes nodes, Set<Integer> rootNodes) {

    this.nodes = nodes;
    this.rootNodes = rootNodes;
  }

  public boolean visit(int nodeID) {

    if (nodeID == GraphTraversal.VISITED_BICOMPONENT) {
      isArticulated = isArticulated || (rootNodesVisited == 0);
      rootNodesVisited = 0;
    }
    else if (nodeID == GraphTraversal.VISITED_COMPONENT) {
      isConnected = isConnected || (visitedNodes.size() == nodes.count());
      visitedNodes.clear();
    }
    // the next node is an articulation point and should not be processed
    else if (nodeID == GraphTraversal.ARTICULATION_POINT) {
      isArticulation = true;
    }
    // visiting a node in the current biconnected component
    else if (nodeID > GraphTraversal.VISITED_NONE) {
      // if it is an articulation point, ignore it
      if (isArticulation) {
        isArticulation = !isArticulation;
      }
      else {
        visitedNodes.add(nodeID);
        if (nodes.getAttributes(nodeID).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
          rootNodesVisited++;
        }
      }
    }
    return true;
  }

  public boolean hasArticulationPoint() {

    // // this graph is not connected
    // if (!isConnected) {
    // return true;
    // }
    // // this graph has no (graph-theoretical) articulation points
    // if ((!isArticulated) || (rootNodes.size() == 0)) {
    // return false;
    // }
    // // the root node of this graph is an articulation point
    // if (rootNodes.size() == 1) {
    // return articulationPoints.containsAll(rootNodes);
    // }
    // // some bicomponent has no root node, or has a single root node that is an
    // articulation point
    // for (Set<Integer> bc : bicomponents) {
    // int rootPoints = 0;
    // int singletonRootNode = -1;
    // for (Integer node : bc) {
    // if (rootNodes.contains(node)) {
    // rootPoints++;
    // singletonRootNode = node;
    // }
    // }
    // // no root point is found in this bicomponent
    // if (rootPoints == 0) {
    // return true;
    // }
    // // the singleton root point in this bicomponent is an articulation point
    // else if ((rootPoints == 1) && (articulationPoints.contains(singletonRootNode))) {
    // return true;
    // }
    // }
    // if we have not found bicomponents matching the articulation point
    // criteria, it is because this graph has no articulation points
    return false;
  }
}

// TODO: respect the return value of the traversal visitors in order
// to enable early breaking from the traversal.
//
//
// TODO: implement a "nested traversal" of every connected component
// induced by removing an articulation point.
//
// 1) Run a BC traversal to compute the articulation points [O(V+E)]
// 2) If the graph is connected, then perform the nested traversal
// 3) For every articulation point P [there are at most V-2 such points]
// While some neighbor V of P has not been visited
// Perform a DF traversal of G starting from V
// If the property has been decided, break

