package etomica.virial.cluster2.graph;

import etomica.virial.cluster2.graph.impl.BCTraversal;
import etomica.virial.cluster2.graph.impl.BFTraversal;
import etomica.virial.cluster2.graph.impl.DFTraversal;

/**
 * This class supports a number of graph algorithms through a simple delegation
 * mechanism. Improved algorithm implementations can be swapped by assigning
 * algorithm instances to the appropriate algorithm static fields. Most of the
 * default implementations of these algorithms are provided in this file as
 * private classes.
 * 
 * @author Demian Lessa
 */
public class Algorithms {

  public static GraphProperty algoArticulationPair = new CheckArticulatonPair();
  public static GraphProperty algoArticulationPoint = new CheckBiconnected();
  public static GraphProperty algoBiconnected = new CheckBiconnected();
  public static GraphProperty algoConnected = new CheckConnected();
  public static GraphProperty algoNodalPoint = new CheckNodalPoint();
  public static GraphPairProperty algoIsomorphic = new CheckIsomorphic();
  public static GraphTraversal algoBFT = new BFTraversal();
  public static GraphTraversal algoDFT = new DFTraversal();
  public static GraphTraversal algoBCT = new BCTraversal();

  public static boolean hasArticulationPair(final Nodes nodes, final Edges edges) {

    return algoArticulationPair.check(nodes, edges);
  }

  public static boolean hasArticulationPoint(final Nodes nodes,
      final Edges edges) {

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

  public static boolean checkIsomorphic(final Nodes nodes, final Edges edges1,
      final Edges edges2) {

    return algoIsomorphic.check(nodes, edges1, edges2);
  }

  public static boolean traverseBF(int nodeID, final Nodes nodes,
      final Edges edges, final NodesVisitor visitor) {

    return algoBFT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseBF(final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    algoBFT.traverseAll(nodes, edges, visitor);
  }

  public static boolean traverseDF(int nodeID, final Nodes nodes,
      final Edges edges, final NodesVisitor visitor) {

    return algoDFT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseDF(final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

    algoDFT.traverseAll(nodes, edges, visitor);
  }

  public static boolean traverseBC(int nodeID, final Nodes nodes,
      final Edges edges, final NodesVisitor visitor) {

    return algoBCT.traverseComponent(nodeID, nodes, edges, visitor);
  }

  public static void traverseBC(final Nodes nodes, final Edges edges,
      final NodesVisitor visitor) {

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
 * Biconnectivity algorithm.
 * 
 * @author Demian Lessa
 */
class CheckBiconnected implements GraphProperty {

// static boolean
// isbiconnected(graph *g, int n)
// /* test if g is biconnected */
// {
// int sp,v,w;
// setword sw;
// setword visited;
// int numvis,num[MAXN],lp[MAXN],stack[MAXN];
// 
// if (n <= 2) return FALSE;
// 
// visited = bit[0];
// stack[0] = 0;
// num[0] = 0;
// lp[0] = 0;
// numvis = 1;
// sp = 0;
// v = 0;
// 
// for (;;)
// {
// if ((sw = g[v] & ~visited)) /* not "==" */
// {
// w = v;
// v = FIRSTBIT(sw); /* visit next child */
// stack[++sp] = v;
// visited |= bit[v];
// lp[v] = num[v] = numvis++;
// sw = g[v] & visited & ~bit[w];
// while (sw)
// {
// w = FIRSTBIT(sw);
// sw &= ~bit[w];
// if (num[w] < lp[v]) lp[v] = num[w];
// }
// }
// else
// {
// w = v; /* back up to parent */
// if (sp <= 1) return numvis == n;
// v = stack[--sp];
// if (lp[w] >= num[v]) return FALSE;
// if (lp[w] < lp[v]) lp[v] = lp[w];
// }
// }
// }
//
  public boolean check(Nodes nodes, Edges edges) {

    if ((nodes == null) || (edges == null)) {
      return false;
    }
    // by definition, a null graph, a singleton graph, and a two node 
    // graph are not biconnected
    if (nodes.count() < 2) {
      return false;
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
 * Connectivity algorithm. A graph is connected if there exists a path from any
 * point to any other point in the graph. A graph that is not connected is said
 * to be disconnected. This definition means that the null graph and singleton
 * graph are both connected.
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
 * Isomorphism delegation algorithm.
 * 
 * @author Demian Lessa
 */
class CheckIsomorphic implements GraphPairProperty {

  public boolean check(Nodes nodes, Edges edges1, Edges edges2) {

    // TODO Auto-generated method stub
    return false;
  }
}

/**
 * Nodal Point algorithm.
 * 
 * @author Demian Lessa
 */
class CheckNodalPoint implements GraphProperty {

  public boolean check(Nodes nodes, Edges edges) {

    // TODO Auto-generated method stub
    return false;
  }
}