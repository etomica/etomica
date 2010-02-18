package etomica.graph.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * TODO: some operations have a *too weak precondition* to guarantee proper operation
 * between graph sets; in particular, the labeling of the nodes to their corresponding
 * colors may not match, so when the final set FIXES the set of nodes, it explicitly
 * ignores that the 'other' graph set may have an entirely different labeling of nodes
 * to colors.
 *
 * In order to fix this problem, there has to be a validation of the graph sets before
 * the operation is computed. Either one of the graph sets must be relabeled so that
 * they are compatible, or the operation cannot be carried. This may be too expensive.
 *
 * The ideal solution is to have a canonical labeling for the graphs and graph sets.
 * This would give us two benefits: easy check of compatibility of the operation and
 * improved performance on most operations; an implicit ordering of the graphs in the
 * graph sets according to their labeling would be useful in may of the algorithms in
 * this class- it is the same principle used for the intersection of two sets: if the
 * sets are given as is, the operation is O(N*M); however, if the sets are ordered, it
 * becomes O(N+M).
 *
 * An alternative solution is to produce a map from nodes to nodes in one of the graphs
 * (basically a relabeling) and use this map before all computations involving nodes
 * and edges. This would only affect performance marginally because it would not require
 * the computation of another representation for the graph. In any case, computing this
 * relabeling is probably graph-isomorphism-complete.
 *
 * The necessity of manipulating edges in here suggests that a change in the graph API
 * is necessary. There is currently no simple way of defining an edges or nodes 'space'
 * and modifying edges or nodes at will. No nodes/edges visitation... The reason for
 * this was essentially the performance gain by a simple and small API backed by an
 * efficient representation (in terms of space). As soon as graph collections support
 * graphs with different node characteristics (colors, labels, etc), it is virtually
 * impossible to use the current approach of a template graph set strategy. The graph
 * sets share the nodes in order to save on a huge family of graphs that share the same
 * node specifications, differing only by their edges.
 *
 * For a simple graph collection (not a set template) having 2 million graphs and 7
 * nodes, the space saving is at least 14 million bytes (using the set template, one
 * encoding for the nodes is used by all graphs, whereas in a collection setting,
 * each graph has its own encoding). If we include the labeling function (a simple
 * permutation taking 7 bytes) the space would double: 28 million bytes. There is
 * also the nodes metadata, which could add up quite a bit to the equation.
 *
 * Now, comparing with the total space utilization, this might not be too bad. For
 * instance, for 7 nodes and using an upper triangle representation, we need about
 * 21 bytes per graph to encode the edges. For 2 million, that's about 43 million
 * bytes. So as a matter of fact, including nodes in the representation would have
 * a significant impact on space utilization: about 70% more space. And then some
 * more overhead for handling one graph object per graph, in addition to the edges
 * and nodes.
 *
 * The benefits of using independent nodes, edges, and graphs are significant. They
 * would allow us to code with more simplicity. Maintenance time would improve. We
 * would also be able to be more scalable. The problem with 'squeezing' performance
 * from an implementation is that several hard choices have to be made and levels
 * of abstraction have to be removed.
 *
 */

public class Operations {

//
//  public static Graph Con(final Graph g1, final Graph g2) {
//
//    // TODO: implement
//    return null;
//  }
//
//  public static GraphSet DifByNode(Graph g, char nodeColor) {
//
//    assert (g != null);
//    List<Edges> list = GraphFactory.createEdgesList(g.getNodes().count());
//    for (int i = 0; i < g.getNodes().count(); i++) {
//      // the nodes are different in the graph set- so cannot use the graph set 'template'
//      // strategy
//    }
//    return null;
//  }
//
//  public static GraphSet DifByEdge(Graph g, char edgeColor) {
//
//    assert (g != null);
//    List<Edges> list = GraphFactory.createEdgesList(g.getNodes().count());
//    for (int i = 0; i < g.getEdges().count(); i++) {
//      // the nodes are different in the graph set- so cannot use the graph set 'template'
//      // strategy
//    }
//    return null;
//  }
//
//  public static GraphSet Int(Graph g, int rootNodeID) {
//
//    // TODO: this interface is incorrect; should return a single graph; once again
//    // we stumble into the problem of having a graph set interface as a template
//    return null;
//  }
//
//  public static GraphSet Spl(Graph g, char edgeColor, char newColor1, char newColor2) {
//
//    // TODO: implement
//    return null;
//  }
//
//  public static GraphSet Mul(final GraphSet gset1, final GraphSet gset2) {
//
//    assert (gset1 != null && gset2 != null);
//    // TODO: *too weak precondition*
//    assert (gset1.getNodes().count() == gset2.getNodes().count());
//    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
//    for (Edges e1 : gset1.getEdgesSet()) {
//      for (Edges e2 : gset2.getEdgesSet()) {
//        Graph g = Mul(GraphFactory.simpleGraph(gset1.getNodes(), e1), GraphFactory.simpleGraph(gset2
//            .getNodes(), e2));
//        list.add(g.getEdges());
//      }
//    }
//    return IsoFree(GraphFactory.simpleGraphSet(gset1.getNodes(), list));
//  }
//
//  public static GraphSet Con(final GraphSet gset1, final GraphSet gset2) {
//
//    assert (gset1 != null && gset2 != null);
//    // TODO: *too weak precondition*
//    assert (gset1.getNodes().count() == gset2.getNodes().count());
//    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
//    for (Edges e1 : gset1.getEdgesSet()) {
//      for (Edges e2 : gset2.getEdgesSet()) {
//        Graph g = Con(GraphFactory.simpleGraph(gset1.getNodes(), e1), GraphFactory.simpleGraph(gset2
//            .getNodes(), e2));
//        list.add(g.getEdges());
//      }
//    }
//    return IsoFree(GraphFactory.simpleGraphSet(gset1.getNodes(), list));
//  }
//
//  public static GraphSet DifByNode(final GraphSet gset, final char nodeColor) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    for (Edges e : gset.getEdgesSet()) {
//      GraphSet gs = DifByNode(GraphFactory.simpleGraph(gset.getNodes(), e), nodeColor);
//      if (result == null) {
//        result = gs;
//      }
//      else {
//        result = IsoFree(Union(result, gs));
//      }
//    }
//    return result;
//  }
//
//  public static GraphSet DifByEdge(final GraphSet gset, final char edgeColor) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    for (Edges e : gset.getEdgesSet()) {
//      GraphSet gs = DifByEdge(GraphFactory.simpleGraph(gset.getNodes(), e), edgeColor);
//      if (result == null) {
//        result = gs;
//      }
//      else {
//        result = IsoFree(Union(result, gs));
//      }
//    }
//    return result;
//  }
//
//  public static GraphSet Int(final GraphSet gset, final int rootNodeID) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    for (Edges e : gset.getEdgesSet()) {
//      GraphSet gs = Int(GraphFactory.simpleGraph(gset.getNodes(), e), rootNodeID);
//      if (result == null) {
//        result = gs;
//      }
//      else {
//        result = IsoFree(Union(result, gs));
//      }
//    }
//    return result;
//  }
//
//  public static GraphSet Pow(final GraphSet gset, final int exponent) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    for (int i = 0; i < exponent; i++) {
//      if (result == null) {
//        result = Mul(gset, gset);
//      }
//      else {
//        result = Mul(gset, result);
//      }
//    }
//    return result;
//  }
//
//  public static GraphSet Exp(final GraphSet gset, final int lowerExponent, final int upperExponent) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    // TODO: fix the coefficient by dividing by k! at each Pow operation
//    for (int i = lowerExponent; i <= upperExponent; i++) {
//      if (result == null) {
//        result = Pow(gset, i);
//      }
//      else {
//        result = Sum(Pow(gset, i), result);
//      }
//    }
//    return result;
//  }
//
//  public static GraphSet Spl(final GraphSet gset, final char edgeColor, final char newColor1,
//      final char newColor2) {
//
//    assert (gset != null);
//    GraphSet result = null;
//    for (Edges e : gset.getEdgesSet()) {
//      GraphSet gs = Spl(GraphFactory.simpleGraph(gset.getNodes(), e), edgeColor, newColor1, newColor2);
//      if (result == null) {
//        result = gs;
//      }
//      else {
//        result = IsoFree(Union(result, gs));
//      }
//    }
//    return result;
//  }
}