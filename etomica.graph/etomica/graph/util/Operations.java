package etomica.graph.util;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.isomorphism.Match;

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

  private static String getZeroVector(int size) {

    String result = "";
    for (int i = 0; i < size; i++) {
      result += "0";
    }
    return result;
  }

  /*
   * PCopy (Positive Copy): GraphSet --> GraphSet
   *
   * Given a graph set S, returns a graph set S' such that every graph g in S has an exact
   * copy g' in S' and coefficient(g') = coefficient(g).
   */
  public static GraphSet PCopy(final GraphSet gset) {

    assert (gset != null);
    List<Edges> list = GraphFactory.createEdgesList(gset.getNodes().count());
    for (Edges e : gset.getEdgesSet()) {
      // e.copy() keeps the same coefficient sign
      list.add(e.copy());
    }
    return GraphFactory.simpleGraphSet(gset.getNodes(), list);
  }

  /*
   * NCopy (Negative Copy): GraphSet --> GraphSet
   *
   * Given a graph set S, returns a graph set S' such that every graph g in S has an exact
   * copy g' in S' and coefficient(g') = -coefficient(g).
   */
  public static GraphSet NCopy(final GraphSet gset) {

    assert (gset != null);
    List<Edges> list = GraphFactory.createEdgesList(gset.getNodes().count());
    for (Edges e : gset.getEdgesSet()) {
      // e.ncopy() reverses the coefficient sign
      list.add(e.ncopy());
    }
    return GraphFactory.simpleGraphSet(gset.getNodes(), list);
  }

  /*
   * Union (GraphSet Union): GraphSet x GraphSet --> GraphSet
   *
   * Given two graph sets S1 and S2, Union(S1, S2) is the set union of S1 and S2.
   */
  public static GraphSet Union(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    // TODO: *too weak precondition*
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    GraphSet c1 = PCopy(gset1);
    GraphSet c2 = PCopy(gset2);
    List<Edges> list = GraphFactory.createEdgesList(c1.getNodes().count());
    for (Edges e : c1.getEdgesSet()) {
      list.add(e);
    }
    for (Edges e2 : c2.getEdgesSet()) {
      list.add(e2);
    }
    return GraphFactory.simpleGraphSet(c1.getNodes(), list);
  }

  /*
   * IsoFree (Isomorph-Free GraphSet): GraphSet --> GraphSet
   *
   * Given a graph set S, IsoFree(S) returns a graph set S' such that no two graphs g1, g2
   * in S' are isomorphic to each other. For every graph g' in S', the coefficient of g'
   * is the sum of the coefficients of every graph g in S such that g is isomorphic to g'.
   */
  public static GraphSet IsoFree(final GraphSet gset) {

    assert (gset != null);
    GraphSet c1 = PCopy(gset);
    List<Edges> list = GraphFactory.createEdgesList(c1.getNodes().count());
    // assumes at first that each graph is unique, then looks for an isomorphic graph in
    // the current result and if one is found, updates the coefficient of the graph in the
    // result otherwise adds the graph to the result
    for (Edges e1 : c1.getEdgesSet()) {
      boolean unique = true;
      for (Edges e2 : list) {
        Graph g1 = GraphFactory.simpleGraph(c1.getNodes(), e1);
        Graph g2 = GraphFactory.simpleGraph(c1.getNodes(), e2);
        if (Match.match(g1, g2)) {
          Coefficient coef1 = e1.getMetadata().getCoefficient();
          Coefficient coef2 = e2.getMetadata().getCoefficient();
          e1.getMetadata().setCoefficient(coef1.add(coef2));
          unique = false;
          break;
        }
      }
      if (unique) {
        list.add(e1);
      }
    }
    return GraphFactory.simpleGraphSet(c1.getNodes(), list);
  }

  /*
   * Delete (GraphSet Deletion): GraphSet x GraphSet --> GraphSet
   *
   * Given two graph sets S1 and S2, Delete(S1, S2) is the set resulting from removing
   * every graph in S2 from the set S1.
   */
  public static GraphSet Delete(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    // TODO: *too weak precondition*
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
    for (Edges e1 : gset1.getEdgesSet()) {
      boolean unique = true;
      for (Edges e2 : gset2.getEdgesSet()) {
        if (e1.equals(e2)) {
          unique = false;
          break;
        }
      }
      if (unique) {
        list.add(e1);
      }
    }
    return GraphFactory.simpleGraphSet(gset1.getNodes(), list);
  }

  /*
   * Sum (GraphSet Sum): GraphSet x GraphSet --> GraphSet
   *
   * Given two graph sets S1 and S2, Sum(S1,S2) = IsoFree(Union(S1,S2)).
   */
  public static GraphSet Sum(final GraphSet gset1, final GraphSet gset2) {

    // TODO: this operation inherits the *too weak precondition* of the union
    return IsoFree(Union(gset1, gset2));
  }

  /*
   * Sub (GraphSet Sum): GraphSet x GraphSet --> GraphSet
   *
   * Given two graph sets S1 and S2, Sub(S1,S2) = IsoFree(Union(S1,NCopy(S2))).
   */
  public static GraphSet Sub(final GraphSet gset1, final GraphSet gset2) {

    // TODO: this operation inherits the *too weak precondition* of the union
    return IsoFree(Union(gset1, NCopy(gset2)));
  }

  /*
   * Mul (Graph Multiplication): Graph x Graph --> Graph
   */
  public static Graph Mul(final Graph g1, final Graph g2) {

    assert (g1 != null && g2 != null);
    Nodes n1 = g1.getNodes();
    Nodes n2 = g2.getNodes();
    // two root nodes with the same label in n1 and n2 must have the same color
    for (byte i = 0; i < n1.count(); i++) {
      if (i < n2.count() && n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)
          && n2.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        assert (n1.getAttributes(i).isSameColor(n2.getAttributes(i)));
      }
    }
    // two edges in g1 and g2 must not connect the same labeled root nodes
    // NOTE: expensive check!
    for (byte i = 0; i < n1.count(); i++) {
      // root nodes of g1
      if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        for (byte j = 0; j < n1.count(); j++) {
          // edges (i,j) between root nodes in g1
          if (i != j && n1.getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)
              && g1.getEdges().hasEdge(i, j)) {
            // corresponding edges (i,j) between nodes in g2
            if (i < n2.count() && j < n2.count() && g2.getEdges().hasEdge(i, j)) {
              // edge (i,j) in g2 must not be between root nodes
              assert (!n2.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES) || !n2
                  .getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES));
            }
          }
        }
      }
    }
    // multiplication: take the union of the nodes
    byte totalNodes = 0;
    byte n1RootCount = 0;
    byte n2RootCount = 0;
    byte totalRootNodes = 0;
    for (byte i = 0; i < n1.count(); i++) {
      if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        n1RootCount++;
        totalRootNodes++;
      }
      totalNodes++;
    }
    for (byte i = 0; i < n2.count(); i++) {
      if (n2.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
        if (i == totalRootNodes) {
          totalNodes++;
          totalRootNodes++;
        }
        n2RootCount++;
      }
      else {
        totalNodes++;
      }
    }
    // VERY IMPORTANT: this implementation is still representation dependent, as it
    // explores the fact that root nodes are the first nodes in the representation of a
    // cluster; in order to alleviate this issue, it is necessary to support labels on
    // top of the node representations.
    //
    // compute the mapping of field nodes in the new graph to nodes in n1 and n2
    byte fieldIndex = totalRootNodes;
    Map<Byte, Byte> fieldMap = new HashMap<Byte, Byte>();
    for (byte i = 0; i < n1.count(); i++) {
      if (n1.getAttributes(i).isSameClass(GraphFactory.FIELD_NODE_ATTRIBUTES)) {
        fieldMap.put(fieldIndex, i);
        fieldIndex++;
      }
    }
    for (byte i = 0; i < n2.count(); i++) {
      if (n2.getAttributes(i).isSameClass(GraphFactory.FIELD_NODE_ATTRIBUTES)) {
        fieldMap.put(fieldIndex, i);
        fieldIndex++;
      }
    }
    char[] rootColors = new char[totalRootNodes];
    char[] fieldColors = new char[totalNodes - totalRootNodes];
    // determine the colors of all root nodes
    for (byte i = 0; i < totalRootNodes; i++) {
      if (i < n1RootCount) {
        rootColors[i] = n1.getAttributes(i).getColor();
      }
      else {
        rootColors[i] = n2.getAttributes(i).getColor();
      }
    }
    // determine the colors of all field nodes
    for (int i = 0; i < totalNodes - totalRootNodes; i++) {
      if (i < n1.count() - n1RootCount) {
        fieldColors[i] = n1.getAttributes(fieldMap.get(i)).getColor();
      }
      else {
        fieldColors[i] = n2.getAttributes(fieldMap.get(i)).getColor();
      }
    }
    // now encode all edges in g1 and g2 to the edges in the new graph
    StringBuffer[] neighbors = new StringBuffer[totalNodes - 1];
    for (int i = 0; i < totalNodes; i++) {
      neighbors[i] = new StringBuffer(getZeroVector(totalNodes - i - 1));
    }
    for (byte i = 0; i < totalNodes; i++) {
      for (byte j = (byte) (i + 1); j < totalNodes; j++) {
        // edges from the first graph
        if (i < n1.count() && j < n1.count() && g1.getEdges().hasEdge(i, j)) {
          int fromNode = -1;
          int toNode = -1;
          if (n1.getAttributes(i).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
            fromNode = i;
          }
          else {
            fromNode = fieldMap.get(i - n1RootCount);
          }
          if (n1.getAttributes(j).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
            toNode = j;
          }
          else {
            toNode = fieldMap.get(j - n1RootCount);
          }
          neighbors[fromNode].setCharAt(toNode, '1');
        }
        // edges from the second graph
        else if (i >= n1.count() && j >= n1.count()
            && g2.getEdges().hasEdge((byte) (i - n1.count()), (byte) (j - n1.count()))) {
          int fromNode = -1;
          int toNode = -1;
          if (n2.getAttributes((byte) (i - n1.count())).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
            fromNode = i - n1.count();
          }
          else {
            fromNode = fieldMap.get(i - n1.count() - n2RootCount);
          }
          if (n2.getAttributes((byte) (j - n1.count())).isSameClass(GraphFactory.ROOT_NODE_ATTRIBUTES)) {
            toNode = j - n1.count();
          }
          else {
            toNode = fieldMap.get(j - n1.count() - n2RootCount);
          }
          neighbors[fromNode].setCharAt(toNode, '1');
        }
      }
    }
    // create the the set of nodes
    Nodes nodes = GraphFactory.defaultNodes(fieldColors, rootColors);
    // create the the set of edges
    StringBuffer strRep = new StringBuffer();
    for (int i = 0; i < totalNodes; i++) {
      strRep.append(neighbors[i]);
    }
    EdgesRepresentationFactory factory = EdgesRepresentationFactory.getFactory(nodes.count());
    Edges edges = GraphFactory.simpleEdges(factory.getRepresentation(BitmapFactory.getBitmap(strRep
        .toString())));
    Coefficient coef1 = g1.getEdges().getMetadata().getCoefficient();
    Coefficient coef2 = g2.getEdges().getMetadata().getCoefficient();
    Coefficient coef = coef1.multiply(coef2);
    edges.getMetadata().setCoefficient(coef);
    return GraphFactory.simpleGraph(nodes, edges);
  }

  public static Graph Con(final Graph g1, final Graph g2) {

    // TODO: implement
    return null;
  }

  public static GraphSet DifByNode(Graph g, char nodeColor) {

    assert (g != null);
    List<Edges> list = GraphFactory.createEdgesList(g.getNodes().count());
    for (int i = 0; i < g.getNodes().count(); i++) {
      // the nodes are different in the graph set- so cannot use the graph set 'template'
      // strategy
    }
    return null;
  }

  public static GraphSet DifByEdge(Graph g, char edgeColor) {

    assert (g != null);
    List<Edges> list = GraphFactory.createEdgesList(g.getNodes().count());
    for (int i = 0; i < g.getEdges().count(); i++) {
      // the nodes are different in the graph set- so cannot use the graph set 'template'
      // strategy
    }
    return null;
  }

  public static GraphSet Int(Graph g, int rootNodeID) {

    // TODO: this interface is incorrect; should return a single graph; once again
    // we stumble into the problem of having a graph set interface as a template
    return null;
  }

  public static GraphSet Spl(Graph g, char edgeColor, char newColor1, char newColor2) {

    // TODO: implement
    return null;
  }

  public static GraphSet Mul(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    // TODO: *too weak precondition*
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
    for (Edges e1 : gset1.getEdgesSet()) {
      for (Edges e2 : gset2.getEdgesSet()) {
        Graph g = Mul(GraphFactory.simpleGraph(gset1.getNodes(), e1), GraphFactory.simpleGraph(gset2
            .getNodes(), e2));
        list.add(g.getEdges());
      }
    }
    return IsoFree(GraphFactory.simpleGraphSet(gset1.getNodes(), list));
  }

  public static GraphSet Con(final GraphSet gset1, final GraphSet gset2) {

    assert (gset1 != null && gset2 != null);
    // TODO: *too weak precondition*
    assert (gset1.getNodes().count() == gset2.getNodes().count());
    List<Edges> list = GraphFactory.createEdgesList(gset1.getNodes().count());
    for (Edges e1 : gset1.getEdgesSet()) {
      for (Edges e2 : gset2.getEdgesSet()) {
        Graph g = Con(GraphFactory.simpleGraph(gset1.getNodes(), e1), GraphFactory.simpleGraph(gset2
            .getNodes(), e2));
        list.add(g.getEdges());
      }
    }
    return IsoFree(GraphFactory.simpleGraphSet(gset1.getNodes(), list));
  }

  public static GraphSet DifByNode(final GraphSet gset, final char nodeColor) {

    assert (gset != null);
    GraphSet result = null;
    for (Edges e : gset.getEdgesSet()) {
      GraphSet gs = DifByNode(GraphFactory.simpleGraph(gset.getNodes(), e), nodeColor);
      if (result == null) {
        result = gs;
      }
      else {
        result = IsoFree(Union(result, gs));
      }
    }
    return result;
  }

  public static GraphSet DifByEdge(final GraphSet gset, final char edgeColor) {

    assert (gset != null);
    GraphSet result = null;
    for (Edges e : gset.getEdgesSet()) {
      GraphSet gs = DifByEdge(GraphFactory.simpleGraph(gset.getNodes(), e), edgeColor);
      if (result == null) {
        result = gs;
      }
      else {
        result = IsoFree(Union(result, gs));
      }
    }
    return result;
  }

  public static GraphSet Int(final GraphSet gset, final int rootNodeID) {

    assert (gset != null);
    GraphSet result = null;
    for (Edges e : gset.getEdgesSet()) {
      GraphSet gs = Int(GraphFactory.simpleGraph(gset.getNodes(), e), rootNodeID);
      if (result == null) {
        result = gs;
      }
      else {
        result = IsoFree(Union(result, gs));
      }
    }
    return result;
  }

  public static GraphSet Pow(final GraphSet gset, final int exponent) {

    assert (gset != null);
    GraphSet result = null;
    for (int i = 0; i < exponent; i++) {
      if (result == null) {
        result = Mul(gset, gset);
      }
      else {
        result = Mul(gset, result);
      }
    }
    return result;
  }

  public static GraphSet Exp(final GraphSet gset, final int lowerExponent, final int upperExponent) {

    assert (gset != null);
    GraphSet result = null;
    // TODO: fix the coefficient by dividing by k! at each Pow operation
    for (int i = lowerExponent; i <= upperExponent; i++) {
      if (result == null) {
        result = Pow(gset, i);
      }
      else {
        result = Sum(Pow(gset, i), result);
      }
    }
    return result;
  }

  public static GraphSet Spl(final GraphSet gset, final char edgeColor, final char newColor1,
      final char newColor2) {

    assert (gset != null);
    GraphSet result = null;
    for (Edges e : gset.getEdgesSet()) {
      GraphSet gs = Spl(GraphFactory.simpleGraph(gset.getNodes(), e), edgeColor, newColor1, newColor2);
      if (result == null) {
        result = gs;
      }
      else {
        result = IsoFree(Union(result, gs));
      }
    }
    return result;
  }
}