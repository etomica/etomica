package etomica.graph.isomorphism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import etomica.graph.model.Graph;
import static etomica.graph.model.Metadata.*;

public class Match {

  public static String DEF_ISOMORPHISM_ALGO = SearchState.ULLMAN_ALGORITHM;
  // known isomorphism counts for N in {1,...,11}
  // for N > 11, the count falls beyond the range of 32 bit integers
  public static final int[] ISMORPHS_COUNT = { 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168 /**
   *
   * , 1018997864
   */
  };
  /**
   * Turns on a pre-matching filter that fast prunes impossible isomorphisms. NOTE:
   * turning this off causes a significant runtime penalty in the client isomorphism
   * testing algorithms.
   */
  public static final boolean PRE_MATCH = true;
  public static int called = 0;
  public static int total = 0;
  public static int graphs = 0;
  public static int oldGraphs = 0;

  /**
   * Initial search state with the default isomorphism test algorithm.
   */
  public static SearchState initialState(Graph g1, Graph g2) {

    if (DEF_ISOMORPHISM_ALGO.equals(SearchState.VF_ALGORITHM)) {
      return new VFSearchState(g1, g2);
    }
    else if (DEF_ISOMORPHISM_ALGO.equals(SearchState.VF2_ALGORITHM)) {
      return new VF2SearchState(g1, g2);
    }
    else if (DEF_ISOMORPHISM_ALGO.equals(SearchState.ULLMAN_ALGORITHM)) {
      return new UllmanSearchState(g1, g2);
    }
    return null;
  }

  /**
   * Execution of an isomorphism test with the default algorithm.
   */
  public static boolean match(Graph g1, Graph g2) {

    // if (calls % 2000 == 0) {
    // System.out.println(">C " + calls);
    // }
    assert (g1 != null && g2 != null);
    total++;
    if (PRE_MATCH && !preMatch(g1, g2)) {
      return false;
    }
    called++;
    SearchState state = initialState(g1, g2);
    return match(state);
  }

  protected static boolean preMatch(Graph g1, Graph g2) {

    // the graphs must have the same number of nodes
    if (g1.getNodeCount() != g2.getNodeCount()) {
      return false;
    }
    // the graphs must have the same number of edges
    if (g1.getEdgeCount() != g2.getEdgeCount()) {
      return false;
    }
    // short circuit nodes tests
    if (g1.getNodeCount() == 0) {
      return true;
    }
    // the graphs must have the same sets of field node colors
    Set<Character> lf1 = g1.getColors(TYPE_NODE_FIELD);
    Set<Character> lf2 = g2.getColors(TYPE_NODE_FIELD);
    if (!lf1.equals(lf2)) {
      return false;
    }
    // the graphs must have the same sets of root node colors
    Set<Character> lr1 = g1.getColors(TYPE_NODE_ROOT);
    Set<Character> lr2 = g2.getColors(TYPE_NODE_ROOT);
    if (!lr1.equals(lr2)) {
      return false;
    }
    // the field nodes of both graphs must be color-partitioned identically
    for (Character color : lf1) {
      List<Byte> p1 = g1.getPartition(TYPE_NODE_FIELD, color);
      List<Byte> p2 = g2.getPartition(TYPE_NODE_FIELD, color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // the root nodes of both graphs must be color-partitioned identically
    for (Character color : lr1) {
      List<Byte> p1 = g1.getPartition(TYPE_NODE_ROOT, color);
      List<Byte> p2 = g2.getPartition(TYPE_NODE_ROOT, color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // short circuit edges tests
    if (g1.getEdgeCount() == 0) {
      return true;
    }
    // the graphs must have the same sets of edge colors
    Set<Character> le1 = g1.getColors(TYPE_EDGE_ANY);
    Set<Character> le2 = g2.getColors(TYPE_EDGE_ANY);
    if (!le1.equals(le2)) {
      return false;
    }
    // the edges of both graphs must be color-partitioned identically
    for (Character color : le1) {
      List<Byte> p1 = g1.getPartition(TYPE_EDGE_ANY, color);
      List<Byte> p2 = g2.getPartition(TYPE_EDGE_ANY, color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // the field nodes of both graphs must have the same degrees in their partitions
    for (Character color : lf1) {
      List<Byte> p1 = g1.getPartition(TYPE_NODE_FIELD, color);
      List<Byte> p2 = g2.getPartition(TYPE_NODE_FIELD, color);
      List<Byte> d1 = new ArrayList<Byte>();
      List<Byte> d2 = new ArrayList<Byte>();
      for (Byte n1 : p1) {
        d1.add(g1.getOutDegree(n1));
      }
      for (Byte n2 : p2) {
        d2.add(g2.getOutDegree(n2));
      }
      Collections.sort(d1);
      Collections.sort(d2);
      if (!d1.equals(d2)) {
        return false;
      }
    }
    // the root nodes of both graphs must have the same degrees in their partitions
    for (Character color : lr1) {
      List<Byte> p1 = g1.getPartition(TYPE_NODE_ROOT, color);
      List<Byte> p2 = g2.getPartition(TYPE_NODE_ROOT, color);
      List<Byte> d1 = new ArrayList<Byte>();
      List<Byte> d2 = new ArrayList<Byte>();
      for (Byte n1 : p1) {
        d1.add(g1.getOutDegree(n1));
      }
      for (Byte n2 : p2) {
        d2.add(g2.getOutDegree(n2));
      }
      Collections.sort(d1);
      Collections.sort(d2);
      if (!d1.equals(d2)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Execution of an isomorphism test with a user provided initial state.
   */
  protected static boolean match(SearchState currentState) {

    assert (currentState != null);
    SearchState state = currentState;
    if (state.isGoal()) {
      // System.out.println("***** MATCH *****");
      // System.out.println(state.getG1().getEdges());
      // System.out.println("======");
      // System.out.println(state.getG2().getEdges());
      // System.out.println("======");
      // for (NodePair p : state.getCoreSet()) {
      // System.out.println(p);
      // }
      // System.out.println("*****************");
      return true;
    }
    if (state.isDead()) {
      return false;
    }
    boolean found = false;
    NodePair pair = state.nextPair(NodePair.nullPair());
    while (!found && pair != null) {
      if (state.isFeasiblePair(pair)) {
        SearchState nextState = state.copy();
        nextState.addPair(pair);
        found = match(nextState);
        if (!found) {
          nextState.backTrack();
        }
      }
      pair = state.nextPair(pair);
    }
    return found;
  }
}