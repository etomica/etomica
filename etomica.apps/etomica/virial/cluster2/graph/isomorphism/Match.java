package etomica.virial.cluster2.graph.isomorphism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import etomica.virial.cluster2.graph.Graph;

public class Match {

  public static String DEF_ISOMORPHISM_ALGO = SearchState.ULLMAN_ALGORITHM;
  // known isomorphism counts for N in {1,...,11}
  // for N > 11, the count falls beyond the range of 32 bit integers
  public static final int[] ISMORPHS_COUNT = { 1, 2, 4, 11, 34, 156, 1044,
      12346, 274668, 12005168 /** , 1018997864 */
  };
  // these are the counts for N in {1,...,11} considering
  // only the generation of graphs with 0..bin(n,2)/2 edges
  // so that the complement operation (much cheaper) can
  // be then applied onto the graph set:
  //
  // edges......: (1,0), (2,0), (3,1), (4,3), (5, 5), (6, 7), (7, 10),
  // (8, 14), (9, 18), (10, 22), (11, 27)
  // graphs.....: (1,1), (2,1), (3,2), (4,7), (5,20), (6,78), (7,522),
  // (8,6996), (9,154354), (10,6002584), ...
  // complement.: (1,1), (2,1), (3,2), (4,4), (5,14), (6,78), (7,522),
  // (8,5350), (9,120314), (10,6002584), ...
  public static final int[] OPTIMAL_ISMORPHS_COUNT = { 1, 1, 2, 7, 20, 78, 522,
      6996, 154354, 6002584 };
  /**
   * Turns on a pre-matching filter that fast prunes impossible isomorphisms.
   * NOTE: turning this off causes a significant runtime penalty in the client
   * isomorphism testing algorithms.
   */
  public static final boolean PRE_MATCH = true;
  private static final boolean PRE_MAP = true;
  public static int calls = 0;
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

    calls++;
// if (calls % 2000 == 0) {
// System.out.println(">C " + calls);
// }
    assert (g1 != null && g2 != null);
    if (PRE_MATCH && !preMatch(g1, g2)) {
      return false;
    }
    SearchState state = initialState(g1, g2);
    /**
     * Try to pre-map all labeled nodes (singleton partitions).
     */
    if (PRE_MAP && g1.getNodes().getPartitionCount() > 1) {
      // the last partition is has unlabeled nodes
      for (int pID = 0; pID < g1.getNodes().getPartitionCount() - 1; pID++) {
        List<Integer> p1 = g1.getNodes().getPartition(pID);
        List<Integer> p2 = g2.getNodes().getPartition(pID);
        assert (p1.size() == 1 && p2.size() == 1);
        NodePair pair = new NodePair(p1.get(0), p2.get(0));
        if (!state.isFeasiblePair(pair)) {
          return false;
        }
        state.addPair(pair);
      }
    }
    /*
     * If we reached this far, we have now to perform the actual match.
     */
    return match(state);
  }

  protected static boolean preMatch(Graph g1, Graph g2) {

    // Isomorphic graphs must have the same number of nodes.
    if (g1.getNodes().count() != g2.getNodes().count()) {
      return false;
    }
    // Isomorphic graphs must have the same number of edges.
    if (g1.getEdges().count() != g2.getEdges().count()) {
      return false;
    }
    // Isomorphic graphs must have the same number of node partitions.
    if (g1.getNodes().getPartitionCount() != g2.getNodes().getPartitionCount()) {
      return false;
    }
    /**
     * Isomorphic graphs must have the same node degree sets within each
     * partition.
     */
    for (int pID = 0; pID < g1.getNodes().getPartitionCount(); pID++) {
      List<Integer> degs1 = new ArrayList<Integer>();
      List<Integer> degs2 = new ArrayList<Integer>();
      for (Integer n1 : g1.getNodes().getPartition(pID)) {
        degs1.add(g1.getEdges().getInDegree(n1));
      }
      for (Integer n2 : g2.getNodes().getPartition(pID)) {
        degs2.add(g2.getEdges().getInDegree(n2));
      }
      Collections.sort(degs1);
      Collections.sort(degs2);
      if (!degs1.equals(degs2)) {
        return false;
      }
    }
// for (int node = 0; node < g1.getNodes().count(); node++) {
// degs1.add(g1.getEdges().getInDegree(node));
// degs2.add(g2.getEdges().getInDegree(node));
// }
// Collections.sort(degs1);
// Collections.sort(degs2);
// if (!degs1.equals(degs2)) {
// return false;
// }
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