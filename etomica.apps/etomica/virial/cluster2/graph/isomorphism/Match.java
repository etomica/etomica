/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.isomorphism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import etomica.virial.cluster2.graph.Graph;
import etomica.virial.cluster2.graph.Nodes;

public class Match {

  public static String DEF_ISOMORPHISM_ALGO = SearchState.ULLMAN_ALGORITHM;
  // known isomorphism counts for N in {1,...,11}
  // for N > 11, the count falls beyond the range of 32 bit integers
  public static final int[] ISMORPHS_COUNT = { 1, 2, 4, 11, 34, 156, 1044,
      12346, 274668, 12005168 /** , 1018997864 */
  };
  /**
   * Turns on a pre-matching filter that fast prunes impossible isomorphisms.
   * NOTE: turning this off causes a significant runtime penalty in the client
   * isomorphism testing algorithms.
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
    if (g1.getNodes().count() != g2.getNodes().count()) {
      return false;
    }
    // the graphs must have the same number of edges
    if (g1.getEdges().count() != g2.getEdges().count()) {
      return false;
    }
    // the graphs must have the same sets of field node colors
    Set<Character> lf1 = g1.getNodes().getColors(Nodes.NODE_CLASS_FIELD);
    Set<Character> lf2 = g2.getNodes().getColors(Nodes.NODE_CLASS_FIELD);
    if (!lf1.equals(lf2)) {
      return false;
    }
    // the graphs must have the same sets of root node colors
    Set<Character> lr1 = g1.getNodes().getColors(Nodes.NODE_CLASS_ROOT);
    Set<Character> lr2 = g2.getNodes().getColors(Nodes.NODE_CLASS_ROOT);
    if (!lr1.equals(lr2)) {
      return false;
    }
    // the graphs must have the same sets of edge colors
    Set<Character> le1 = g1.getEdges().getColors();
    Set<Character> le2 = g2.getEdges().getColors();
    if (!le1.equals(le2)) {
      return false;
    }
    // the field nodes of both graphs must be color-partitioned identically
    for (Character color : lf1) {
      List<Integer> p1 = g1.getNodes().getPartition(Nodes.NODE_CLASS_FIELD, color);
      List<Integer> p2 = g2.getNodes().getPartition(Nodes.NODE_CLASS_FIELD, color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // the root nodes of both graphs must be color-partitioned identically
    for (Character color : lr1) {
      List<Integer> p1 = g1.getNodes().getPartition(Nodes.NODE_CLASS_ROOT, color);
      List<Integer> p2 = g2.getNodes().getPartition(Nodes.NODE_CLASS_ROOT, color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // the edges of both graphs must be color-partitioned identically
    for (Character color : le1) {
      List<Integer> p1 = g1.getEdges().getPartition(color);
      List<Integer> p2 = g2.getEdges().getPartition(color);
      if (p1.size() != p2.size()) {
        return false;
      }
    }
    // the field nodes of both graphs must have the same degrees in their partitions
    for (Character color : lf1) {
      List<Integer> p1 = g1.getNodes().getPartition(Nodes.NODE_CLASS_FIELD, color);
      List<Integer> p2 = g2.getNodes().getPartition(Nodes.NODE_CLASS_FIELD, color);
      List<Integer> d1 = new ArrayList<Integer>();
      List<Integer> d2 = new ArrayList<Integer>();
      for (Integer n1 : p1) {
        d1.add(g1.getEdges().getInDegree(n1));
      }
      for (Integer n2 : p2) {
        d2.add(g2.getEdges().getInDegree(n2));
      }
      Collections.sort(d1);
      Collections.sort(d2);
      if (!d1.equals(d2)) {
        return false;
      }
    }
    // the root nodes of both graphs must have the same degrees in their partitions
    for (Character color : lr1) {
      List<Integer> p1 = g1.getNodes().getPartition(Nodes.NODE_CLASS_ROOT, color);
      List<Integer> p2 = g2.getNodes().getPartition(Nodes.NODE_CLASS_ROOT, color);
      List<Integer> d1 = new ArrayList<Integer>();
      List<Integer> d2 = new ArrayList<Integer>();
      for (Integer n1 : p1) {
        d1.add(g1.getEdges().getInDegree(n1));
      }
      for (Integer n2 : p2) {
        d2.add(g2.getEdges().getInDegree(n2));
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