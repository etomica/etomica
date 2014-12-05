/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.isomorphism;

import etomica.graph.model.Graph;

public class Match {

  public static String DEF_ISOMORPHISM_ALGO = SearchState.ULLMAN_ALGORITHM; // SearchState.VF2_ALGORITHM;
  // known isomorphism counts for N in {1,...,11}
  // for N > 11, the count falls beyond the range of 32 bit integers
  public static final int[] ISMORPHS_COUNT = { 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168 /**
   *
   * , 1018997864
   */
  };

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
  public static boolean match(Graph g1, Graph g2, boolean preMatch) {

    // if (calls % 2000 == 0) {
    // System.out.println(">C " + calls);
    // }
    assert (g1 != null && g2 != null);
    total++;
    if (preMatch && !g1.getSignature().equals(g2.getSignature())) {
      return false;
    }
    called++;
    SearchState state = initialState(g1, g2);
    return match(state);
  }

  /**
   * Execution of an isomorphism test with a provided initial state.
   */
  protected static boolean match(SearchState currentState) {

    assert (currentState != null);
    SearchState state = currentState;
    if (state.isGoal()) {
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