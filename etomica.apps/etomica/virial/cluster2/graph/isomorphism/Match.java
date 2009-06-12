package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;

public class Match {

  public static boolean match(SearchState initialState) {

    assert (initialState != null);
    return match(initialState, NodePair.nullPair());
  }

  protected static boolean match(SearchState currentState, NodePair currentPair) {

    assert (currentState != null);
    SearchState state = currentState;
    assert (currentPair != null);
    assert (state.getG1() != null);
    assert (state.getG2() != null);
    Graph g1 = state.getG1();
    Graph g2 = state.getG2();
    assert (g1.getNodes().count() == g2.getNodes().count());
    if (state.isGoal()) {
      return true;
    }
    if (state.isDead()) {
      return false;
    }
    boolean found = false;
    while (!found) {
      NodePair pair = state.nextPair(currentPair);
      if (pair != null && state.isFeasiblePair(pair)) {
        SearchState nextState = state.copy();
        nextState.addPair(pair);
        found = match(nextState, pair);
        nextState.backTrack();
      }
    }
    return found;
  }
}