package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;


public class Match {
  

  public static boolean match(Graph g1, Graph g2) {

    VFSearchState state = new VFSearchState(g1, g2);
    if (state.isDead()) {
      return false;
    }
    if (state.isGoal()) {
      maps.add(state.getMap());
      return true;
    }
    boolean found = false;
    while (!found && state.hasNextCandidate()) {
      NodePair candidate = state.nextCandidate();
      if (state.isMatchFeasible(candidate)) {
        VFSearchState nextState = state.nextState(candidate);
        found = mapFirst(nextState);
        nextState.backTrack();
      }
    }
    return found;
  }
}
