package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;

public interface SearchState {

  public void addPair(int node1, int node2);

  public void backTrack();

  public SearchState copy();

  public NodePair[] getCoreSet();

  public int getCoreLen();

  public Graph getG1();

  public Graph getG2();

  public boolean isDead();

  public boolean isFeasiblePair(int node1, int node2);

  public boolean isGoal();

  public NodePair nextPair(int prev_n1, int prev_n2);
}