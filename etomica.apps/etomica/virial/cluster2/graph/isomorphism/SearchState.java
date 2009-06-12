package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;

public interface SearchState {

  public final static int NULL_NODE = 0xFFFF;

  public void addPair(NodePair pair);

  public void backTrack();

  public SearchState copy();

  public NodePair[] getCoreSet();

  public int getCoreLen();

  public Graph getG1();

  public Graph getG2();

  public boolean isDead();

  public boolean isFeasiblePair(NodePair pair);

  public boolean isGoal();

  public NodePair nextPair(NodePair prev);
}