package etomica.graph.isomorphism;

import etomica.graph.model.Graph;

public interface SearchState {

  public final static byte NULL_NODE = (byte) 0xFF;
  public static final String ULLMAN_ALGORITHM = "Ullman";
  public static final String VF_ALGORITHM = "VF";
  public static final String VF2_ALGORITHM = "VF2";

  public void addPair(NodePair pair);

  public void backTrack();

  public SearchState copy();

  public void copy(SearchState fromState);

  public NodePair[] getCoreSet();

  public byte getCoreLen();

  public Graph getG1();

  public Graph getG2();

  public boolean isDead();

  public boolean isFeasiblePair(NodePair pair);

  public boolean isGoal();

  public NodePair nextPair(NodePair prev);
}