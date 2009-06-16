package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;
import etomica.virial.cluster2.graph.Tagged;

public interface SearchState extends Tagged {

  public final static int NULL_NODE = 0xFFFF;
  public static final String ULLMAN_ALGORITHM = "Ullman";
  public static final String VF_ALGORITHM = "VF";
  public static final String VF2_ALGORITHM = "VF2";

  public void addPair(NodePair pair);

  public void backTrack();

  public SearchState copy();

  public void copy(SearchState fromState);

  public NodePair[] getCoreSet();

  public int getCoreLen();

  public Graph getG1();

  public Graph getG2();

  public boolean isDead();

  public boolean isFeasiblePair(NodePair pair);

  public boolean isGoal();

  public NodePair nextPair(NodePair prev);
}