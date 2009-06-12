package etomica.virial.cluster2.graph.isomorphism;

import java.util.ArrayList;

import etomica.virial.cluster2.graph.Graph;

public class VF2SearchState extends AbstractSearchState {

  private int core_len, orig_core_len;
  private int added_node1;
  private int t1both_len, t2both_len;
  private int t1in_len, t1out_len;
  private int t2in_len, t2out_len;
  private int n1, n2;
  private int[] core_1, core_2;
  private int[] in_1, in_2;
  private int[] out_1, out_2;
  private long share_count;

  public VF2SearchState(Graph g1, Graph g2) {

    super(g1, g2);
    n1 = g1.getNodes().count();
    n2 = g2.getNodes().count();
    core_len = orig_core_len = 0;
    t1both_len = t1in_len = t1out_len = 0;
    t2both_len = t2in_len = t2out_len = 0;
    added_node1 = NULL_NODE;
    core_1 = new int[n1];
    core_2 = new int[n2];
    in_1 = new int[n1];
    in_2 = new int[n2];
    out_1 = new int[n1];
    out_2 = new int[n2];
    share_count = 0;
    int i;
    for (i = 0; i < n1; i++) {
      core_1[i] = NULL_NODE;
      in_1[i] = 0;
      out_1[i] = 0;
    }
    for (i = 0; i < n2; i++) {
      core_2[i] = NULL_NODE;
      in_2[i] = 0;
      out_2[i] = 0;
    }
    share_count = 1;
  }

  public VF2SearchState(VF2SearchState state) {

    super(state.getG1(), state.getG2());
    n1 = state.n1;
    n2 = state.n2;
    core_len = orig_core_len = state.core_len;
    t1in_len = state.t1in_len;
    t1out_len = state.t1out_len;
    t1both_len = state.t1both_len;
    t2in_len = state.t2in_len;
    t2out_len = state.t2out_len;
    t2both_len = state.t2both_len;
    added_node1 = NULL_NODE;
    core_1 = state.core_1;
    core_2 = state.core_2;
    in_1 = state.in_1;
    in_2 = state.in_2;
    out_1 = state.out_1;
    out_2 = state.out_2;
    share_count = state.share_count;
    share_count++;
  }

  public void addPair(NodePair pair) {

    assert (pair != null);
    int node1 = pair.getFirstNode();
    int node2 = pair.getSecondNode();
    assert (node1 < n1);
    assert (node2 < n2);
    assert (core_len < n1);
    assert (core_len < n2);
    core_len++;
    added_node1 = node1;
    if (in_1[node1] == 0) {
      in_1[node1] = core_len;
      t1in_len++;
      if (out_1[node1] != 0) {
        t1both_len++;
      }
    }
    if (out_1[node1] == 0) {
      out_1[node1] = core_len;
      t1out_len++;
      if (in_1[node1] != 0) {
        t1both_len++;
      }
    }
    if (in_2[node2] == 0) {
      in_2[node2] = core_len;
      t2in_len++;
      if (out_2[node2] != 0) {
        t2both_len++;
      }
    }
    if (out_2[node2] == 0) {
      out_2[node2] = core_len;
      t2out_len++;
      if (in_2[node2] != 0) {
        t2both_len++;
      }
    }
    // update the matched core set
    core_1[node1] = node2;
    core_2[node2] = node1;
    int i, other;
    for (i = 0; i < getG1().getEdges().getInDegree(node1); i++) {
      other = getG1().getEdges().getInNode(node1, i);
      if (in_1[other] == 0) {
        in_1[other] = core_len;
        t1in_len++;
        if (out_1[other] != 0) {
          t1both_len++;
        }
      }
    }
    for (i = 0; i < getG1().getEdges().getOutDegree(node1); i++) {
      other = getG1().getEdges().getOutNode(node1, i);
      if (out_1[other] == 0) {
        out_1[other] = core_len;
        t1out_len++;
        if (in_1[other] != 0) {
          t1both_len++;
        }
      }
    }
    for (i = 0; i < getG2().getEdges().getInDegree(node2); i++) {
      other = getG2().getEdges().getInNode(node2, i);
      if (in_2[other] == 0) {
        in_2[other] = core_len;
        t2in_len++;
        if (out_2[other] != 0) {
          t2both_len++;
        }
      }
    }
    for (i = 0; i < getG2().getEdges().getOutDegree(node2); i++) {
      other = getG2().getEdges().getOutNode(node2, i);
      if (out_2[other] == 0) {
        out_2[other] = core_len;
        t2out_len++;
        if (in_2[other] != 0) {
          t2both_len++;
        }
      }
    }
  }

  public void backTrack() {

    assert (core_len - orig_core_len <= 1);
    assert (added_node1 != NULL_NODE);
    if (orig_core_len < core_len) {
      int i, node2;
      if (in_1[added_node1] == core_len)
        in_1[added_node1] = 0;
      for (i = 0; i < getG1().getEdges().getInDegree(added_node1); i++) {
        int other = getG1().getEdges().getInNode(added_node1, i);
        if (in_1[other] == core_len) {
          in_1[other] = 0;
        }
      }
      if (out_1[added_node1] == core_len)
        out_1[added_node1] = 0;
      for (i = 0; i < getG1().getEdges().getOutDegree(added_node1); i++) {
        int other = getG1().getEdges().getOutNode(added_node1, i);
        if (out_1[other] == core_len) {
          out_1[other] = 0;
        }
      }
      node2 = core_1[added_node1];
      if (in_2[node2] == core_len)
        in_2[node2] = 0;
      for (i = 0; i < getG2().getEdges().getInDegree(node2); i++) {
        int other = getG2().getEdges().getInNode(node2, i);
        if (in_2[other] == core_len) {
          in_2[other] = 0;
        }
      }
      if (out_2[node2] == core_len)
        out_2[node2] = 0;
      for (i = 0; i < getG2().getEdges().getOutDegree(node2); i++) {
        int other = getG2().getEdges().getOutNode(node2, i);
        if (out_2[other] == core_len) {
          out_2[other] = 0;
        }
      }
      core_1[added_node1] = NULL_NODE;
      core_2[node2] = NULL_NODE;
      core_len = orig_core_len;
      added_node1 = NULL_NODE;
    }
  }

  public SearchState copy() {

    return new VF2SearchState(this);
  }

  public int getCoreLen() {

    return core_len;
  }

  public NodePair[] getCoreSet() {

    ArrayList<NodePair> pairList = new ArrayList<NodePair>();
    for (int i = 0; i < n1; i++) {
      if (core_1[i] != NULL_NODE) {
        pairList.add(new NodePair(i, core_1[i]));
      }
    }
    return pairList.toArray(new NodePair[] {});
  }

  public boolean isDead() {

    return (n1 != n2) || (t1both_len != t2both_len) || (t1out_len != t2out_len)
        || (t1in_len != t2in_len);
  }

  public boolean isFeasiblePair(NodePair pair) {

    assert (pair != null);
    int node1 = pair.getFirstNode();
    int node2 = pair.getSecondNode();
    assert (node1 < n1);
    assert (node2 < n2);
    assert (core_1[node1] == NULL_NODE);
    assert (core_2[node2] == NULL_NODE);
    if (!getG1().getNodes().getAttributes(node1).isCompatible(
        getG2().getNodes().getAttributes(node2))) {
      return false;
    }
    int i, other1, other2;
    int termout1 = 0, termout2 = 0, termin1 = 0, termin2 = 0, new1 = 0, new2 = 0;
    // Check the 'out' edges of node1
    for (i = 0; i < getG1().getEdges().getOutDegree(node1); i++) {
      other1 = getG1().getEdges().getOutNode(node1, i);
      if (core_1[other1] != NULL_NODE) {
        other2 = core_1[other1];
        if (!getG2().getEdges().hasEdge(node2, other2)
            || !getG1().getEdges().getAttributes(node1, other1).isCompatible(
                getG2().getEdges().getAttributes(node2, other2))) {
          return false;
        }
      }
      else {
        if (in_1[other1] != 0) {
          termin1++;
        }
        if (out_1[other1] != 0) {
          termout1++;
        }
        if (in_1[other1] == 0 && out_1[other1] == 0) {
          new1++;
        }
      }
    }
    // Check the 'in' edges of node1
    for (i = 0; i < getG1().getEdges().getInDegree(node1); i++) {
      other1 = getG1().getEdges().getInNode(node1, i);
      if (core_1[other1] != NULL_NODE) {
        other2 = core_1[other1];
        if (!getG2().getEdges().hasEdge(other2, node2)
            || !getG1().getEdges().getAttributes(node1, other1).isCompatible(
                getG2().getEdges().getAttributes(other2, node2))) {
          return false;
        }
      }
      else {
        if (in_1[other1] != 0) {
          termin1++;
        }
        if (out_1[other1] != 0) {
          termout1++;
        }
        if (in_1[other1] == 0 && out_1[other1] == 0) {
          new1++;
        }
      }
    }
    // Check the 'out' edges of node2
    for (i = 0; i < getG2().getEdges().getOutDegree(node2); i++) {
      other2 = getG2().getEdges().getOutNode(node2, i);
      if (core_2[other2] != NULL_NODE) {
        other1 = core_2[other2];
        if (!getG1().getEdges().hasEdge(node1, other1)) {
          return false;
        }
      }
      else {
        if (in_2[other2] != 0) {
          termin2++;
        }
        if (out_2[other2] != 0) {
          termout2++;
        }
        if (in_2[other2] == 0 && out_2[other2] == 0) {
          new2++;
        }
      }
    }
    // Check the 'in' edges of node2
    for (i = 0; i < getG2().getEdges().getInDegree(node2); i++) {
      other2 = getG2().getEdges().getInNode(node2, i);
      if (core_2[other2] != NULL_NODE) {
        other1 = core_2[other2];
        if (!getG1().getEdges().hasEdge(other1, node1)) {
          return false;
        }
      }
      else {
        if (in_2[other2] != 0) {
          termin2++;
        }
        if (out_2[other2] != 0) {
          termout2++;
        }
        if (in_2[other2] == 0 && out_2[other2] == 0) {
          new2++;
        }
      }
    }
    return termin1 == termin2 && termout1 == termout2 && new1 == new2;
  }

  public boolean isGoal() {

    return (core_len == n1) && (core_len == n2);
  }

  public NodePair nextPair(NodePair prev) {

    assert (prev != null);
    int prev_n1 = prev.getFirstNode();
    int prev_n2 = prev.getSecondNode();
    if (prev_n1 == NULL_NODE)
      prev_n1 = 0;
    if (prev_n2 == NULL_NODE)
      prev_n2 = 0;
    else
      prev_n2++;
    if (t1both_len > core_len && t2both_len > core_len) {
      while (prev_n1 < n1
          && (core_1[prev_n1] != NULL_NODE || out_1[prev_n1] == 0 || in_1[prev_n1] == 0)) {
        prev_n1++;
        prev_n2 = 0;
      }
    }
    else if (t1out_len > core_len && t2out_len > core_len) {
      while (prev_n1 < n1
          && (core_1[prev_n1] != NULL_NODE || out_1[prev_n1] == 0)) {
        prev_n1++;
        prev_n2 = 0;
      }
    }
    else if (t1in_len > core_len && t2in_len > core_len) {
      while (prev_n1 < n1
          && (core_1[prev_n1] != NULL_NODE || in_1[prev_n1] == 0)) {
        prev_n1++;
        prev_n2 = 0;
      }
    }
// else if (prev_n1==0 && order!=NULL)
// { int i=0;
// while (i<n1 && core_1[prev_n1=order[i]]!=NULL_NODE)
// i++;
// if (i==n1)
// prev_n1=n1;
// }
    else {
      while (prev_n1 < n1 && core_1[prev_n1] != NULL_NODE) {
        prev_n1++;
        prev_n2 = 0;
      }
    }
    if (t1both_len > core_len && t2both_len > core_len) {
      while (prev_n2 < n2
          && (core_2[prev_n2] != NULL_NODE || out_2[prev_n2] == 0 || in_2[prev_n2] == 0)) {
        prev_n2++;
      }
    }
    else if (t1out_len > core_len && t2out_len > core_len) {
      while (prev_n2 < n2
          && (core_2[prev_n2] != NULL_NODE || out_2[prev_n2] == 0)) {
        prev_n2++;
      }
    }
    else if (t1in_len > core_len && t2in_len > core_len) {
      while (prev_n2 < n2
          && (core_2[prev_n2] != NULL_NODE || in_2[prev_n2] == 0)) {
        prev_n2++;
      }
    }
    else {
      while (prev_n2 < n2 && core_2[prev_n2] != NULL_NODE) {
        prev_n2++;
      }
    }
    if (prev_n1 < n1 && prev_n2 < n2) {
      return new NodePair(prev_n1, prev_n2);
    }
    return null;
  }
}