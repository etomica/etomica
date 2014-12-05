/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.isomorphism;

import etomica.virial.cluster2.graph.Graph;

public class UllmanSearchState extends AbstractSearchState {

  private boolean[][] M; // Matrix encoding the compatibility of the nodes

  public UllmanSearchState(Graph g1, Graph g2) {

    super(g1, g2);
    n1 = getN1().count();
    n2 = getN2().count();
    core_len = 0;
    core_1 = new int[n1];
    core_2 = new int[n2];
    M = new boolean[n1][];
    for (int i = 0; i < n1; i++) {
      M[i] = new boolean[n2];
    }
    for (int i = 0; i < n1; i++) {
      core_1[i] = NULL_NODE;
    }
    for (int i = 0; i < n2; i++) {
      core_2[i] = NULL_NODE;
    }
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        M[i][j] = getE1().getInDegree(i) == getE2().getInDegree(j)
            && getE1().getOutDegree(i) == getE2().getOutDegree(j)
            && getN1Attrs(i).isCompatible(getN2Attrs(j));
      }
    }
  }

  public UllmanSearchState(UllmanSearchState state) {

    copy(state);
  }

  public void addPair(NodePair pair) {

    assert (pair != null);
    int node1 = pair.getN1();
    int node2 = pair.getN2();
    assert (node1 < n1);
    assert (node2 < n2);
    assert (core_len < n1);
    assert (core_len < n2);
    core_1[node1] = node2;
    core_2[node2] = node1;
    core_len++;
    for (int k = core_len; k < n1; k++) {
      M[k][node2] = false;
    }
    refine();
  }

  public SearchState copy() {

    return new UllmanSearchState(this);
  }

  public void copy(SearchState fromState) {

    assert (fromState != null);
    assert (fromState instanceof UllmanSearchState);
    UllmanSearchState state = (UllmanSearchState) fromState;
    setG1(fromState.getG1());
    setG2(fromState.getG2());
    n1 = state.n1;
    n2 = state.n2;
    core_len = state.core_len;
    core_1 = new int[n1];
    core_2 = new int[n2];
    M = new boolean[n1][];
    for (int i = 0; i < core_len; i++) {
      M[i] = null;
    }
    for (int i = core_len; i < n1; i++) {
      M[i] = new boolean[n2];
    }
    for (int i = 0; i < n1; i++) {
      core_1[i] = state.core_1[i];
    }
    for (int i = 0; i < n2; i++) {
      core_2[i] = state.core_2[i];
    }
    for (int i = core_len; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        M[i][j] = state.M[i][j];
      }
    }
  }

  public boolean isDead() {

    if (n1 != n2) {
      return true;
    }
    outer: for (int i = core_len; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        if (M[i][j]) {
          continue outer;
        }
      }
      return true;
    }
    return false;
  }

  public boolean isFeasiblePair(NodePair pair) {

    assert (pair != null);
    int node1 = pair.getN1();
    int node2 = pair.getN2();
    assert (node1 < n1);
    assert (node2 < n2);
    return M[node1][node2];
  }

  public boolean isGoal() {

    return (core_len == n1) && (core_len == n2);
  }

  public NodePair nextPair(NodePair prev) {

    assert (prev != null);
    int prev_n1 = prev.getN1();
    int prev_n2 = prev.getN2();
    if (prev_n1 == NULL_NODE) {
      prev_n1 = core_len;
      prev_n2 = 0;
    }
    else if (prev_n2 == NULL_NODE) {
      prev_n2 = 0;
    }
    else {
      prev_n2++;
    }
    if (prev_n2 >= n2) {
      prev_n1++;
      prev_n2 = 0;
    }
    if (prev_n1 != core_len) {
      return null;
    }
    while (prev_n2 < n2 && !M[prev_n1][prev_n2]) {
      prev_n2++;
    }
    if (prev_n2 < n2) {
      return new NodePair(prev_n1, prev_n2);
    }
    else {
      return null;
    }
  }

  private void refine() {

    int i, j, k, l;
    for (i = core_len; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (M[i][j]) {
          boolean edge_ik, edge_ki, edge_jl, edge_lj;
          // The following (commented out) for wasn't necessary...
          // for(k=0; k<core_len; k++)
          for (k = core_len - 1; k < core_len; k++) {
            l = core_1[k];
            assert (l != NULL_NODE);
            edge_ik = getE1().hasEdge(i, k);
            edge_ki = getE1().hasEdge(k, i);
            edge_jl = getE2().hasEdge(j, l);
            edge_lj = getE2().hasEdge(l, j);
            if (edge_ik != edge_jl || edge_ki != edge_lj) {
              M[i][j] = false;
              break;
            }
            else if (edge_ik
                && !getE1Attrs(i, k).isCompatible(getE2Attrs(j, l))) {
              M[i][j] = false;
              break;
            }
            else if (edge_ki
                && !getE1Attrs(k, i).isCompatible(getE2Attrs(l, j))) {
              M[i][j] = false;
              break;
            }
          }
        }
      }
    }
  }

  @Override
  protected String getTag() {

    return SearchState.ULLMAN_ALGORITHM;
  }
}