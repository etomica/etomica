/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.isomorphism;

import etomica.graph.model.Graph;
import etomica.graph.model.Node;

public class UllmanSearchState extends AbstractSearchState {

  private boolean[][] M; // Matrix encoding the compatibility of the nodes

  public UllmanSearchState(Graph g1, Graph g2) {

    super(g1, g2);
    n1 = g1.nodeCount();
    n2 = g1.nodeCount();
    core_len = 0;
    core_1 = new byte[n1];
    core_2 = new byte[n2];
    M = new boolean[n1][];
    for (byte i = 0; i < n1; i++) {
      M[i] = new boolean[n2];
    }
    for (byte i = 0; i < n1; i++) {
      core_1[i] = NULL_NODE;
    }
    for (byte i = 0; i < n2; i++) {
      core_2[i] = NULL_NODE;
    }
    for (byte i = 0; i < n1; i++) {
      byte odi = g1.getOutDegree(i);
      Node iNode = getN1(i);
      for (byte j = 0; j < n2; j++) {
        M[i][j] = odi == g2.getOutDegree(j) && iNode.isCompatible(getN2(j));
      }
    }
  }

  public UllmanSearchState(UllmanSearchState state) {

    copy(state);
  }

  public void addPair(NodePair pair) {

    assert (pair != null);
    byte node1 = pair.getN1();
    byte node2 = pair.getN2();
    assert (node1 < n1);
    assert (node2 < n2);
    assert (core_len < n1);
    assert (core_len < n2);
    core_1[node1] = node2;
    core_2[node2] = node1;
    core_len++;
    for (byte k = core_len; k < n1; k++) {
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
    core_1 = new byte[n1];
    core_2 = new byte[n2];
    M = new boolean[n1][];
    for (byte i = 0; i < core_len; i++) {
      M[i] = null;
    }
    for (byte i = core_len; i < n1; i++) {
      M[i] = new boolean[n2];
    }
    for (byte i = 0; i < n1; i++) {
      core_1[i] = state.core_1[i];
    }
    for (byte i = 0; i < n2; i++) {
      core_2[i] = state.core_2[i];
    }
    for (byte i = core_len; i < n1; i++) {
      for (byte j = 0; j < n2; j++) {
        M[i][j] = state.M[i][j];
      }
    }
  }

  public boolean isDead() {

    if (n1 != n2) {
      return true;
    }
    outer: for (byte i = core_len; i < n1; i++) {
      for (byte j = 0; j < n2; j++) {
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
    byte node1 = pair.getN1();
    byte node2 = pair.getN2();
    assert (node1 < n1);
    assert (node2 < n2);
    return M[node1][node2];
  }

  public boolean isGoal() {

    return (core_len == n1) && (core_len == n2);
  }

  public NodePair nextPair(NodePair prev) {

    assert (prev != null);
    byte prev_n1 = prev.getN1();
    byte prev_n2 = prev.getN2();
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

    byte i, j, k, l;
    for (i = core_len; i < n1; i++) {
      for (j = 0; j < n2; j++) {
        if (M[i][j]) {
          boolean edge_ik, edge_jl; // , edge_ki, edge_lj;
          // The following (commented out) for wasn't necessary...
          // for(k=0; k<core_len; k++)
          for (k = (byte) (core_len - 1); k < core_len; k++) {
            l = core_1[k];
            assert (l != NULL_NODE);
            edge_ik = getG1().hasEdge(i, k);
            // edge_ki = getG1().hasEdge(k, i);
            edge_jl = getG2().hasEdge(j, l);
            // edge_lj = getG2().hasEdge(l, j);
            if (edge_ik != edge_jl) {// || edge_ki != edge_lj) {
              M[i][j] = false;
              break;
            }
            else if (edge_ik && !getE1(i, k).isCompatible(getE2(j, l))) {
              M[i][j] = false;
              break;
            }
            // else if (edge_ki && !getE1Attrs(k, i).isCompatible(getE2Attrs(l, j))) {
            // M[i][j] = false;
            // break;
            // }
          }
        }
      }
    }
  }
}