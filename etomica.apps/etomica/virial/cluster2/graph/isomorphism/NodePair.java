/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.isomorphism;

public final class NodePair {

  private int firstNode;
  private int secondNode;
  private static NodePair NODE_PAIR_NULL = new NodePair(SearchState.NULL_NODE,
      SearchState.NULL_NODE);

  public NodePair(Integer n1, Integer n2) {

    firstNode = n1;
    secondNode = n2;
  }

  public int getN1() {

    return firstNode;
  }

  public int getN2() {

    return secondNode;
  }

  public static NodePair nullPair() {

    return NODE_PAIR_NULL;
  }

  @Override
  public String toString() {

    return "(" + firstNode + "," + secondNode + ")";
  }
}