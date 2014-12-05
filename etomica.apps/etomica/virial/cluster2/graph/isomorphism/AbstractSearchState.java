/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.isomorphism;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import etomica.virial.cluster2.graph.EdgeAttributes;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.Graph;
import etomica.virial.cluster2.graph.NodeAttributes;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.util.TagsList;

public abstract class AbstractSearchState implements SearchState {

  private Graph firstGraph;
  private Graph secondGraph;
  protected int core_len;
  protected int n1, n2;
  protected int[] core_1, core_2;
  private TagsList tags = new TagsList();

  protected AbstractSearchState() {

    // default constructor for use by descendant classes
  }

  public AbstractSearchState(Graph g1, Graph g2) {

    assert (g1 != null);
    assert (g2 != null);
    firstGraph = g1;
    secondGraph = g2;
    computeTags();
  }

  public void backTrack() {

    // backtracking does nothing by default
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

  public Edges getE1() {

    return getG1().getEdges();
  }

  public EdgeAttributes getE1Attrs(int fromNodeID, int toNodeID) {

    return getE1().getAttributes(fromNodeID, toNodeID);
  }

  public Edges getE2() {

    return getG2().getEdges();
  }

  public EdgeAttributes getE2Attrs(int fromNodeID, int toNodeID) {

    return getE2().getAttributes(fromNodeID, toNodeID);
  }

  public Graph getG1() {

    return firstGraph;
  }

  public Graph getG2() {

    return secondGraph;
  }

  public Nodes getN1() {

    return getG1().getNodes();
  }

  public NodeAttributes getN1Attrs(int nodeID) {

    return getN1().getAttributes(nodeID);
  }

  public Nodes getN2() {

    return getG2().getNodes();
  }

  public NodeAttributes getN2Attrs(int nodeID) {

    return getN2().getAttributes(nodeID);
  }

  public final List<String> getTags() {

    return Collections.unmodifiableList(tags);
  }

  protected void setG1(Graph g1) {

    firstGraph = g1;
  }

  protected void setG2(Graph g2) {

    secondGraph = g2;
  }

  protected abstract String getTag();

  protected void computeTags() {
    
    tags.clear();
    tags.add(getTag());
  }
}