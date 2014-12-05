/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesSetVisitor;
import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.graph.Nodes;
import etomica.virial.cluster2.util.TagsList;

public class SimpleGraphSet implements GraphSet {

  private Nodes nodes;
  private List<Edges> edgesList;
  private TagsList tags;

  public SimpleGraphSet(Nodes nodes, List<Edges> edgesList) {

    this.nodes = nodes;
    this.edgesList = edgesList;
  }

  public void addComplements() {

    List<Edges> complements = new ArrayList<Edges>(edgesList.size());
    for (int i = 0; i < edgesList.size(); i++) {
      complements.add(edgesList.get(i).complement());
    }
    edgesList.addAll(complements);
  }

  public Set<Edges> getEdgesSet() {

    return new HashSet<Edges>(edgesList);
  }

  public int getSize() {

    return edgesList.size();
  }

  public Nodes getNodes() {

    return nodes;
  }

  @Override
  public String toString() {

    String result = "";
    for (int i = 0; i < getSize(); i++) {
      result += i + ": " + edgesList.get(i).toString() + "\n";
    }
    return result;
  }

  public void visitEdgesSet(EdgesSetVisitor visitor) {

    for (int i = 0; i < getSize(); i++) {
      if (!visitor.visit(edgesList.get(i))) {
        break;
      }
    }
  }

  public void setTags(List<String> tags) {

    this.tags = new TagsList(tags);
  }

  public TagsList getTags() {

    return tags;
  }
}