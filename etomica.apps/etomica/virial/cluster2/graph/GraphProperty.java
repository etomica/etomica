/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

/**
 * This interface generalizes algorithms that decide whether a graph G has a
 * particular property P. The property checker algorithm only needs to look at
 * the graph nodes and vertices to decide whether G has property P. An example
 * of an implementing class is a class that checks if G is connected.
 */
public interface GraphProperty {

  public boolean check(Nodes nodes, Edges edges);
}