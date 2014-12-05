/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import etomica.virial.cluster2.graph.impl.SimpleNodes;


public class GraphSetFactory {
  
  public static GraphSet completeGraphSet(char[] fieldColors, char[] rootColors) {

    return completeGraphSet(fieldColors, rootColors, false);

  }

  public static GraphSet completeGraphSet(char[] fieldColors, char[] rootColors, boolean filterIsomorphs) {

    Nodes nodes = new SimpleNodes(fieldColors, rootColors);
    return GraphFactory.naiveGraphSet(nodes, null, filterIsomorphs);

  }

  public static GraphSet isomorphFreeGraphSet(char[] fieldColors, char[] rootColors) {

    Nodes nodes = new SimpleNodes(fieldColors, rootColors);
    return GraphFactory.storedGraphSet(nodes, null, true);
  }
}