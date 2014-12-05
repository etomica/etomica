/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.Nodes;

public abstract class AbstractNodes implements Nodes {

  private byte nodeCount = 0;

  public AbstractNodes(byte numNodes) {

    if (numNodes == 0) {
      throw new RuntimeException("Invalid node count (N=0)");
    }
    nodeCount = numNodes;
  }

  public byte count() {

    return nodeCount;
  }
}