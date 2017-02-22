/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

public class IntParameters implements Parameters {

  private byte nodeId;

  public IntParameters(byte nodeId) {

    this.nodeId = nodeId;
  }

  public byte nodeId() {

    return nodeId;
  }
}
