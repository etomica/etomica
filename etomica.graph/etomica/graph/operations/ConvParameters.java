/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import etomica.graph.operations.Mul.MulParameters;

public class ConvParameters implements Parameters {

  private byte nodeId;
  private MulParameters mulParameters;

  public ConvParameters(byte nodeId, MulParameters mulParameters) {
    this.nodeId = nodeId;
    this.mulParameters = mulParameters;
  }

  public byte nodeId() {
    return nodeId;
  }
  
  public MulParameters mulParameters() {
    return mulParameters;
  }
}
