/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

public class PowParameters implements Parameters {

  private byte exponent;

  public PowParameters(byte exponent) {

    this.exponent = exponent;
  }

  public byte exponent() {

    return exponent;
  }
}
