/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

public class ExpParameters implements Parameters {

  private byte expLo;
  private byte expHi;

  public ExpParameters(byte expLo, byte expHi) {

    this.expLo = expLo;
    this.expHi = expHi;
  }

  public byte expLo() {

    return expLo;
  }

  public byte expHi() {

    return expHi;
  }
}
