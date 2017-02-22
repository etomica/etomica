/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

public final class RelabelParameters implements Parameters {

  private byte[] permutation;

  public RelabelParameters(byte[] permutation) {

    this.permutation = permutation;
  }

  public byte map(byte index) {

    return permutation[index];
  }

  public int size() {

    return permutation.length;
  }
}