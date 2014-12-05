/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

public class DeleteEdgeParameters implements Parameters {

  private char color;

  public DeleteEdgeParameters(char color) {

    this.color = color;
  }

  public char color() {

    return color;
  }
}
