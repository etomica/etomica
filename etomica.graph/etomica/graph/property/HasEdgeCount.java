/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;

public class HasEdgeCount implements Property {

  private byte edgeCount;

  public HasEdgeCount(byte edgeCount) {

    this.edgeCount = edgeCount;
  }

  public boolean check(Graph graph) {

    return graph.edgeCount() == edgeCount;
  }
}
