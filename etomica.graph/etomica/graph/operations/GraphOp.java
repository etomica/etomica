/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import etomica.graph.model.Graph;

/**
 * Interface for a simple operation that takes a Graph and returns a Graph.
 * 
 * @author Andrew Schultz
 */
public interface GraphOp {
  public Graph apply(Graph g);

  /**
   * A dummy implementation of GraphOp that does nothing.
   * 
   * @author Andrew Schultz
   */
  public static class GraphOpNull implements GraphOp {
    public Graph apply(Graph g) {return g;}
  }
}