/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

import static etomica.graph.model.Metadata.*;

import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.model.impl.EdgeImpl;
import etomica.graph.model.impl.GraphImpl;
import etomica.graph.model.impl.NodeImpl;

public class GraphFactory {

  public static Coefficient createCoefficient() {

    return new CoefficientImpl(1);
  }

  public static Edge createEdge(byte edgeId) {

    return EdgeImpl.createEdge(edgeId, COLOR_CODE_DEFAULT);
  }

  public static Node createNode(byte nodeId) {

    return NodeImpl.createFieldNode(nodeId, COLOR_CODE_DEFAULT);
  }

  public static Node createNode(byte nodeId, boolean isRootNode) {

    if (isRootNode) {
      return NodeImpl.createRootNode(nodeId, COLOR_CODE_DEFAULT);
    }
    return NodeImpl.createFieldNode(nodeId, COLOR_CODE_DEFAULT);
  }

  public static Node createNode(byte nodeId, char color, char type) {

    if (type == TYPE_NODE_FIELD) {
      return NodeImpl.createFieldNode(nodeId, color);
    }
    else if (type == TYPE_NODE_ROOT) {
      return NodeImpl.createRootNode(nodeId, color);
    }
    return null;
  }

  public static Graph createGraph(byte nodeCount, Bitmap store) {

    return new GraphImpl(nodeCount, store);
  }

  public static Graph createGraph(byte nodeCount, byte rootNodeCount, Bitmap store) {

    return new GraphImpl(nodeCount, rootNodeCount, store);
  }

  public static Graph createGraph(byte nodeCount) {

    return new GraphImpl(nodeCount);
  }

  public static Graph createGraph(Node[] nodes) {

    return new GraphImpl(nodes);
  }
}