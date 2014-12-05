/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import static etomica.graph.model.Metadata.*;

import java.util.Map;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.Node;

public class EdgelessIterator implements GraphIterator {

  private byte rootNodeCount;
  private char[] rootColors;
  private int[] rootPartitions;
  private byte fieldNodeCount;
  private char[] fieldColors;
  private int[] fieldPartitions;
  private CartesianPermutator permutator;

  public EdgelessIterator(Map<Character, Byte> rootMap, Map<Character, Byte> fieldMap) {

    // root partitions
    this.rootNodeCount = 0;
    this.rootColors = new char[rootMap.size()];
    this.rootPartitions = new int[rootMap.size()];
    int prindex = 0;
    for (Character color : rootMap.keySet()) {
      rootNodeCount += rootMap.get(color);
      rootPartitions[prindex] = rootMap.get(color);
      rootColors[prindex] = color;
      prindex++;
    }
    // field partitions
    this.fieldNodeCount = 0;
    this.fieldColors = new char[fieldMap.size()];
    this.fieldPartitions = new int[fieldMap.size()];
    int pfindex = 0;
    for (Character color : fieldMap.keySet()) {
      fieldNodeCount += fieldMap.get(color);
      fieldPartitions[pfindex] = fieldMap.get(color);
      fieldColors[pfindex] = color;
      pfindex++;
    }

    this.permutator = new CartesianPermutator(rootPartitions, fieldPartitions);
  }

  public boolean hasNext() {

    return permutator.hasNext();
  }

  public Graph next() {

    if (hasNext()) {
      return permutationToGraph(permutator.next());
    }
    return null;
  }

  private Graph permutationToGraph(byte[] permutation) {

    // System.out.print(Arrays.toString(permutation) + " :: ");
    Node[] nodes = new Node[permutation.length];
    for (byte nodeId = 0; nodeId < nodes.length; nodeId++) {
      if (nodeId < rootNodeCount) {
        nodes[nodeId] = GraphFactory.createNode(nodeId, rootColors[permutation[nodeId]], TYPE_NODE_ROOT);
      }
      else {
        nodes[nodeId] = GraphFactory.createNode(nodeId, fieldColors[permutation[nodeId]], TYPE_NODE_FIELD);
      }
    }
    return GraphFactory.createGraph(nodes);
  }

  public void remove() {

    // no-op
  }
}