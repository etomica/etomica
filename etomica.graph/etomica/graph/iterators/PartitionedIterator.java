package etomica.graph.iterators;

import java.util.Map;

import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class PartitionedIterator extends CartesianIterator {

  private byte nodeCount;
  private Map<Character, Byte> fieldMap;
  private Map<Character, Byte> rootMap;

  public PartitionedIterator(Map<Character, Byte> rootMap, Map<Character, Byte> fieldMap) {

    this.rootMap = rootMap;
    this.fieldMap = fieldMap;
    this.nodeCount = computeNodeCount();
    bootstrap();
  }

  private byte computeNodeCount() {

    byte result = (byte) 0;
    for (Byte partitionSize : rootMap.values()) {
      result += partitionSize;
    }
    for (Byte partitionSize : fieldMap.values()) {
      result += partitionSize;
    }
    return result;
  }

  @Override
  protected Graph combineGraphs(Graph outer, Graph inner) {

//    System.out.println(outer.getStore().toString() + " x " + inner.nodesToString());
    for (Edge edge : outer.edges()) {
      inner.putEdge(edge.getId());
      inner.getEdge(edge.getId()).setColor(edge.getColor());
    }
    inner.coefficient().multiply(outer.coefficient());
    return inner;
  }

  @Override
  public GraphIterator createInnerIterator() {

    return new EdgelessIterator(rootMap, fieldMap);
  }

  @Override
  public GraphIterator createOuterIterator() {

    return new DefaultIterator(nodeCount);
  }
}