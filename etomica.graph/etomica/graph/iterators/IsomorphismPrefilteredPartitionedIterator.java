package etomica.graph.iterators;

import java.util.Map;

import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.GraphIterator;


public class IsomorphismPrefilteredPartitionedIterator extends PartitionedIterator {

  public IsomorphismPrefilteredPartitionedIterator(Map<Character, Byte> rootMap, Map<Character, Byte> fieldMap) {

    super(rootMap, fieldMap);
  }

  @Override
  public GraphIterator createOuterIterator() {

    return new IsomorphismFilter(super.createOuterIterator());
  }
}