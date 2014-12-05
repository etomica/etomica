/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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