/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

public class IteratorToSet {

  public Set<Graph> getSet(GraphIterator iterator) {

    Set<Graph> result = new HashSet<Graph>();
    while (iterator.hasNext()) {
      result.add(iterator.next());
    }
    return result;
  }
}