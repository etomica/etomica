/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import etomica.graph.iterators.IteratorWrapper;
import etomica.graph.iterators.filters.GlobalFilter;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;

import java.util.HashSet;
import java.util.Set;

public class IsoFree implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {

    IteratorWrapper wrapper = new IteratorWrapper(argument.iterator(), true);
    GraphIterator isomorphs = new IsomorphismFilter(wrapper, null);
    Set<Graph> result = new HashSet<Graph>();
    while (isomorphs.hasNext()) {
      result.add(isomorphs.next());
    }
    return result;
  }

  public static class IsoFreeParams extends GlobalFilter.SignatureMaker implements Parameters {
  }
}