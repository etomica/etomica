/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators.filters;


import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.Property;

public class PropertyFilter extends LocalFilter {

  private Property property;

  public PropertyFilter(GraphIterator iterator, Property property) {

    super(iterator);
    this.property = property;
  }

  @Override
  protected boolean accept(Graph g) {

    return property.check(g);
  }
}