/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import java.util.HashSet;
import java.util.Set;

import etomica.graph.engine.Parser.Constructor;

public class ConstructorImpl implements Constructor {

  private Set<String> filters = new HashSet<String>();

  public void addFilter(String filterName) {

    this.filters.add(filterName);
  }

  public boolean hasFilter(String filterName) {

    return filters.contains(filterName);
  }
}