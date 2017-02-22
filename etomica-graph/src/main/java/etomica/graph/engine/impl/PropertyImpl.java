/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Property;

public class PropertyImpl implements Property {

  private String property;

  public PropertyImpl(String property) {

    this.property = property;
  }

  public String getName() {

    return this.property;
  }
}
