/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.Value;

public class ValueImpl implements Value {

  private String value;

  public ValueImpl(String value) {

    this.value = value;
  }

  public String getValue() {

    return this.value;
  }
}