/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.CommandPropertyValue;
import etomica.graph.engine.Parser.Property;
import etomica.graph.engine.Parser.Value;

public class CommandPropertyValueImpl extends CommandImpl implements CommandPropertyValue {

  private Property property;
  private Value value;

  public CommandPropertyValueImpl(String command, Property property, Value value) {

    super(command);
    this.property = property;
    this.value = value;
  }

  public Property getProperty() {

    return this.property;
  }

  public Value getValue() {

    return this.value;
  }
}
