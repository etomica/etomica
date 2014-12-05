/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import java.util.Map;

import etomica.graph.engine.Parser.Expression;
import etomica.graph.engine.Parser.UnaryOpByteMap;

public class UnaryOpByteMapImpl extends UnaryOpImpl implements UnaryOpByteMap {

  private Map<Byte, Byte> colorMap;

  public UnaryOpByteMapImpl(String operation, Expression expression, Map<Byte, Byte> colorMap) {

    super(operation, expression);
    this.colorMap = colorMap;
  }

  public Map<Byte, Byte> getColorMap() {

    return this.colorMap;
  }
}
