/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import java.util.Map;

import etomica.graph.engine.Parser.ConstructorColored;

public class ConstructorColoredImpl extends ConstructorImpl implements ConstructorColored {

  private Map<Character, Byte> fieldColorMap;
  private boolean isoFree;
  private Map<Character, Byte> rootColorMap;

  public ConstructorColoredImpl(Map<Character, Byte> rootColorMap, Map<Character, Byte> fieldColorMap,
      boolean isoFree) {

    this.rootColorMap = rootColorMap;
    this.fieldColorMap = fieldColorMap;
    this.isoFree = isoFree;
  }

  public Map<Character, Byte> getFieldColorMap() {

    return fieldColorMap;
  }

  public Map<Character, Byte> getRootColorMap() {

    return rootColorMap;
  }

  public boolean isIsoFree() {

    return isoFree;
  }
}
