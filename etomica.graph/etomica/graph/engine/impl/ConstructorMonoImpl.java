/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import etomica.graph.engine.Parser.ConstructorMono;

public class ConstructorMonoImpl extends ConstructorImpl implements ConstructorMono {

  private byte fieldNodes;
  private boolean isoFree;
  private byte rootNodes;

  public ConstructorMonoImpl(byte rootNodes, byte fieldNodes, boolean isoFree) {

    this.rootNodes = rootNodes;
    this.fieldNodes = fieldNodes;
    this.isoFree = isoFree;
  }

  public byte getFieldNodes() {

    return fieldNodes;
  }

  public byte getRootNodes() {

    return rootNodes;
  }

  public boolean isIsoFree() {

    return isoFree;
  }
}