/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import etomica.graph.property.Property;

public final class SplitParameters implements Parameters {

  private char edgeColor;
  private char newColor0;
  private char newColor1;
  private final Property discardCriteria;

  public SplitParameters(char edgeColor, char newColor0, char newColor1) {
    this(edgeColor, newColor0, newColor1, null);
  }
  
  public SplitParameters(char edgeColor, char newColor0, char newColor1, Property discardCriteria) {
    this.edgeColor = edgeColor;
    this.newColor0 = newColor0;
    this.newColor1 = newColor1;
    this.discardCriteria = discardCriteria;
  }

  public char edgeColor() {

    return edgeColor;
  }

  public char newColor0() {

    return newColor0;
  }

  public char newColor1() {

    return newColor1;
  }
  
  public Property getDiscardProperty() {
    return discardCriteria;
  }
}