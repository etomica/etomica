/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Color;

public class ColorEntry {

  private Color color;
  private String text;

  public ColorEntry(Color color, String text) {

    this.color = color;
    this.text = text;
  }

  public String getText() {

    return text;
  }

  public Color getColor() {

    return color;
  }
}