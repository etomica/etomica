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