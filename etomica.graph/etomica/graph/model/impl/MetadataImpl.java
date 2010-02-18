package etomica.graph.model.impl;

import java.util.HashMap;
import java.util.Map;

import etomica.graph.model.Metadata;

public class MetadataImpl implements Metadata {

  private static final Map<String, Metadata> stock = new HashMap<String, Metadata>();
  private char type;
  private char color;

  protected MetadataImpl(char type, char color) {

    this.type = type;
    this.color = color;
  }

  public static Metadata getStockComponent(char type, char color) {

    String key = "" + type + color;
    if (!stock.containsKey(key)) {
      stock.put(key, new MetadataImpl(type, color));
    }
    return stock.get(key);
  }

  public char getColor() {

    return color;
  }

  public char getType() {

    return type;
  }

  public boolean isCompatible(Metadata other) {

    return (this == other);
  }

  public boolean isSameColor(Metadata other) {

    return (other != null && getColor() == other.getColor());
  }

  public boolean isSameType(Metadata other) {

    return (other != null && getType() == other.getType());
  }

  @Override
  public String toString() {

    return "" + type + color;
  }
}