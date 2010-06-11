package etomica.graph.model.impl;

import java.util.HashMap;
import java.util.Map;

import etomica.graph.model.Metadata;

public class MetadataImpl implements Metadata {

  private static final Map<String, Metadata> stock = new HashMap<String, Metadata>();
  private char type;
  private char color;

  // FIXME there has to be a better way to do this
  // For heterogenous systems, turn this on so that root and field points are not
  // interchangeable.
  // even for homogenous fluids, two root points within one component of a graph
  // should also be special since they are rooted at a specific distance from
  // each other.
  public static boolean rootPointsSpecial = false;

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
    return (this == other) || (!rootPointsSpecial && isSameColor(other));
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

  public int compareTo(Metadata other) {

    // it only makes sense to compare objects of the same type
    if (other == null) {
      return 1;
    }
    if (other == this) {
      return 0;
    }
    if (color == other.getColor()) {
      return 0;
    }
    return color > other.getColor() ? 1 : -1;
  }
}