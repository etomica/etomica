package etomica.virial.cluster2.graph.impl;

import java.util.HashMap;
import java.util.Map;

import etomica.virial.cluster2.graph.EdgeAttributes;

public class SimpleEdgeAttributes implements EdgeAttributes {

  private char colorID;

  private SimpleEdgeAttributes(char color) {
    
    colorID = color;
  }
  
  private static Map<Character, EdgeAttributes> stock = new HashMap<Character, EdgeAttributes>();
  
  public static EdgeAttributes getAttributes(char color) {
    
    if (!stock.containsKey(color)) {
      stock.put(color, new SimpleEdgeAttributes(color));
    }
    return stock.get(color);
  }
  
  public boolean isCompatible(final EdgeAttributes attr) {

    return isSameColor(attr);
  }

  public char getColor() {

    return colorID;
  }

  public boolean isSameColor(EdgeAttributes attr) {

    return attr.getColor() == getColor();
  }
}