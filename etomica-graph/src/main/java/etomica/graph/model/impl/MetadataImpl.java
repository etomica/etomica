/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.impl;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import etomica.graph.model.Metadata;

public class MetadataImpl implements Metadata {

  private static final Map<String, Metadata> stock = new HashMap<String, Metadata>();
  public static Comparator<Metadata> metaDataComparator = null;
  private char type;
  private char color;
  public static final List<ArrayList<Character>> edgeColorPairs = new ArrayList<ArrayList<Character>>();

  // FIXME there has to be a better way to do this
  // For heterogenous systems, turn this on so that root and field points are not
  // interchangeable.
  // even for homogenous fluids, two root points within one component of a graph
  // should also be special since they are rooted at a specific distance from
  // each other.
  public static boolean rootPointsSpecial = false;//only true for Wertheim diagram generator
  
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

    return (other != null && color == other.getColor());
  }

  public boolean isSameType(Metadata other) {

    return (other != null && type == other.getType());
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
    if (color == other.getColor() && type == other.getType()) {
      return 0;
    }
    if (metaDataComparator != null) {
      return metaDataComparator.compare(this, other);
    }
    if (color != other.getColor()) {
      return color > other.getColor() ? 1 : -1;
    }
    return type > other.getType() ? 1 : -1;
  }
  
  public static char getReverseEdgeColor(char color){//change color of bond to reverse color
	  for(List<Character> list:edgeColorPairs ){
		  char c0 = list.get(0);
		  char c1 = list.get(1);
		  if (color == c0){
			  return c1;
		  }
		  else if (color == c1){
			  return c0;
		  }
	  }
	  return  color;
  }
}