/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import etomica.virial.cluster2.graph.NodeAttributes;
import etomica.virial.cluster2.graph.Nodes;

/**
 * Because a new constructor is needed for this Nodes class, a nested class
 * had to be defined instead of using an anonymous class, which does not
 * support non-default constructors. This class supports both field and root
 * nodes. Field nodes are compatible with field nodes, and root nodes are
 * incompatible with all nodes. We assume that the first fieldNodeCount nodes
 * are field nodes, and the remaining are root nodes.
 *
 * @author Demian Lessa
 */
public class SimpleNodes implements Nodes {

  private NodeAttributes[] fieldAttributes;
  private NodeAttributes[] rootAttributes;
  private Map<Character,List<Integer>> fieldPartitions;
  private Map<Character,List<Integer>> rootPartitions;

  /*
   * The parameters are vectors of colors, one for each field node
   * and one for each root node
   */
  public SimpleNodes(final char[] fieldColors, final char[] rootColors) {

    assert (fieldColors != null);
    assert (rootColors != null);
    assert (fieldColors.length + rootColors.length >= 0);
    assert (fieldColors.length + rootColors.length <= 255);
    configure(fieldColors, rootColors);
  }

  public SimpleNodes(byte fieldNodeCount, byte rootNodeCount) {

    assert (fieldNodeCount >= 0);
    assert (rootNodeCount >= 0);
    assert (fieldNodeCount + rootNodeCount >= 0);
    assert (fieldNodeCount + rootNodeCount <= 255);
    char[] fieldColors = new char[fieldNodeCount];
    char[] rootColors = new char[rootNodeCount];
    for (int i = 0; i < fieldNodeCount; i++) {
      fieldColors[i] = Nodes.NODE_COLOR_DEFAULT;
    }
    for (int i = 0; i < rootNodeCount; i++) {
      rootColors[i] = Nodes.NODE_COLOR_DEFAULT;
    }
    configure(fieldColors, rootColors);
  }

  public static NodeAttributes fieldNodeAttributes(char color) {

    return new FieldNodeAttributes(color);
  }

  public static NodeAttributes rootNodeAttributes(char color) {

    return new RootNodeAttributes(color);
  }

  private void configure(final char[] fieldColors, final char[] rootColors) {

    rootAttributes = new NodeAttributes[rootColors.length];
    rootPartitions = new HashMap<Character,List<Integer>>();
    for (int i = 0; i < rootAttributes.length; i++) {
      rootAttributes[i] = new RootNodeAttributes(rootColors[i]);
      if (rootPartitions.containsKey(rootColors[i])) {
        rootPartitions.get(rootColors[i]).add(i);
      }
      else {
        List<Integer> colorsIndices = new ArrayList<Integer>();
        colorsIndices.add(i);
        rootPartitions.put(rootColors[i], colorsIndices);
      }
    }
    fieldAttributes = new NodeAttributes[fieldColors.length];
    fieldPartitions = new HashMap<Character,List<Integer>>();
    for (int i = 0; i < fieldAttributes.length; i++) {
      fieldAttributes[i] = new FieldNodeAttributes(fieldColors[i]);
      if (fieldPartitions.containsKey(fieldColors[i])) {
        fieldPartitions.get(fieldColors[i]).add(rootAttributes.length + i);
      }
      else {
        List<Integer> colorsIndices = new ArrayList<Integer>();
        colorsIndices.add(rootAttributes.length + i);
        fieldPartitions.put(fieldColors[i], colorsIndices);
      }
    }
  }

  public byte count() {

    return (byte) (fieldAttributes.length + rootAttributes.length);
  }

  // from 0..rootNodeCount-1 => root node
  // from rootNodeCount..count() => field node
  public NodeAttributes getAttributes(int nodeID) {

    assert (nodeID >= 0 && nodeID < count());
    if (nodeID < rootAttributes.length) {
      return rootAttributes[nodeID];
    }
    else {
      return fieldAttributes[nodeID - rootAttributes.length];
    }
  }

  public List<Integer> getPartition(byte nodeClass, char nodeColor) {

    assert (nodeClass == NODE_CLASS_FIELD || nodeClass == NODE_CLASS_ROOT);
    if (nodeClass == NODE_CLASS_FIELD) {
      assert(fieldPartitions.containsKey(nodeColor));
      return fieldPartitions.get(nodeColor);
    }
    assert(rootPartitions.containsKey(nodeColor));
    return rootPartitions.get(nodeColor);
  }

  public Set<Character> getColors(byte nodeClass) {

    assert (nodeClass == NODE_CLASS_FIELD || nodeClass == NODE_CLASS_ROOT);
    if (nodeClass == NODE_CLASS_FIELD) {
      return Collections.unmodifiableSet(fieldPartitions.keySet());
    }
     return Collections.unmodifiableSet(rootPartitions.keySet());
  }
}

abstract class AbstractNodeAttributes implements NodeAttributes {

  private char colorID;

  public AbstractNodeAttributes(char color) {

    colorID = color;
  }

  public char getColor() {

    return colorID;
  }

  public boolean isCompatible(NodeAttributes attr) {

    return isSameClass(attr) && isSameColor(attr);
  }

  public boolean isSameColor(NodeAttributes attr) {

    return attr.getColor() == colorID;
  }
}

class FieldNodeAttributes extends AbstractNodeAttributes {

  public FieldNodeAttributes(char colorID) {

    super(colorID);
  }

  public boolean isSameClass(NodeAttributes attr) {

    return (attr instanceof FieldNodeAttributes);
  }
}

class RootNodeAttributes extends AbstractNodeAttributes {

  public RootNodeAttributes(char colorID) {

    super(colorID);
  }

  public boolean isCompatible(NodeAttributes attr) {

    return equals(attr);
  }

  public boolean isSameClass(NodeAttributes attr) {

    return (attr instanceof RootNodeAttributes);
  }
}
//static class SimpleNodes implements Nodes {
//
//  private byte fieldNodeCount = 0;
//  private byte rootNodeCount = 0;
//  private NodeAttributes[] rootNodeAttributes;
//  private List<List<Integer>> partition;
//
//  private NodeAttributes[] fieldAttributes;
//  private NodeAttributes[] rootAttributes;
//  private Map<Integer,List<Integer>> fieldPartitions;
//
//  /*
//   * The parameters are vectors of colors, one for each field node
//   * and one for each root node
//   */
//  public SimpleNodes(final byte[] fieldColors, final byte[] rootColors) {
//
//    assert (fieldColors != null);
//    assert (rootColors != null);
//    assert (fieldColors.length + rootColors.length > 0);
//    assert (fieldColors.length + rootColors.length <= 255);
//    fieldAttributes = new NodeAttributes[fieldColors.length];
//    rootAttributes = new NodeAttributes[rootColors.length];
//    for (int i = 0; i < fieldAttributes.length; i++) {
//      fieldAttributes[i] = new FieldNodeAttributes(fieldColors[i]);
//      if (fieldPartitions.containsKey(fieldColors[i])) {
//        fieldPartitions.get(fieldColors[i]).add(i);
//      }
//      else {
//        List<Integer> colorsIndices = new ArrayList<Integer>();
//        colorsIndices.add(i);
//        fieldPartitions.put(new Integer(fieldColors[i]), colorsIndices);
//      }
//    }
//    for (int i = 0; i < rootAttributes.length; i++) {
//      rootAttributes[i] = new RootNodeAttributes(rootColors[i]);
//    }
//  }
//
//  public SimpleNodes(byte fieldNodeCount, byte rootNodeCount) {
//
//    assert (fieldNodeCount >= 0);
//    assert (rootNodeCount >= 0);
//    assert (fieldNodeCount + rootNodeCount > 0);
//    assert (fieldNodeCount + rootNodeCount <= 255);
//    this.fieldNodeCount = fieldNodeCount;
//    this.rootNodeCount = rootNodeCount;
//    this.rootNodeAttributes = new RootNodeAttributes[rootNodeCount];
//    for (byte i = 0; i < rootNodeCount; i++) {
//      rootNodeAttributes[i] = new RootNodeAttributes(i);
//    }
//    int partitionCount = rootNodeCount + (fieldNodeCount > 0 ? 1 : 0);
//    List<List<Integer>> ps = new ArrayList<List<Integer>>(partitionCount);
//    for (int i = 0; i < rootNodeCount; i++) {
//      ArrayList<Integer> p = new ArrayList<Integer>(1);
//      p.add(i);
//      ps.add(Collections.unmodifiableList(p));
//    }
//    if (fieldNodeCount > 0) {
//      ArrayList<Integer> p = new ArrayList<Integer>(fieldNodeCount);
//      for (int i = 0; i < fieldNodeCount; i++) {
//        p.add(rootNodeCount + i);
//      }
//      ps.add(Collections.unmodifiableList(p));
//    }
//    partition = Collections.unmodifiableList(ps);
//  }
//
//  public byte count() {
//
//    return (byte) (fieldNodeCount + rootNodeCount);
//  }
//
//  // from 0..rootNodeCount-1 => root node
//  // from rootNodeCount..count() => field node
//  public NodeAttributes getAttributes(int nodeID) {
//
//    assert (nodeID >= 0 && nodeID < count());
//    if (nodeID < rootNodeCount) {
//      return rootNodeAttributes[nodeID];
//    }
//    else {
//      return FIELD_NODE_ATTRIBUTES;
//    }
//  }
//
//  public List<Integer> getPartition(int partitionID) {
//
//    assert (partitionID >= 0 && partitionID < getPartitionCount());
//    return partition.get(partitionID);
//  }
//
//  public int getPartitionCount() {
//
//    return partition.size();
//  }
//}