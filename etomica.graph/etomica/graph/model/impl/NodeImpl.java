/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model.impl;

import static etomica.graph.model.Metadata.*;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;

public class NodeImpl implements Node {

  private byte id;
  private Metadata stockMetadata;

  public NodeImpl(byte id, char type, char color) {

    this.id = id;
    stockMetadata = MetadataImpl.getStockComponent(type, color);
  }

  public static Node createFieldNode(byte id, char color) {

    return createNode(id, TYPE_NODE_FIELD, color);
  }

  protected static Node createNode(byte id, char type, char color) {

    assert (type == TYPE_NODE_FIELD || type == TYPE_NODE_ROOT);
    return new NodeImpl(id, type, color);
  }

  public static Node createRootNode(byte id, char color) {

    return createNode(id, TYPE_NODE_ROOT, color);
  }

  public Node copy() {

    return new NodeImpl(getId(), getType(), getColor());
  }

  @Override
  public boolean equals(Object obj) {

    if (obj instanceof Node) {
      Node other = (Node) obj;
      return id == other.getId() && stockMetadata.equals(other.getMetadata());
    }
    return false;
  }

  public char getColor() {

    return stockMetadata.getColor();
  }

  public byte getId() {

    return id;
  }

  public Metadata getMetadata() {

    return stockMetadata;
  }

  public char getType() {

    return stockMetadata.getType();
  }

  public boolean isCompatible(Node other) {

    return stockMetadata.isCompatible(other.getMetadata());
  }

  public boolean isSameColor(Node other) {

    return stockMetadata.isSameColor(other.getMetadata());
  }

  public boolean isSameId(Node other) {

    return (id == other.getId());
  }

  public boolean isSameType(Metadata other) {

    return stockMetadata.isSameType(other);
  }

  public void setColor(char color) {

    stockMetadata = MetadataImpl.getStockComponent(getType(), color);
  }

  public void setType(char type) {

    stockMetadata = MetadataImpl.getStockComponent(type, getColor());
  }

  @Override
  public String toString() {

    return stockMetadata.toString() + id;
  }

  public int compareTo(Node other) {

    if (other == null) {
      return 1;
    }
    if (id != other.getId()) {
      return id > other.getId() ? 1 : -1;
    }
    return stockMetadata.compareTo(other.getMetadata());
  }
}