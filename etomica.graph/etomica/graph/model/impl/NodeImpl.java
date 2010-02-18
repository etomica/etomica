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

    return stockMetadata.isCompatible(other.getMetadata())
        && (getType() == TYPE_NODE_FIELD || (getType() == TYPE_NODE_ROOT && isSameId(other)));
  }

  public boolean isSameColor(Metadata other) {

    return stockMetadata.isSameColor(other);
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
}