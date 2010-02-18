package etomica.graph.model.impl;

import etomica.graph.model.Edge;
import etomica.graph.model.Metadata;
import static etomica.graph.model.Metadata.*;

public class EdgeImpl implements Edge {

  private byte id;
  private Metadata stockMetadata;

  public EdgeImpl(byte id, char color) {

    this.id = id;
    stockMetadata = MetadataImpl.getStockComponent(TYPE_EDGE_ANY, color);
  }

  public static Edge createEdge(byte id, char color) {

    return new EdgeImpl(id, color);
  }

  public Edge copy() {

    return new EdgeImpl(getId(), getColor());
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

  public boolean isCompatible(Edge other) {

    return stockMetadata.isCompatible(other.getMetadata());
  }

  public boolean isSameColor(Edge other) {

    return stockMetadata.isSameColor(other.getMetadata());
  }

  public boolean isSameId(Edge other) {

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

    return stockMetadata.toString();
  }
}