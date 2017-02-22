/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

public interface Node extends Comparable<Node> {

  public Node copy();

  public char getColor();

  public byte getId();

  public Metadata getMetadata();

  public char getType();

  public boolean isCompatible(Node other);

  public boolean isSameColor(Node other);

  public boolean isSameId(Node other);

  public void setColor(char color);

  public void setType(char type);
}