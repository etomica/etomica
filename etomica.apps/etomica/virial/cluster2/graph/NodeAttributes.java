/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

public interface NodeAttributes {

  public char getColor();

  public boolean isCompatible(NodeAttributes attr);

  public boolean isSameClass(NodeAttributes attr);

  public boolean isSameColor(NodeAttributes attr);
}