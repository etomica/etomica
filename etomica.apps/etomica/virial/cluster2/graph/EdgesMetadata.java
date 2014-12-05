/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

public interface EdgesMetadata {

  public EdgesMetadata copy();

  public EdgesMetadata ncopy();

  public GraphCoefficient getCoefficient();

  public void setCoefficient(GraphCoefficient value);
}