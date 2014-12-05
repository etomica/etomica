/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

public interface GraphCoefficient {

  public int getValue1();

  public int getValue2();

  public void setValue1(int newValue);

  public void setValue2(int newValue);
  
  public int getSign();

  public GraphCoefficient switchSign();

  public void inc();

  public GraphCoefficient add(GraphCoefficient value);
}
