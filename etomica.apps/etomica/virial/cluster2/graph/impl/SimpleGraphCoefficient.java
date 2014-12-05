/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import etomica.virial.cluster2.graph.GraphCoefficient;

public class SimpleGraphCoefficient implements GraphCoefficient {

  private int value1;
  private int value2;
  private int sign;

  public SimpleGraphCoefficient(int v) {
    
    this(v,1,1);
  }
  
  public SimpleGraphCoefficient(int v1, int v2) {
    
    this(v1,v2,1);
  }
  
  public SimpleGraphCoefficient(int v1, int v2, int sgn) {
    
    value1 = v1;
    value2 = v2;
    sign = sgn;
  }

  public void inc() {
  
    value1 += value2;
  }
  
  public int getSign() {

    return sign;
  }

  public int getValue1() {

    return value1;
  }

  public int getValue2() {

    return value2;
  }

  public void setValue1(int newValue) {

    value1 = newValue;
  }

  public void setValue2(int newValue) {

    value2 = newValue;
  }

  public GraphCoefficient switchSign() {

    return new SimpleGraphCoefficient(getValue1(), getValue2(), -getSign());
  }
  
  public GraphCoefficient add(GraphCoefficient value) {

    return this;
  }
}