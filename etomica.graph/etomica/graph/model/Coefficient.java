/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.model;

public interface Coefficient {

  public void add(Coefficient value);

  public Coefficient copy();

  public int getDenominator();

  public int getNumerator();

  public double getValue();

  public void multiply(Coefficient value);

  public void divide(Coefficient value);

  public void setDenominator(int value);

  public void setNumerator(int value);
  
  public boolean hasOverflow();
}
