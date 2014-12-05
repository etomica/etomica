/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

public interface Console {

  public void clear();

  public void setReader(ConsoleReader reader);

  public void writeLn();

  public void write(Exception e);

  public void write(String value);

  public void updateBegin();

  public void updateDone();
}