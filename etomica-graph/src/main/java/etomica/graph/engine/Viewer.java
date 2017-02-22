/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import java.util.Set;

import etomica.graph.model.Graph;

public interface Viewer {

  public void close();

  public void open();

  public void update(Set<Graph> graphs);
}