/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph;

import java.util.List;

public interface EdgesFilter extends Tagged {

  // reserved for simple algorithms- O(1), O(N), etc
  public boolean preAccept(Edges edges, List<Edges> edgesList);

  // reserved for more complex algorithms
  public boolean accept(Edges e, List<Edges> edgesList);

  // chain a filter
  public void chain(EdgesFilter filter);
}