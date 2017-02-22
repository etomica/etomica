/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.tmmc;

import etomica.api.IBox;

/**
 * Interface for class that defines the macrostates.  Provides
 * methods that tell the number of expected macrostates in the
 * box, and that give an index indicating which macrostate
 * the box presently is in.
 */
public interface MacrostateManager {
 /**
  * Number of expected macrostates in the box.
  */
    public int numberOfStates(IBox p);
 /**
  * Returns an index identifying the current macrostate occupied by the box.
  */
    public int stateIndex(IBox p);
    
 /**
  * Returns the value of the macrostate variable corresponding to the given index.
  */
    public double state(int i);
}//end of MacrostateManager
