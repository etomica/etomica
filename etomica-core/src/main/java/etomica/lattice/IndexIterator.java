/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.lattice;


/**
 * Loops through the set of indexes appropriate to a lattice of given size.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on Jan 3, 2005 by kofke
 */
public interface IndexIterator {

    public void reset();
    
    public boolean hasNext();
    
    public int[] next();
    
    public int getD();
}
