/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom.iterator;


import etomica.potential.IteratorDirective;

/**
 * 
 * Interface for an iterator that can interpret specification of
 * direction UP or DOWN.
 *
 * @author David Kofke
 *
 */
public interface AtomsetIteratorDirectable extends AtomsetIterator {

	public void setDirection(IteratorDirective.Direction direction);
	
}
