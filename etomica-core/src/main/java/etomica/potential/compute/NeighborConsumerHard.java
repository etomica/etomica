/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.IAtom;
import etomica.space.Vector;

/**
 * Interface for neighbor iteration callback for hard MD (includes state for
 * the relevant atom pair).
 */
public interface NeighborConsumerHard extends NeighborIterator.NeighborConsumer {
    void acceptHard(IAtom jAtom, Vector rij, int state);
}
