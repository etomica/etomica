/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IPotentialAtomic;
import etomica.atom.AtomType;


/**
 * Interface for a potential that is artificially truncated, and thus which can
 * provide a zero-body potential that (approximately) corrects for the
 * truncation. The PotentialMaster gets a LRC potential from a potential
 * implementing this interface when the potential is given to the
 * PotentialMaster's setSpecies method.  PotentialGroup also recognizes
 * this interface and adds a LRC potential to the PotentialMaster when a
 * PotentialTruncated is added to it via the type-specifying addPotential method.
 * 
 * @see PotentialMaster
 * @see PotentialGroup
 */
public interface PotentialTruncated extends IPotentialAtomic {

    /**
     * Returns a class that calculates the long-range contribution to the potential
     * that becomes neglected by the truncation.  May return null.
     */
    Potential0Lrc makeLrcPotential(AtomType[] types);
}
