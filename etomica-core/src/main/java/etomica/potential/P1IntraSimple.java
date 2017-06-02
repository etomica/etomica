/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.iterator.ApiBuilder;


/**
 * Generic intramolecular potential group, having one potential for bonded
 * atoms, and a different potential for unbonded ones.
 *
 * @author David Kofke
 */
public class P1IntraSimple {
    
    public static PotentialGroup makeP1IntraSimple(PotentialMaster potentialMaster, Potential2 bonded, Potential2 nonbonded) {
        PotentialGroup pGroup = potentialMaster.makePotentialGroup(1);
        pGroup.addPotential(bonded, ApiBuilder.makeAdjacentPairIterator());
        pGroup.addPotential(nonbonded, ApiBuilder.makeNonAdjacentPairIterator());
        return pGroup;
    }
}
   
