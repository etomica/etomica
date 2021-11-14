/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.space.Space;

public class PotentialNonAdditiveDifference extends PotentialMolecular {

    protected final IPotentialMolecular p2, pFull;
    protected final MoleculePair pair;
    
    public PotentialNonAdditiveDifference(Space space, IPotentialMolecular p2, IPotentialMolecular pFull) {
        super(Integer.MAX_VALUE, space);
        this.p2 = p2;
        this.pFull = pFull;
        pair = new MoleculePair();
    }

    public double energy(IMoleculeList molecules) {
        double u = pFull.energy(molecules);
        for (int i = 0; i<molecules.size(); i++) {
            pair.mol0 = molecules.get(i);
            for (int j = i+1; j<molecules.size(); j++) {
                pair.mol1 = molecules.get(j);
                u -= p2.energy(pair);
            }
        }
        return u;
    }
}
