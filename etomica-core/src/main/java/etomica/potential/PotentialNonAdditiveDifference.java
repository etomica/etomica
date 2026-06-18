/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;

public class PotentialNonAdditiveDifference implements IPotentialMolecular {

    protected final IPotentialMolecular p2, pFull;
    protected final MoleculeArrayList pair;
    
    public PotentialNonAdditiveDifference(IPotentialMolecular p2, IPotentialMolecular pFull) {
        super();
        this.p2 = p2;
        this.pFull = pFull;
        pair = new MoleculeArrayList(2);
    }

    public double energy(IMoleculeList molecules) {
        double u = pFull.energy(molecules);

        for (int i = 0; i<molecules.size(); i++) {
            pair.set(0, molecules.get(i));
            for (int j = i+1; j<molecules.size(); j++) {
                pair.set(1, molecules.get(j));
                u -= p2.energy(pair);
            }
        }
        return u;
    }
}
