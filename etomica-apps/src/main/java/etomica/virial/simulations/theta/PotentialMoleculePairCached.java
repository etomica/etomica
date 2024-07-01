/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial.simulations.theta;

import etomica.molecule.IMolecule;
import etomica.potential.PotentialMoleculePair;
import etomica.space.Space;
import etomica.species.SpeciesManager;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;

public class PotentialMoleculePairCached extends PotentialMoleculePair {

    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected BoxCluster box;

    public PotentialMoleculePairCached(Space space, SpeciesManager sm) {
        super(space, sm);
    }

    public void setBox(BoxCluster box) {
        this.box = box;
    }

    public double energy(IMolecule molecule1, IMolecule molecule2) {
        CoordinatePairSet cPairs = box.getCPairSet();
        long thisCPairID = cPairs.getID();
//            System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID-1000000) {
//                System.out.println("clusterSum "+cPairID+" returning recent "+value);
            return value;
        }
        if (thisCPairID == lastCPairID-10000000) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
//                System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
            return value;
        }

        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        value = super.energy(molecule1, molecule2);
        return value;
    }
}
