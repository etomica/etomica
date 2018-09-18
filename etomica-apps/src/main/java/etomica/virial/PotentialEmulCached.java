/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialEmul;
import etomica.space.Space;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Potential class whose sole purpose is to remember previously calculated
 * energies and return them if no change has been made to the configuration.
 * No attempt is made to remember energies from any previous configuration.
 * 
 * @author Andrew Schultz
 */
public class PotentialEmulCached extends PotentialEmul {

    protected int lastCPairID = -1;
    protected final int maxNumAtoms;
    protected final double[] lastEnergies;
    protected final int[] ids;

    public PotentialEmulCached(Space space, String templateName, int maxNumAtoms) {
        this(space, templateName, 2.5, maxNumAtoms);
    }

    public PotentialEmulCached(Space space, String templateName, double rCore, int maxNumAtoms) {
        this(space, templateName, countMolecules(templateName), rCore, maxNumAtoms);
    }

    public PotentialEmulCached(Space space, String templateName, int nBody,
                               double rCore, int maxNumAtoms) {
        super(space, templateName, nBody, rCore);
        this.maxNumAtoms = maxNumAtoms;
        ids = new int[nBody];
        for (int j=0; j<nBody; j++) {
            ids[j] = maxNumAtoms-nBody+j;
        }
        lastEnergies = new double[VirialDiagrams.getGroupID(ids, maxNumAtoms)+1];
    }
    
    public double energy(IMoleculeList molecules) {
        for (int i=0; i<nBody; i++) {
            ids[i] = molecules.get(i).getIndex();
        }
        int gid = VirialDiagrams.getGroupID(ids, maxNumAtoms);
        if (Double.isNaN(lastEnergies[gid])) {
            lastEnergies[gid] = super.energy(molecules);
        }
        return lastEnergies[gid];
    }

    public void setBox(Box newBox) {
        boolean discard = true;
        if (newBox instanceof BoxCluster) {
            discard = lastCPairID != ((BoxCluster)newBox).getCPairID();
        }
        if (discard) {
            for (int j=0; j<lastEnergies.length; j++) {
                lastEnergies[j] = Double.NaN;
            }
        }
        super.setBox(newBox);
    }
}
