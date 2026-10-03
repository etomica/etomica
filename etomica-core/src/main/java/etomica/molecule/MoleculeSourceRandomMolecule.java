/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.molecule;

import etomica.box.Box;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * MoleculeSource that returns a completely random molecule.
 */
public class MoleculeSourceRandomMolecule implements MoleculeSource, java.io.Serializable {

    protected Box box = null;
    protected IRandom random;
    protected ISpecies species;

    public MoleculeSourceRandomMolecule() {

    }

    public MoleculeSourceRandomMolecule(Box box, IRandom random) {
        this.box = box;
        this.random = random;
    }

    public void setSpecies(ISpecies s) {
        this.species = s;
    }

    /**
     * Sets the random number generator used to pick molecules
     */
    public void setRandomNumberGenerator(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick molecules
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setBox(Box p) {
        box = p;
    }
    
    /**
     * returns a random molecule from the box
     */
    public IMolecule getMolecule() {
        if (species==null) return box.getMoleculeList().get(random.nextInt(box.getMoleculeList().size()));
        int n = box.getNMolecules(species);
        if (n == 0) return null;
        return box.getMoleculeList(species).get(random.nextInt(n));
    }
    
}
