/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.util.random.IRandom;
import etomica.api.ISpecies;

/**
 * AtomSource that returns a completely random leaf atom.
 */
public class AtomSourceRandomSpecies implements AtomSource {

    public AtomSourceRandomSpecies() {}
    
    public AtomSourceRandomSpecies(IRandom random, ISpecies s) {
        setRandomNumberGenerator(random);
        setSpecies(s);
    }
    
    /**
     * Sets the random number generator used to pick atoms
     */
    public void setRandomNumberGenerator(IRandom newRandom) {
        random = newRandom;
    }
    
    /**
     * Returns the random number generator used to pick atoms
     */
    public IRandom getRandomNumberGenerator() {
        return random;
    }
    
    public void setSpecies(ISpecies s) {
        species = s;
        if (box != null) list = box.getMoleculeList(species);
    }
    
    public void setBox(Box p) {
        box = p;
        if (species != null) list = box.getMoleculeList(species);
    }
    
    /**
     * returns a random atom from the box's leaf atom list
     */
    public IAtom getAtom() {
        int n = list.getMoleculeCount();
        if (n == 0) return null;
        int r = random.nextInt(n);
        IMolecule m = list.getMolecule(r);
        IAtomList atoms = m.getChildList();
//        IAtomList atoms = list.getMolecule(random.nextInt(n)).getChildList();
        n = atoms.getAtomCount();
        if (n == 0) return null;
        return atoms.getAtom(random.nextInt(n));
    }
    
    protected IMoleculeList list = null;
    protected IRandom random;
    protected Box box;
    protected ISpecies species;
}
