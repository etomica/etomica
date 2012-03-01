package etomica.atom;

import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IRandom;

/**
 * MoleculeSource that returns a completely random molecule.
 */
public class MoleculeSourceRandomMolecule implements MoleculeSource, java.io.Serializable {

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
    
    public void setBox(IBox p) {
        box = p;
    }
    
    /**
     * returns a random molecule from the box
     */
    public IMolecule getMolecule() {
        return box.getMoleculeList().getMolecule(random.nextInt(box.getMoleculeList().getMoleculeCount()));
    }
    
    private static final long serialVersionUID = 1L;
    protected IBox box = null;
    protected IRandom random;
}
