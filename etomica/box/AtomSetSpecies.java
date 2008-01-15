package etomica.box;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.species.Species;
import etomica.util.Debug;

/**
 * AtomSet that wraps a list of all molecules and makes it look like the
 * molecules of a single Species.
 * 
 * @author Andrew Schultz
 */
public class AtomSetSpecies implements AtomSet {

    public AtomSetSpecies(AtomArrayList moleculeList, Species species) {
        this.moleculeList = moleculeList;
        this.species = species;
    }
    
    public void setMoleculeStartIndices(int[] moleculeStartIndices) {
        int speciesIndex = species.getIndex();
        startIndex = moleculeStartIndices[speciesIndex];
        atomCount = moleculeStartIndices[speciesIndex+1] - startIndex;
    }

    public IAtom getAtom(int i) {
        if (Debug.ON && i < 0 || i > atomCount) {
            throw new ArrayIndexOutOfBoundsException("Invalid index");
        }
        return moleculeList.getAtom(i + startIndex);
    }

    public int getAtomCount() {
        return atomCount;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomArrayList moleculeList;
    protected final Species species;
    protected int startIndex, atomCount;
}
