package etomica.atom;

import etomica.phase.Phase;
import etomica.simulation.Simulation;

public class AtomSourceRandomMoleculeSeq extends AtomSourceRandomMolecule {

    public void setPhase(Phase p) {
        speciesMaster = p.getSpeciesMaster();
        reset();
    }
    
    public Atom getAtom() {
        if (prevLinker == null || !prevLinker.atom.node.inSimulation()) {
            // no suitable previous atom to step forward from 
            reset();
            // no atoms in the phase
            if (prevLinker == null) return null;
        }
        int lookAhead = Simulation.random.nextInt(maxLookAhead+1);

        int i = 0;
        while (i < lookAhead) {
            for ( ; i<lookAhead && prevLinker.atom != null; i++) {
                prevLinker = prevLinker.next;
            }
            if (prevLinker.atom == null) {
                // wrapped back to first molecule of the species
                // ==> advance to next species
                SpeciesAgent s = prevLinker.next.atom.node.parentSpeciesAgent().nextSpecies();
                // last species ==> go to first species
                if (s == null)
                    s = speciesMaster.firstSpecies();
                prevLinker = ((AtomTreeNodeGroup)s.node).childList.firstEntry();
                // no molecules in that species, find one with molecules
                while (prevLinker == null) {
                    s = s.nextSpecies();
                    prevLinker = ((AtomTreeNodeGroup)s.node).childList.firstEntry();
                }
            }
        }
        return prevLinker.atom;
    }

    /**
     * Reset the atom used to step from to a random molecule
     */
    public void reset() {
        int size = speciesMaster.moleculeCount();
        if (size == 0) {
            prevLinker = null;
            return;
        }
        int jump = Simulation.random.nextInt(size);
        int sum = 0;
        SpeciesAgent s;
        for(s=speciesMaster.firstSpecies(); s!=null; s=s.nextSpecies()) {
            sum += s.node.childAtomCount();
            if(sum > jump) break;
        }
        prevLinker = ((AtomTreeNodeGroup)s.node).childList.entry(jump-(sum-((AtomTreeNodeGroup)s.node).childAtomCount()));
    }

    /**
     * Returns the maximum number of molecules the source will advance for
     * each call to getAtom
     */
    public int getMaxLookAhead() {
        return maxLookAhead;
    }
    /**
     * Sets the maximum number of molecules the source will advance for
     * each call to getAtom
     */
    public void setMaxLookAhead(int maxLookAhead) {
        this.maxLookAhead = maxLookAhead;
    }
    
    protected int maxLookAhead = 10;
    protected SpeciesMaster speciesMaster;
    protected AtomLinker prevLinker;
}