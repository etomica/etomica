package etomica.box;

import etomica.api.IBox;
import etomica.api.ISpecies;


/**
 * Event that conveys that the maximum global index in a Box has changed.
 */
public class BoxNumMoleculesEvent extends BoxEvent {

    public BoxNumMoleculesEvent(IBox box, ISpecies species, int newNumMolecules) {
        super(box);
        this.species = species;
        this.newNumMolecules = newNumMolecules;
    }

    public ISpecies getSpecies() {
        return species;
    }
    
    public int getNewNumMolecules() {
        return newNumMolecules;
    }
    
    protected final ISpecies species;
    protected final int newNumMolecules;
    private static final long serialVersionUID = 1L;
}
