package etomica.box;

import etomica.api.IBox;
import etomica.api.IBoxMoleculeCountEvent;
import etomica.api.ISpecies;

public class BoxMoleculeCountEvent extends BoxEvent implements IBoxMoleculeCountEvent {

    public BoxMoleculeCountEvent(IBox box, ISpecies _species, int _count) {
        super(box);
        this.species = _species;
        this.count = _count;
    }
    
    public ISpecies getSpecies() {
        return species;
    }
    
    public int getCount() {
        return count;
    }
    
    protected int count = -1;
    protected ISpecies species = null;
}
