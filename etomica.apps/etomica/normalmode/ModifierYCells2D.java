package etomica.normalmode;

import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.modifier.Modifier;
import etomica.units.Dimension;
import etomica.units.Null;

/**
 * Modifier class that enables change of the number of cells in 2D
 * 
 */
public class ModifierYCells2D implements Modifier, java.io.Serializable {

    /**
     * @param speciesAgent Agent of the affected species in the affected box.
     * Cannot be changed after construction.
     */
    public ModifierYCells2D(IBox box, ISpecies species, int x) {
        this.box = box;
        this.species = species;
        this.xCell = x;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        previousValue = mostRecentValue;
        mostRecentValue = (int)d;
        
        box.setNMolecules(species, 2*(int)d*xCell);
    }

    public double getValue() {
        return box.getNMolecules(species)/(2*xCell);
    }

    public Dimension getDimension() {
        return Null.DIMENSION;
    }
    
    public String getLabel() {
        return "y-Cell Number";
    }
    
 
    private static final long serialVersionUID = 1L;
    protected final IBox box;
    protected final ISpecies species;
    protected int mostRecentValue, previousValue, xCell;
}
