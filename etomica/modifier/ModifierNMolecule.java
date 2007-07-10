package etomica.modifier;

import etomica.box.Box;
import etomica.species.Species;
import etomica.units.Dimension;
import etomica.units.Quantity;

/**
 * Modifier class that enables change of the number of molecules of a particular species
 * in a particular box.
 */
public class ModifierNMolecule implements Modifier, java.io.Serializable {

    /**
     * @param speciesAgent Agent of the affected species in the affected box.
     * Cannot be changed after construction.
     */
    public ModifierNMolecule(Box box, Species species) {
        this.box = box;
        this.species = species;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        previousValue = mostRecentValue;
        mostRecentValue = (int)d;
        box.setNMolecules(species, (int) d);
    }

    public double getValue() {
        return box.getNMolecules(species);
    }

    public Dimension getDimension() {
        return Quantity.DIMENSION;
    }
    
    public String getLabel() {
        return species + " molecules";
    }
    
    public String toString() {
        return "Change number of "+species+
                " molecules from " + previousValue + " to " + mostRecentValue;
    }
    private static final long serialVersionUID = 1L;
    protected final Box box;
    protected final Species species;
    protected int mostRecentValue, previousValue;
}
