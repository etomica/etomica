package etomica.modifier;

import etomica.atom.ISpeciesAgent;
import etomica.atom.SpeciesAgent;
import etomica.units.Dimension;
import etomica.units.Quantity;

/**
 * Modifier class that enables change of the number of molecules of a particular species
 * in a particular phase.
 */
public class ModifierNMolecule implements Modifier, java.io.Serializable {

    /**
     * @param speciesAgent Agent of the affected species in the affected phase.
     * Cannot be changed after construction.
     */
    public ModifierNMolecule(ISpeciesAgent speciesAgent) {
        this.speciesAgent = speciesAgent;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        previousValue = mostRecentValue;
        mostRecentValue = (int)d;
        speciesAgent.setNMolecules((int) d);
//        if (this.selector.display != null)
//            this.selector.display.repaint();
//        this.selector.integrator.reset();
    }

    public double getValue() {
        return (speciesAgent != null) ? (double)speciesAgent.getNMolecules() : 0;
    }

    public Dimension getDimension() {
        return Quantity.DIMENSION;
    }
    
    public String getLabel() {
        return speciesAgent.getType().getSpecies().getName() + " molecules";
    }
    
    public String toString() {
        return "Change number of "+speciesAgent.getType().getSpecies().getName()+
                " molecules from " + previousValue + " to " + mostRecentValue;
    }
    private static final long serialVersionUID = 1L;
    private final ISpeciesAgent speciesAgent;
    private int mostRecentValue, previousValue;
}
