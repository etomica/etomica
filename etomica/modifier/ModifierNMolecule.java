package etomica.modifier;

import etomica.Modifier;
import etomica.SpeciesAgent;
import etomica.units.Dimension;

/**
 * Modifier class that enables change of the number of molecules of a particular species
 * in a particular phase.
 */
/*
 * History Created on Jan 31, 2005 by kofke
 */
public class ModifierNMolecule implements Modifier {

    /**
     * @param speciesAgent Agent of the affected species in the affected phase.
     * Cannot be changed after construction.
     */
    public ModifierNMolecule(SpeciesAgent speciesAgent) {
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
        return (speciesAgent != null) ? (double)speciesAgent.moleculeCount() : 0;
    }

    public Dimension getDimension() {
        return Dimension.QUANTITY;
    }
    
    public String getLabel() {
        return speciesAgent.node.parentSpecies().getName() + " molecules";
    }
    
    public String toString() {
        return "Change number of "+speciesAgent.node.parentSpecies().getName()+
                " molecules from " + previousValue + " to " + mostRecentValue;
    }
    private final SpeciesAgent speciesAgent;
    private int mostRecentValue, previousValue;
}