package etomica.modifier;

import etomica.Modifier;
import etomica.SpeciesAgent;
import etomica.units.Dimension;

/*
 * History Created on Jan 31, 2005 by kofke
 */
public class ModifierNMolecule extends Modifier {

    public ModifierNMolecule(SpeciesAgent speciesAgent) {
        super(Dimension.QUANTITY);
        this.speciesAgent = speciesAgent;
    }

    public void setValue(double d) {
        if (d < 0) d = 0;
        speciesAgent.setNMolecules((int) d);
//        if (this.selector.display != null)
//            this.selector.display.repaint();
//        this.selector.integrator.reset();
    }

    public double getValue() {
        return (speciesAgent != null) ? (double)speciesAgent.moleculeCount() : 0;
    }
    
    private final SpeciesAgent speciesAgent;
}