/*
 * Created on May 18, 2005
 */
package etomica.models.hexane;

import etomica.chem.elements.ElementSimple;
import etomica.simulation.Simulation;

/**
 * Species used to create hexane molecules per Dr. Monson's data. Hydrogen
 * molecules are ignored.
 * 
 * @author nancycribbin
 */

public class SpeciesHexane extends etomica.species.SpeciesSpheres {
    public SpeciesHexane(Simulation sim){
        super(sim, 6, new ElementSimple("M", sim.getDefaults().atomMass), 
              new ConformationHexane(sim.space));
    }
    
}
