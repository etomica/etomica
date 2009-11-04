package etomica.virial;

import etomica.models.water.SpeciesWater3P;
import etomica.api.ISpecies;
import etomica.space.ISpace;


/**
 * SpeciesFactory that makes SpeciesWater
 */
public class SpeciesFactoryWater3P implements SpeciesFactory, java.io.Serializable {
    public ISpecies makeSpecies(ISpace space) {
        return new SpeciesWater3P(space);
    }
}
