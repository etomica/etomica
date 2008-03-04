package etomica.species;

import etomica.api.ISpecies;


/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface SpeciesResolver {

    public ISpecies whichOneDoYouLike(ISpecies[] candidates);
    
}
