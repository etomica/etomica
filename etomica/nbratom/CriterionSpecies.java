/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.nbratom;

import etomica.AtomPair;
import etomica.Species;
import etomica.nbr.NeighborCriterion;

/**
 * Filters atoms pairs to match a given pair of species.
 * 
 * @author andrew
 */
public class CriterionSpecies extends CriterionAdapter {

    public CriterionSpecies(NeighborCriterion criterion, 
            Species species0, Species species1) {
        super(criterion);
        this.species0 = species0;
        this.species1 = species1;
        isIntraSpecies = species0 == species1;
    }
    
    /**
     * Returns true if the species for the pair of atoms match 
     * the species given at construction (without regard to the
     * order of the pair), and if the wrapped criterion accept
     * also returns true.
     */
    public boolean accept(AtomPair pair) {
        Species atom0Species = pair.atom0.type.getSpecies();
        Species atom1Species = pair.atom1.type.getSpecies();
        if( (atom0Species == species0 && atom1Species == species1) 
               || (atom0Species == species1 && atom1Species == species0) ) {
            return subCriterion.accept(pair);
        }
        return false;
    }
    
    private final boolean isIntraSpecies;
    private final Species species0;
    private final Species species1;

}
