/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.nbr.molecule;

import etomica.molecule.IMoleculeList;
import etomica.species.ISpecies;

/**
 * Filters molecule pairs to match a given pair of Species.
 * 
 * @author Tai Boon Tan
 */
public class CriterionSpeciesPair extends CriterionAdapterMolecular {

    public CriterionSpeciesPair(NeighborCriterionMolecular criterion, 
            ISpecies species0, ISpecies species1) {
        super(criterion);
        this.species0 = species0;
        this.species1 = species1;
    }
    
    /**
     * Returns true if the Species for the pair of molecules match the Species 
     * given at construction (without regard to the order of the pair), and if 
     * the wrapped criterion also accepts the pair.
     */
    public boolean accept(IMoleculeList pair) {
        ISpecies molecule0Species = pair.getMolecule(0).getType();
        ISpecies molecule1Species = pair.getMolecule(1).getType();
        if ( (molecule0Species == species0 && molecule1Species == species1) ||
             (molecule0Species == species1 && molecule1Species == species0) ) {
            return subCriterion.accept(pair);
        }
        return false;
    }
    
    /**
     * Returns the Species accepted by this NeighborCriterionMolecular
     */
    public ISpecies[] getSpecies() {
        return new ISpecies[]{species0, species1};
    }
    
    private static final long serialVersionUID = 1L;
    private final ISpecies species0;
    private final ISpecies species1;
}
