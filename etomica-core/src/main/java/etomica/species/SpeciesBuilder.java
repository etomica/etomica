

package etomica.species;


import etomica.atom.AtomType;
import etomica.config.ConformationGeneric;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Creates a species using parameters.
 *
 * @author navneeth
 */

public class SpeciesBuilder {

    public static Species SpeciesBuilder(Space space, AtomType[] atomTypes, int[] atomCount, Vector[] position){

        ConformationGeneric conformation = new ConformationGeneric(position);
        Species species;
        if(atomTypes.length == 1) {
            if(atomCount[0] > 1){
                species = new SpeciesSpheres(space, atomCount[0],atomTypes[0],conformation);
            }
            else {
                species = new SpeciesSpheresMono(space,atomTypes[0]);
            }
        }
        else {
            species = new SpeciesSpheresHetero(space, atomTypes);
            species.setConformation(conformation);
            ((SpeciesSpheresHetero)species).setChildCount(atomCount);
        }
        return species;
    }
}
