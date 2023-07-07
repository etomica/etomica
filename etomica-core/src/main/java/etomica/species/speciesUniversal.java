package etomica.species;

import etomica.atom.IAtom;
import etomica.molecule.IMolecule;
import etomica.space3d.Space3D;
import etomica.molecule.IMolecule;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpBuilder;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

public class speciesUniversal {
    //static String confName = "F:/ethane";
    public static ISpecies buildUniversal(boolean isDynamic, String confName) {
        // String confName = getConfName();
        ISpecies species = SpBuilder.buildSpecies(confName);
        IMolecule molecule = species.makeMolecule();
        SpeciesBuilder speciesBuilder = new SpeciesBuilder(Space3D.getInstance());
        for (IAtom atom : molecule.getChildList()) {
            System.out.println(atom.getType() + " " + atom.getPosition());
        }

        for (IAtom atom : molecule.getChildList()) {
            speciesBuilder.addAtom(atom.getType(), atom.getPosition());
            //speciesBuilder.withConformation(new ConfirmationMethane(Space3D.getInstance()));

        }
        speciesBuilder.setDynamic(isDynamic);
        return speciesBuilder.build();
    }
}
