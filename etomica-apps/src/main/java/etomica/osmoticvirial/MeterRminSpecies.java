package etomica.osmoticvirial;

import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;

/**
 * Calculates minimum pair distance from all the possible pairs of atoms of a particular species
 * as required by Ashton and Wilding's method
 */
public class MeterRminSpecies extends DataSourceScalar {
    protected final Vector dr;
    protected final Box box;
    protected final ISpecies species;

    public MeterRminSpecies(Space space, Box box, ISpecies species){
        super("RminSpecies", Length.DIMENSION);
        dr = space.makeVector();
        this.box = box;
        this.species = species;
    }

    @Override
    public double getDataAsScalar() {
        double rminSq = Double.POSITIVE_INFINITY;
        IMoleculeList molecules = box.getMoleculeList(species);
        for (int i=0; i<molecules.size(); i++) {
            Vector pi = molecules.get(i).getChildList().get(0).getPosition();
            for (int j=i+1; j<molecules.size(); j++) {
                dr.Ev1Mv2(molecules.get(j).getChildList().get(0).getPosition(), pi);
                box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                if (rminSq > r2) rminSq = r2;
            }
        }
        return Math.sqrt(rminSq);
    }
}
