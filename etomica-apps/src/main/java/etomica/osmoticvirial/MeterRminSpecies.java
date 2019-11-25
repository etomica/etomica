package etomica.osmoticvirial;

import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.molecule.iterator.MpiIntraspeciesAA;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Length;

/**
 * Calculates minimum pair distance from all the possible pairs of atoms of a particular species
 * as required by Ashton and Wilding's method
 */
public class MeterRminSpecies extends DataSourceScalar {
    protected AtomType type1, type2;
    protected MpiIntraspeciesAA iterator;
    protected Vector dr;
    protected Boundary boundary;

    public MeterRminSpecies(Space space, Box box, ISpecies species){
        super("RminSpecies", Length.DIMENSION);
        dr = space.makeVector();
        iterator = new MpiIntraspeciesAA(species);
        iterator.setBox(box);
        boundary = box.getBoundary();
    }

    @Override
    public double getDataAsScalar() {
        double rminSq = Double.POSITIVE_INFINITY;
        iterator.reset();
        for (IMoleculeList pair = iterator.next(); pair != null;
             pair = iterator.next()) {
            if (type1 != null && (pair.getMolecule(0).getType() != type1 || pair.getMolecule(1).getType() != type2)) continue;
            dr.Ev1Mv2(pair.getMolecule(1).getChildList().getAtom(0).getPosition(),pair.getMolecule(0).getChildList().getAtom(0).getPosition());
            boundary.nearestImage(dr);
            double r2 = dr.squared();
            if (rminSq > r2) rminSq = r2;
        }
        return Math.sqrt(rminSq);
    }
}
