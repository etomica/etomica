package etomica.osmoticvirial;

import com.fasterxml.jackson.annotation.JsonAnyGetter;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiLeafAtoms;
import etomica.atom.iterator.AtomsetIteratorBoxDependent;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.DataSourceUniform;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.action.IAction;
import etomica.units.Length;


/**
 * Calculates minimum pair distance of all the possible pairs
 *
 * Created by aksharag on 07-17-2017.
 */

public class MeterRmin extends DataSourceScalar {

    protected AtomType type1, type2;
    protected AtomsetIteratorBoxDependent iterator;
    protected Vector dr;
    protected Boundary boundary;


    public MeterRmin(Space space, Box box) {
        super("Rmin", Length.DIMENSION);
        dr = space.makeVector();
        iterator = new ApiLeafAtoms();
        iterator.setBox(box);
        boundary = box.getBoundary();
    }

    @Override
    public double getDataAsScalar() {
        double rminSq = Double.POSITIVE_INFINITY;
        iterator.reset();
        for (IAtomList pair = iterator.next(); pair != null;
             pair = iterator.next()) {
            if (type1 != null && (pair.getAtom(0).getType() != type1 || pair.getAtom(1).getType() != type2)) continue;
            dr.Ev1Mv2(pair.getAtom(1).getPosition(),pair.getAtom(0).getPosition());
            boundary.nearestImage(dr);
            double r2 = dr.squared();
            if (rminSq > r2) rminSq = r2;

        }
        return Math.sqrt(rminSq);
    }

}

