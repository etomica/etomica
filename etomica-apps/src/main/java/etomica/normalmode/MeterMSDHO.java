package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMSDHO extends DataSourceScalar {
    protected Box box;

    public MeterMSDHO(int nBeads, Box box) {
        super("PI MSD", Null.DIMENSION);
        this.box = box;
    }

    @Override
    public double getDataAsScalar() {
        double r2sum = 0;
        Vector dr = box.getSpace().makeVector();
        for (IMolecule molecule : box.getMoleculeList()) {
            Vector com = CenterOfMass.position(box, molecule);
            for (IAtom atom : molecule.getChildList()) {
                dr.Ev1Mv2(atom.getPosition(), com);
                box.getBoundary().nearestImage(dr);
                r2sum += dr.squared();
            }
        }
        return r2sum/box.getLeafList().size();
    }
}
