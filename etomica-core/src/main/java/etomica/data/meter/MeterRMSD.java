package etomica.data.meter;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Length;

/**
 * Measured root mean squared displacement of all atoms in the box.  This meter
 * will properly track atoms across periodic boundaries so long as atoms do not
 * move more than L/2 between calls to the meter.
 * <p>
 * Created by andrew on 6/15/17.
 */
public class MeterRMSD extends DataSourceScalar {
    
    protected final Box box;
    protected final Vector[] originalPosition;
    protected final Vector[] lastPosition;
    protected final Vector dr, drTmp;
    
    public MeterRMSD(Box box, Space space) {
        super("RMSD", Length.DIMENSION);
        this.box = box;
        IAtomList atoms = box.getLeafList();
        originalPosition = new Vector[atoms.getAtomCount()];
        lastPosition = new Vector[atoms.getAtomCount()];
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            originalPosition[i] = space.makeVector();
            lastPosition[i] = space.makeVector();
            originalPosition[i].E(atoms.getAtom(i).getPosition());
            lastPosition[i].E(originalPosition[i]);
        }
        dr = space.makeVector();
        drTmp = space.makeVector();
    }
    
    @Override
    public double getDataAsScalar() {
        double sum = 0;
        IAtomList atoms = box.getLeafList();
        Boundary boundary = box.getBoundary();
        for (int i = 0; i < atoms.getAtomCount(); i++) {
            Vector p = atoms.getAtom(i).getPosition();
            dr.Ev1Mv2(p, lastPosition[i]);
            drTmp.E(dr);
            boundary.nearestImage(drTmp);
            if (!dr.equals(drTmp)) {
                originalPosition[i].PE(dr);
                originalPosition[i].ME(drTmp);
            }
            sum += p.Mv1Squared(originalPosition[i]);
            lastPosition[i].E(p);
        }
        return Math.sqrt(sum / atoms.getAtomCount());
    }
}
