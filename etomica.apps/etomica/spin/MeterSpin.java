package etomica.spin;

import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.DataSource;
import etomica.data.DataSourceScalar;
import etomica.space.Space;
import etomica.units.Undefined;


/**
 * Meter that provides the x-component of the vector average of
 * spin values (which is represented by the atom's position
 * vector).
 *
 * @author David Kofke
 *
 */
public class MeterSpin extends DataSourceScalar implements DataSource {

    /**
     * 
     */
    public MeterSpin(Space space) {
        super("Spin",Undefined.DIMENSION);
        sum = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.MeterScalar#getDataAsScalar(etomica.Box)
     */
    public double getDataAsScalar() {
        sum.E(0.0);
        int count = 0;
        iterator.setBox(box);
        iterator.reset();
        for (IAtomPositioned atom = (IAtomPositioned)iterator.nextAtom(); atom != null;
             atom = (IAtomPositioned)iterator.nextAtom()) {
            sum.PE(atom.getPosition());
            count++;
        }
        return sum.x(0)/count;
    }

    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private final IVector sum;
}
