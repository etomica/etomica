package etomica.spin;

import etomica.Atom;
import etomica.Phase;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorListTabbed;
import etomica.data.DataSourceScalar;
import etomica.space.Vector;
import etomica.units.Dimension;


/**
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 25, 2005 by kofke
 */
public class MeterSpin extends DataSourceScalar {

    /**
     * 
     */
    public MeterSpin(Space space) {
        super();
        sum = space.makeVector();
    }

    /* (non-Javadoc)
     * @see etomica.data.meter.MeterScalar#getDataAsScalar(etomica.Phase)
     */
    public double getDataAsScalar(Phase phase) {
        sum.E(0.0);
        int count = 0;
        iterator.setList(phase.speciesMaster.atomList);
        iterator.reset();
        while(iterator.hasNext()) {
            Atom atom = iterator.nextAtom();
            sum.PE(atom.coord.position());
            count++;
        }
        return sum.x(0)/count;
    }

    /* (non-Javadoc)
     * @see etomica.DataSource#getDimension()
     */
    public Dimension getDimension() {
        return Dimension.UNDEFINED;
    }

    private final AtomIteratorListTabbed iterator = new AtomIteratorListTabbed();
    private final Vector sum;
}
