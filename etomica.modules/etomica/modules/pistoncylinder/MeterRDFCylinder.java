package etomica.modules.pistoncylinder;

import etomica.atom.iterator.ApiLeafAtoms;
import etomica.data.Data;
import etomica.data.meter.MeterRDF;
import etomica.modules.pistoncylinder.ApiFilteredCylinder.AtomFilterInCylinder;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * MeterRDF sublcass that properly calculates the RDF for the Piston/Cylinder
 * apparatus.  Including all pairs would yield a very non-bulk RDF since atoms
 * near the wall have open space next tho them.  To account for this, if
 * "padding" is the collision radius + RDF cutoff, then the RDF calculation
 * includes all pairs where both are a distance "padding" away from all walls
 * and half the pairs where only one is a distance "padding" away from all
 * walls.  The RDF is normalized to account for this.
 *
 * @author Andrew Schultz
 */
public class MeterRDFCylinder extends MeterRDF {

    public MeterRDFCylinder(Space space) {
        super(space);
    }
    
    public void setPotential(P1HardMovingBoundary newPistonPotential) {
        pistonPotential = newPistonPotential;
        reset();
    }
    
    public void reset() {
        super.reset();
        // make a new iterator with a new filter.  xMax might have changed
        AtomFilterInCylinder filter = new AtomFilterInCylinder(phase.getBoundary(), pistonPotential, xDataSource.getXMax());
        iterator = new ApiFilteredCylinder(new ApiLeafAtoms(), filter);
    }

    public Data getData() {
        super.getData();

        // renormalize the RDF to account for the excluded pairs
        double pistonRatio = 1;
        IVector dimensions = phase.getBoundary().getDimensions();
        double radius = pistonPotential.getCollisionRadius();
        for (int i=0; i<space.D(); i++) {
            if (i != 1) {
                pistonRatio *= (dimensions.x(i) - 2*(xMax + radius)) / (dimensions.x(i) - 2*radius);
            }
            else {
                pistonRatio *= ((dimensions.x(i)*0.5 - pistonPotential.getWallPosition()) - 2*(xMax + radius)) / 
                               ((dimensions.x(i)*0.5 - pistonPotential.getWallPosition()) - 2*radius);
            }
        }
        data.TE(1/pistonRatio);
        return data;
    }
    
    private static final long serialVersionUID = 1L;
    protected P1HardMovingBoundary pistonPotential;
}
