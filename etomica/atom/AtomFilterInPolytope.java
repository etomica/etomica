package etomica.atom;

import etomica.Atom;
import etomica.math.geometry.Polytope;


/**
 * Filter that accepts atom if it is inside a specified polytope.
 * Position of atom is determined by an AtomPositionDefinition, which
 * if unspecified default's for each atom to that given by the atom's
 * type. 
 *
 * @author David Kofke
 */

/*
 * History
 * Created on Jun 14, 2005 by kofke
 */
public class AtomFilterInPolytope implements AtomFilter {

    /**
     * Create filter in which position definition is given by atom's type.
     */
    public AtomFilterInPolytope(Polytope polytope) {
        super();
        this.polytope = polytope;
    }

    /**
     * Returns true if the atom's position is inside the polytope.
     */
    public boolean accept(Atom atom) {
        if(positionDefinition == null) {
            return polytope.contains(atom.type.getPositionDefinition().position(atom));
        }
        return polytope.contains(positionDefinition.position(atom));
    }

    /**
     * @return Returns the polytope.
     */
    public Polytope getPolytope() {
        return polytope;
    }
    /**
     * @param polytope The polytope to set.
     */
    public void setPolytope(Polytope polytope) {
        this.polytope = polytope;
    }
    /**
     * @return Returns the positionDefinition.
     */
    public AtomPositionDefinition getPositionDefinition() {
        return positionDefinition;
    }
    /**
     * Sets position definition.  If given null, positionDefinition from
     * each atom's type will be used.
     * @param positionDefinition The positionDefinition to set.
     */
    public void setPositionDefinition(AtomPositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }
    
    private Polytope polytope;
    private AtomPositionDefinition positionDefinition;
}
