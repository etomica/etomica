package etomica.virial;

import etomica.api.IAtomSet;
import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IVector;

import etomica.atom.AtomPair;
import etomica.potential.Potential2;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.Potential2Spherical;
import etomica.space.Space;

/**
 * Transformed potential.  For distances less than some cutoff, the potential
 * is untransformed.  For very long distances, the potential energy returned
 * is the subpotential's energy divided by r2.  The transformation is applied
 * slowly, beginning at the cutoff.
 * @author Andrew Schultz
 */
public class Potential2TransformedR2 extends Potential2 {

    /**
     * @param space
     */
    public Potential2TransformedR2(Space space, Potential2 p) {
        super(space);
        potential = p;
        dr = space.makeVector();
        setTransformDistance(1.0);
    }
    
    public void setTransformDistance(double newDistance) {
        r2Transform = newDistance*newDistance;
    }

    public double getTransformDistance() {
        return Math.sqrt(r2Transform);
    }

    public double energy(IAtomSet atoms) {
        double r2;
        if (((AtomPair)atoms).atom0 instanceof IAtomPositioned) {
            dr.Ev1Mv2(((IAtomPositioned)((AtomPair)atoms).atom1).getPosition(),((IAtomPositioned)((AtomPair)atoms).atom0).getPosition());
            r2 = dr.squared();
        }
        else {
            IAtom atom0 = ((AtomPair)atoms).atom0;
            IAtom atom1 = ((AtomPair)atoms).atom1;
            dr.E(atom0.getType().getPositionDefinition().position(atom0));
            dr.ME(atom1.getType().getPositionDefinition().position(atom1));
            r2 = dr.squared();
        }
        double e;
        if (potential instanceof Potential2SoftSpherical) {
            e = ((Potential2Spherical)potential).u(r2);
        }
        else {
            e = potential.energy(atoms);
        }
        if (r2 < r2Transform) {
            return e;
        }
        return e / (1 + (r2 - r2Transform));
    }

    public void setBox(IBox b) {
    }
    
    public double getRange() {return potential.getRange();}
    
    private static final long serialVersionUID = 1L;
    protected final Potential2 potential;
    protected final IVector dr;
    protected double r2Transform;
}
