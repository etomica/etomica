package etomica.normalmode;

import etomica.EtomicaInfo;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.box.Box;
import etomica.potential.Potential2;
import etomica.potential.Potential2HardSpherical;
import etomica.potential.Potential2Spherical;
import etomica.potential.PotentialHard;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Hard potential that enforces ordering of the x-coordinates of the
 * pairs.  Returns infinite energy if the difference in atom indexes
 * and difference in x coordinates are of opposite sign; returns
 * zero otherwise.  Designed for use in 1D simulations.
 * 
 * @author David Kofke
 * @author Jhumpa Adhikari
 */
public class P2XOrder extends Potential2 implements Potential2Spherical, PotentialHard {
    
    private static final long serialVersionUID = 1L;
    protected final IVector dr;
    protected Box box;
    protected Potential2HardSpherical wrappedPotential;
    
    public P2XOrder(Space space, Potential2HardSpherical wrappedPotential) {
        super(space);
        dr = space.makeVector();
        this.wrappedPotential = wrappedPotential;
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Potential that enforces ordering in x coordinate (meant for 1D simulations)");
        return info;
    }

    /**
     * Interaction energy of the pair.
     * Zero if x coordinates are ordered differently from atom indexes.
     */
    public double energy(AtomSet pair) {
        IAtom atom0 = pair.getAtom(0);
        IAtom atom1 = pair.getAtom(1);
        dr.Ev1Mv2(((IAtomPositioned)atom1).getPosition(), ((IAtomPositioned)atom0).getPosition());
        int dI = atom1.getIndex() - atom0.getIndex();
        // assume 1 species
        if (Math.abs(dI) == box.getMoleculeList().getAtomCount()-1) {
            dr.PEa1Tv1(dI > 0 ? -1 : 1, box.getBoundary().getDimensions());
            return (dr.x(0) * dI > 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else if (dI == 1 || dI == -1) {
            return (dr.x(0) * dI < 0.0) ? Double.POSITIVE_INFINITY : wrappedPotential.u(dr.squared());
        }
        else {
            return 0;
        }
    }
    
    /**
     * Returns infinity.
     */
    public double getRange() {
        return wrappedPotential.getRange();
    }
    
    public Potential2Spherical getWrappedPotential() {
        return wrappedPotential;
    }

    public void setBox(Box newBox) {
        box = newBox;
        wrappedPotential.setBox(newBox);
    }

    public void bump(AtomSet atom, double falseTime) {
        wrappedPotential.bump(atom, falseTime);
    }

    public double collisionTime(AtomSet atom, double falseTime) {
        return wrappedPotential.collisionTime(atom, falseTime);
    }

    public double energyChange() {
        return wrappedPotential.energyChange();
    }

    public double lastCollisionVirial() {
        return wrappedPotential.lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return wrappedPotential.lastCollisionVirialTensor();
    }
    public double u(double r2) {
        return wrappedPotential.u(r2);
    }
}