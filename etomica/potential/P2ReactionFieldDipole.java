package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.AtomType;
import etomica.atom.DipoleSource;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.api.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.IVector3D;

public class P2ReactionFieldDipole extends Potential2 implements PotentialSoft, IPotentialTorque {

    public P2ReactionFieldDipole(Space space) {
        super(space);
        iDipole = (IVector3D)space.makeVector();
        cavityDipole = (IVector3D)space.makeVector();
        dr = space.makeVector();
        gradientAndTorque = new IVector[2][2];
        gradientAndTorque[0][0] = space.makeVector();
        gradientAndTorque[0][1] = space.makeVector();
        gradientAndTorque[1][0] = space.makeVector();
        gradientAndTorque[1][1] = space.makeVector();
    }

    /**
     * Returns the dipole source used by this object.
     */
    public DipoleSource getDipoleSource() {
        return dipoleSource;
    }

    /**
     * Sets the dipole source used by this object should use.
     */
    public void setDipoleSource(DipoleSource newDipoleSource) {
        dipoleSource = newDipoleSource;
    }
    
    public double getRange() {
        return cutoff;
    }
    
    public void setRange(double newRange) {
        cutoff = newRange;
        cutoff2 = newRange * newRange;
        fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
    }
    
    /**
     * Returns the dielectric constant of the fluid surrounding the cavity.
     */
    public double getDielectric() {
        return epsilon;
    }
    
    /**
     * Sets the dielectric constant of the fluid surrounding the cavity.
     */
    public void setDielectric(double newDielectric) {
        epsilon = newDielectric;
        if (cutoff > 0) {
            fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff2*cutoff);
        }
    }

    public void setBox(Box box) {
        nearestImageTransformer = box.getBoundary();
    }

    public double energy(AtomSet atoms) {
        iDipole.E(dipoleSource.getDipole((IMolecule)atoms.getAtom(0)));
        double idotj = iDipole.dot(dipoleSource.getDipole((IMolecule)atoms.getAtom(1)));

        return -fac*idotj;
    }
    
    public IVector[][] gradientAndTorque(AtomSet atoms) {
        iDipole.E(dipoleSource.getDipole((IMolecule)atoms.getAtom(0)));

        iDipole.XE((IVector3D)dipoleSource.getDipole((IMolecule)atoms.getAtom(1)));
        iDipole.TE(fac);
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(iDipole);
        gradientAndTorque[1][1].Ea1Tv1(-1,iDipole);

        return gradientAndTorque;
    }

    public IVector[] gradient(AtomSet atoms) {
        return gradientAndTorque[0];
    }

    public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(AtomSet atoms) {
        return 0;
    }

    /**
     * Returns a 0-body potential that should be added along with this
     * potential.
     */
    public IPotential makeP0() {
        return new P0ReactionField(this);
    }

    private static final long serialVersionUID = 1L;
    protected final IVector3D iDipole, cavityDipole;
    protected final IVector dr;
    protected DipoleSource dipoleSource;
    protected NearestImageTransformer nearestImageTransformer;
    protected double cutoff2, cutoff;
    protected double epsilon;
    protected final IVector[][] gradientAndTorque;
    protected double fac;
    
    /**
     * A 0-body potential that should be added along with this potential.  The
     * 0-body potential includes the effective self-interaction of the
     * molecules (the molecule induces a dipole in the surrounding fluid, which
     * has an interaction energy with the molecule).  This part of the
     * potential does not result in a gradient or torque on the molecule and is
     * independent of position or orientation.
     */
    public static class P0ReactionField extends Potential0Lrc {

        public P0ReactionField(P2ReactionFieldDipole p) {
            super(p.getSpace(), new AtomType[2], p);
            this.potential = p;
            gradient = new IVector[0];
        }
        
        public double energy(AtomSet atoms) {
            double epsilon = potential.getDielectric();
            double cutoff = potential.getRange();
            DipoleSource dipoleSource = potential.getDipoleSource();
            double fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff*cutoff*cutoff);
            double u = 0;
            if (targetAtom != null) {
                IVector iDipole = dipoleSource.getDipole(targetAtom);
                u = -0.5 * fac * iDipole.squared();
            }
            else {
                AtomSet moleculeList = box.getMoleculeList();
                for (int i=0; i<moleculeList.getAtomCount(); i++) {
                    IVector iDipole = dipoleSource.getDipole((IMolecule)moleculeList.getAtom(i));
                    u += -0.5 * fac * iDipole.squared();
                }
            }
            return u;
        }
        
        public void setBox(Box newBox) {
            box = newBox;
        }
        
        public void setTargetAtoms(AtomSet atoms) {
            if (atoms == null || atoms.getAtom(0) == null || !(atoms.getAtom(0) instanceof IMolecule)) {
                targetAtom = null;
                return;
            }
            targetAtom = (IMolecule)atoms.getAtom(0);
        }
        
        public IVector[] gradient(AtomSet atoms) {
            return gradient;
        }
        
        public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
        
        public double virial(AtomSet atoms) {
            return 0;
        }

        private static final long serialVersionUID = 1L;
        protected final P2ReactionFieldDipole potential;
        protected final IVector[] gradient;
        protected IMolecule targetAtom;
    }
}
