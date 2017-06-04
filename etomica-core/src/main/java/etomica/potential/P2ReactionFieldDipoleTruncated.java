package etomica.potential;

import etomica.atom.IAtom;
import etomica.space.Boundary;
import etomica.atom.IMoleculePositionDefinition;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMolecular;
import etomica.space.Vector;
import etomica.atom.DipoleSource;
import etomica.space.Space;
import etomica.space.Tensor;

public class P2ReactionFieldDipoleTruncated extends PotentialMolecular implements PotentialMolecularSoft, IPotentialMolecularTorque {

    public P2ReactionFieldDipoleTruncated(Space space, IMoleculePositionDefinition positionDefinition) {
        super(2, space);
        this.positionDefinition = positionDefinition;
        iDipole = space.makeVector();
        cavityDipole = space.makeVector();
        dr = space.makeVector();
        gradientAndTorque = new Vector[2][2];
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
        boundary = box.getBoundary();
        if (cutoffRatio > 0){
        	Vector vectorBox = boundary.getBoxSize();
        	double minBoxSize = vectorBox.getX(0);
        	for (int i = 1;i<vectorBox.getD();i++){
        		if (vectorBox.getX(i) < minBoxSize){
        			minBoxSize = vectorBox.getX(i);
        		}
        	}
        	setRange(minBoxSize * cutoffRatio);
        }
    }

    public void setCutoffRatio(double newCutoffRatio) {
    	cutoffRatio = newCutoffRatio;
    }
    public double energy(IMoleculeList atoms) {
        dr.E(positionDefinition.position(atoms.getMolecule(1)));
        dr.ME(positionDefinition.position(atoms.getMolecule(0)));
        boundary.nearestImage(dr);
        if (dr.squared() > cutoff2) {
            return 0;
        }
        iDipole.E(dipoleSource.getDipole(atoms.getMolecule(0)));
        double idotj = iDipole.dot(dipoleSource.getDipole(atoms.getMolecule(1)));

        return -fac*idotj;
    }
    
    public Vector[][] gradientAndTorque(IMoleculeList atoms) {
        iDipole.E(dipoleSource.getDipole(atoms.getMolecule(0)));

        iDipole.XE(dipoleSource.getDipole(atoms.getMolecule(1)));
        iDipole.TE(fac);
        gradientAndTorque[0][0].E(0);
        gradientAndTorque[0][1].E(0);
        gradientAndTorque[1][0].E(iDipole);
        gradientAndTorque[1][1].Ea1Tv1(-1,iDipole);

        return gradientAndTorque;
    }

    public Vector[] gradient(IMoleculeList atoms) {
        return gradientAndTorque[0];
    }

    public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IMoleculeList atoms) {
        return 0;
    }

    /**
     * Returns a 0-body potential that should be added along with this
     * potential.
     */
    public IPotentialMolecular makeP0() {
        return new P0ReactionField(this.space, this);
    }

    private static final long serialVersionUID = 1L;
    protected final Vector iDipole, cavityDipole;
    protected final Vector dr;
    protected DipoleSource dipoleSource;
    protected Boundary boundary;
    protected double cutoff2, cutoff;
    protected double epsilon;
    protected final Vector[][] gradientAndTorque;
    protected double fac;
    protected double cutoffRatio;
    protected final IMoleculePositionDefinition positionDefinition;
    
    /**
     * A 0-body potential that should be added along with this potential.  The
     * 0-body potential includes the effective self-interaction of the
     * molecules (the molecule induces a dipole in the surrounding fluid, which
     * has an interaction energy with the molecule).  This part of the
     * potential does not result in a gradient or torque on the molecule and is
     * independent of position or orientation.
     */
    public static class P0ReactionField extends PotentialMolecular implements IPotential0Lrc, PotentialMolecularSoft {

        public P0ReactionField(Space space, P2ReactionFieldDipoleTruncated p) {
            super(0,space);
            this.potential = p;
            gradient = new Vector[0];
        }
        
        public double energy(IMoleculeList atoms) {
            double epsilon = potential.getDielectric();
            double cutoff = potential.getRange();
            DipoleSource dipoleSource = potential.getDipoleSource();
            double fac = 2*(epsilon-1)/(2*epsilon+1)/(cutoff*cutoff*cutoff);
            double u = 0;
            if (targetAtom != null) {
                Vector iDipole = dipoleSource.getDipole(targetAtom);
                u = -0.5 * fac * iDipole.squared();
            }
            else {
                IMoleculeList moleculeList = box.getMoleculeList();
                for (int i=0; i<moleculeList.getMoleculeCount(); i++) {
                    Vector iDipole = dipoleSource.getDipole(moleculeList.getMolecule(i));
                    u += -0.5 * fac * iDipole.squared();
                }
            }
            return u;
        }
        
        public void setBox(Box newBox) {
            box = newBox;
        }
        
        public void setTargetMolecule(IMolecule atom) {
            if (atom == null) {
                targetAtom = null;
                return;
            }
            targetAtom = atom;
        }
        
        public void setTargetAtom(IAtom targetAtom) {
            throw new RuntimeException("Can't provide correction for an individual atom");
        }

        public double getRange() {
            return 0;
        }

        public Vector[] gradient(IMoleculeList atoms) {
            return gradient;
        }
        
        public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
            return gradient(atoms);
        }
        
        public double virial(IMoleculeList atoms) {
            return 0;
        }

        private static final long serialVersionUID = 1L;
        protected final P2ReactionFieldDipoleTruncated potential;
        protected final Vector[] gradient;
        protected IMolecule targetAtom;
        protected Box box;
        

    }
}
