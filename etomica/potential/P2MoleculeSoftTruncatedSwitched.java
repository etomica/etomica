package etomica.potential;

import etomica.atom.AtomSet;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.box.Box;
import etomica.models.water.P2WaterSPCSoft;
import etomica.api.IVector;
import etomica.space.NearestImageTransformer;
import etomica.space.Tensor;
import etomica.space3d.Space3D;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  The energy
 * is switched from fully-on to 0 over a short range (i.e. from 0.95*rC to rC).
 */
public class P2MoleculeSoftTruncatedSwitched extends Potential2 implements IPotentialTorque {
    
    public P2MoleculeSoftTruncatedSwitched(IPotentialTorque potential, double truncationRadius) {
        super(potential.getSpace());
        this.potential = potential;
        setTruncationRadius(truncationRadius);
        zeroGradientAndTorque = new IVector[2][0];
        dr = space.makeVector();
    }
    
    /**
     * Returns the wrapped potential.
     */
    public IPotentialTorque getWrappedPotential() {
        return potential;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {return rCutoff;}
    
    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }

    public IVector[][] gradientAndTorque(AtomSet atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getAtom(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getAtom(0)).getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 < r2Cutoff) {
            IVector[][] gradientAndTorque = potential.gradientAndTorque(atoms);
            if (r2 > 0.95*0.95*r2Cutoff) {
                double r = Math.sqrt(r2);
                double fac = getF(r);
                // G = f G + u df/dr
                // tau = f tau
                gradientAndTorque[0][0].TE(fac);
                gradientAndTorque[0][1].TE(fac);
                gradientAndTorque[1][0].TE(fac);
                gradientAndTorque[1][1].TE(fac);
                // (df/dr)/r
                fac = getdFdr(r)/r;
                double u = potential.energy(atoms);
                gradientAndTorque[0][0].PEa1Tv1(-fac*u, dr);
                gradientAndTorque[0][1].PEa1Tv1(+fac*u, dr);
            }
            return gradientAndTorque;
        }
        if (zeroGradientAndTorque[0].length < atoms.getAtomCount()) {
            zeroGradientAndTorque[0] = new IVector[atoms.getAtomCount()];
            zeroGradientAndTorque[1] = new IVector[atoms.getAtomCount()];
            for (int i=0; i<atoms.getAtomCount(); i++) {
                zeroGradientAndTorque[0][i] = potential.getSpace().makeVector();
                zeroGradientAndTorque[1][i] = potential.getSpace().makeVector();
            }
        }
        return zeroGradientAndTorque;
    }
    
    protected double getF(double r) {
        switch (taperOrder) {
            case 1:
                return (rCutoff-r)/(rCutoff*0.05);
            case 2:
                return (r2Cutoff-2*rCutoff*r+r*r)/(r2Cutoff*0.05*0.05)+1e-7;
            case 3:
                double rt = 0.95*rCutoff;
                double a = (r-rt)/(rCutoff-rt);
                return (rCutoff-r)/(rCutoff-rt)*(1-a*a) + a*(1-a)*(1-a);
            default:
                throw new RuntimeException("oops");
        }
    }

    protected double getdFdr(double r) {
        switch (taperOrder) {
            case 1:
                return -1.0/(rCutoff*0.05);
            case 2:
                return -2 * (rCutoff - r) / (r2Cutoff*0.05*0.05);
            case 3:
                double rt = 0.95*rCutoff;
                double a = (r-rt)/(rCutoff-rt);
                double b = rCutoff-rt;
                double c = rCutoff-r;
                return (-(1.0-a*a)/b - 2*c*a/(b*b)) + (1-a)*(1-a)/b - 2*a/b + 2*a*a/b;
            default:
                throw new RuntimeException("oops");
        }
    }
    
    public static void main(String[] args) {
        P2MoleculeSoftTruncatedSwitched p = new P2MoleculeSoftTruncatedSwitched(new P2WaterSPCSoft(Space3D.getInstance()), 2);
        for (double x = 1.900001; x<1.999999; x+=0.001) {
            System.out.println(x+" "+p.getF(x)+" "+p.getdFdr(x));
        }
    }
    
    public IVector[] gradient(AtomSet atoms) {
        return gradientAndTorque(atoms)[0];
    }

    public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
        return gradientAndTorque(atoms)[0];
    }
    
    public double energy(AtomSet atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getAtom(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getAtom(0)).getPosition());
        nearestImageTransformer.nearestImage(dr);
        double r2 = dr.squared();
        if (dr.squared() > r2Cutoff) {
            return 0;
        }
        double u = potential.energy(atoms);
        if (r2 > 0.95*0.95*r2Cutoff) {
            u *= getF(Math.sqrt(r2));
        }
        return u;
    }
    
    public double virial(AtomSet atoms) {
        dr.Ev1Mv2(((MoleculeOrientedDynamic)atoms.getAtom(1)).getPosition(),((MoleculeOrientedDynamic)atoms.getAtom(0)).getPosition());
        nearestImageTransformer.nearestImage(dr);
        if (dr.squared() < r2Cutoff) {
            return potential.virial(atoms);
        }
        return 0;
    }
    
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public etomica.units.Dimension getTruncationRadiusDimension() {return etomica.units.Length.DIMENSION;}
    
    public void setBox(Box newBox) {
        potential.setBox(newBox);
        nearestImageTransformer = newBox.getBoundary();
    }
    
    private static final long serialVersionUID = 1L;
    protected double rCutoff, r2Cutoff;
    protected final IPotentialTorque potential;
    protected final IVector dr;
    protected NearestImageTransformer nearestImageTransformer;
    protected final IVector[][] zeroGradientAndTorque;
    protected int taperOrder = 2;
}
