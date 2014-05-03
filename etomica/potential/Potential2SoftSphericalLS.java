package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.space.BoundaryRectangularNonperiodic;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public class Potential2SoftSphericalLS extends Potential2 implements PotentialSoft{
   
    public Potential2SoftSphericalLS(ISpace space, double rCut, double a0, Potential2Soft p2Soft) {
         super(space);
        gradient = new IVectorMutable[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        this.rCut2 = rCut*rCut;
        this.a0 = a0;
        this.p2Soft = p2Soft;
    	Lxyz = space.makeVector();
		drtmp = space.makeVector();
		nShells = (int) Math.ceil(rCut/a0 - 0.49);
    }
        
    public double energy(IAtomList atoms) {
		double tmpU = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        for(int nx = -nShells; nx <= nShells; nx++) {
        	Lxyz.setX(0, nx*a0);
            for(int ny = -nShells; ny <= nShells; ny++) {
            	Lxyz.setX(1, ny*a0);
                for(int nz = -nShells; nz <= nShells; nz++) {
                	Lxyz.setX(2, nz*a0);
					drtmp.Ev1Pv2(dr, Lxyz);
					double dr2 = drtmp.squared();
					if(dr2 > rCut2 ) continue;
					if(dr2 != 0.0){
						tmpU += p2Soft.u(dr2);
					}
                }
            }
        }
        return tmpU;
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        double tmpVir = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
		for(int nx = -nShells; nx <= nShells; nx++) {
        	Lxyz.setX(0, nx*a0);
            for(int ny = -nShells; ny <= nShells; ny++) {
            	Lxyz.setX(1, ny*a0);
                for(int nz = -nShells; nz <= nShells; nz++) {
                	Lxyz.setX(2, nz*a0);
					drtmp.Ev1Pv2(dr, Lxyz);
					tmpVir += p2Soft.du(drtmp.squared());
                }
            }
        }
      return tmpVir;
    }
    
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public IVector[] gradient(IAtomList atoms) {
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        gradient[0].E(0);
        gradient[1].E(0);
        for(int nx = -nShells; nx <= nShells; nx++) {
        	Lxyz.setX(0, nx*a0);
            for(int ny = -nShells; ny <= nShells; ny++) {
            	Lxyz.setX(1, ny*a0);
                for(int nz = -nShells; nz <= nShells; nz++) {
                	Lxyz.setX(2, nz*a0);
					drtmp.Ev1Pv2(dr, Lxyz);
			        gradient[1].PEa1Tv1(p2Soft.du(drtmp.squared())/drtmp.squared(),drtmp);
                }
            }
        }
        gradient[0].PEa1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        gradient(atoms);
        pressureTensor.PEv1v2(gradient[0],dr);
        return gradient;
    }
    
    
    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
        p2Soft.setBox(box);
    }

    protected final IVectorMutable[] gradient;
    protected IBoundary boundary;
    protected final int nShells;
    protected final double a0;
    protected final Potential2Soft p2Soft;
    protected final IVectorMutable Lxyz;
    protected final IVectorMutable dr;
    protected final IVectorMutable drtmp;
    protected final double rCut2;
    

}//end of Potential2SoftSpherical
