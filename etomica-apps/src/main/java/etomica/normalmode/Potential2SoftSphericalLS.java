package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialSoft;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Methods for a soft (non-impulsive), spherically-symmetric pair potential.
 * Subclasses must provide concrete definitions for the energy (method
 * u(double)) and its derivatives.
 * 
 * @author David Kofke
 */
 
public class Potential2SoftSphericalLS extends Potential2 implements PotentialSoft{
	public Potential2SoftSphericalLS(Space space, double rCut, double[] a0, Potential2Soft p2Soft){
		this(space, rCut, a0,p2Soft,null);
		}
    public Potential2SoftSphericalLS(Space space, double rCut, double[] a0, Potential2Soft p2Soft,MoleculeAgentManager latticeCoordinates) {
         super(space);
        gradient = new Vector[2];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        dr = space.makeVector();
        ldr = space.makeVector();
        shift = space.makeVector();
        r0 = space.makeVector();
        this.rCut2 = rCut*rCut;
        this.a0 = a0;
        this.p2Soft = p2Soft;
        this.latticeCoordinates = latticeCoordinates;
        this.positionDefinition = new MoleculePositionCOM(space);
    	Lxyz = space.makeVector();
		drtmp = space.makeVector();
		nShells = new int[] {(int) Math.ceil(rCut/a0[0] - 0.49999), (int) Math.ceil(rCut/a0[1] - 0.49999), (int) Math.ceil(rCut/a0[2] - 0.49999)};
	}
        
    public double energy(IAtomList atoms) {
    	
    	boolean isSelf = (atoms.getAtom(1) == atoms.getAtom(0));
		double u_LJ = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        r0.E(latticeCoordinates == null?positionDefinition.position(atoms.getAtom(0).getParentGroup()):
        	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atoms.getAtom(0).getParentGroup())).position);
        Vector r1 = latticeCoordinates == null?positionDefinition.position(atoms.getAtom(1).getParentGroup()):
        	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atoms.getAtom(1).getParentGroup())).position;
        
        ldr.Ev1Mv2(r1, r0);
        shift.Ea1Tv1(-1, ldr);
        boundary.nearestImage(ldr);
        shift.PE(ldr);
		dr.PE(shift);
        
        
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	Lxyz.setX(2, nz*a0[2]);
					drtmp.Ev1Pv2(ldr, Lxyz);
					double ldr2 = drtmp.squared();
					if(ldr2 > rCut2 ) continue;
                	boolean centerImage = (nx*nx+ny*ny+nz*nz == 0);
                	if(isSelf && centerImage) continue;
                	drtmp.Ev1Pv2(dr, Lxyz);
					double dr2 = drtmp.squared();
					
//					if(Math.abs((isSelf ? 0.5 : 1.0)*p2Soft.u(dr2) )< 0.0001){
//						System.out.println(p2Soft.u(dr2));
//						throw new RuntimeException();
//					}
//					
//					System.out.println(atoms.getAtom(0).getLeafIndex() + " " +atoms.getAtom(1).getLeafIndex() + " " + (isSelf ? 0.5 : 1.0)*p2Soft.u(dr2));
					
					
                	u_LJ += (isSelf ? 0.5 : 1.0)*p2Soft.u(dr2);
                }
            }
        }
        return u_LJ;
    }
    
    public static boolean debug = false;
    
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(IAtomList atoms) {
        double tmpVir = 0;
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        boundary.nearestImage(dr);
		for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	Lxyz.setX(2, nz*a0[2]);
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
    public Vector[] gradient(IAtomList atoms) {
    	boolean isSelf = (atoms.getAtom(1) == atoms.getAtom(0));
        dr.Ev1Mv2(atoms.getAtom(1).getPosition(),atoms.getAtom(0).getPosition());
        r0.E(latticeCoordinates == null?positionDefinition.position(atoms.getAtom(0).getParentGroup()):
        	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atoms.getAtom(0).getParentGroup())).position);
        Vector r1 = latticeCoordinates == null?positionDefinition.position(atoms.getAtom(1).getParentGroup()):
        	((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(atoms.getAtom(1).getParentGroup())).position;
     
        
        ldr.Ev1Mv2(r1, r0);
        shift.Ea1Tv1(-1, ldr);
        boundary.nearestImage(ldr);
        shift.PE(ldr);
        dr.PE(shift);
   
        
        gradient[0].E(0);
        gradient[1].E(0);
        for(int nx = -nShells[0]; nx <= nShells[0]; nx++) {
        	Lxyz.setX(0, nx*a0[0]);
            for(int ny = -nShells[1]; ny <= nShells[1]; ny++) {
            	Lxyz.setX(1, ny*a0[1]);
                for(int nz = -nShells[2]; nz <= nShells[2]; nz++) {
                	boolean nonCenterImage = (nx*nx+ny*ny+nz*nz > 0);
                	if(isSelf && nonCenterImage) continue;
                	Lxyz.setX(2, nz*a0[2]);
                	drtmp.Ev1Pv2(ldr, Lxyz);
					double ldr2 = drtmp.squared();
					if(ldr2 > rCut2 ) continue;
					drtmp.Ev1Pv2(dr, Lxyz);
					double dr2 = drtmp.squared();
			        gradient[1].PEa1Tv1(p2Soft.du(dr2)/dr2,drtmp);
                }
            }
        }
        gradient[0].PEa1Tv1(-1,gradient[1]);
        return gradient;
    }
    
    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
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

    public void setBox(Box box) {
        boundary = box.getBoundary();
        p2Soft.setBox(box);
    }

    protected final Vector[] gradient;
    protected Boundary boundary;
    protected final int[] nShells;
    protected final double[] a0;
    protected final Potential2Soft p2Soft;
    protected final Vector Lxyz;
    protected final Vector dr;
    protected final Vector ldr;
    protected final Vector shift;
    protected final Vector drtmp;
    protected final double rCut2;
    protected final MoleculeAgentManager latticeCoordinates;
    protected IMoleculePositionDefinition positionDefinition;
    protected final Vector r0;

}//end of Potential2SoftSpherical
