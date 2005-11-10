package etomica.potential;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.EmbeddedAtomMethodP2.Wrapper;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 * This class defines the many-body term of the Embedded-Atom Method (EAM) interatomic 
 * potential as given by Baskes (1992).
 * 
 * It was created for use with the EmbeddedAtomMethodPInitial and EmbeddedAtomMethodPMany
 * classes in the EAMMd3D simulation class by A. Schultz and K.R. Schadel July 2005.
 * 
 * The potential is treated as one-body potential.  The pair-wise terms required for
 * its calculation are all stored in the allatomAgents array.  The potential can be
 * treated as one-body because to calculate it, one only needs to call the information
 * stored for that particular atom in the array.
 */

public final class EmbeddedAtomMethodPMany extends Potential1 {

	public EmbeddedAtomMethodPMany(Space space, ParameterSetEAM p, EmbeddedAtomMethodP2 eamP2) {
		super(space);
        this.p = p;
        gradient = space.makeVector(); //initializes the vector "gradient"
        agentSource = eamP2;
    }

    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}
    private ParameterSetEAM p;
    private double r2Last = -1.0;
    
    private final Vector gradient;
	private double radius;
	
    public void setPhase(Phase phase) {
        super.setPhase(phase);
        agents = agentSource.getAgents(phase);
    }
    
	/**
     * The following field returns the EAM model's many-body-potential term, "energy",
     * which approximates the energy required to embedd atom i in the surrounding 
     * electron cloud.  Only the electron densities of the nearest neighbors, which are 
     * approximated as atomic densities, are considered.  Just as for the pairwise
     * term of the potential, the only variable in the many-body energy is the distance 
     * between i and j.  The parameters, which are dependent upon the type of metallic
     * system to be modeled, are provided in the class ParameterSetEAM.
     * 
     * Rho, the atomic electron density, also appears in the EAM model's pairwise 
     * potential term. The sum over rho for all of atom i's nearest neighbors is 
     * calculated in EmbeddedAtomMethodP2.  The embedding energy, F, 
     * for each atom i equals A*E/Z*rohSummed(ln(rohSummed)-ln(Z)).
     */
	
	public double rhoSummed(AtomSet a) {
		return agents[((Atom)a).getGlobalIndex()].x;
	}
	
	public double energy(AtomSet a) {
		double rhoSummed = rhoSummed(a);
		return p.A*p.E/p.Z*rhoSummed*
			(Math.log(rhoSummed)-Math.log(p.Z));
	}
	
	public Vector A(AtomSet a) {
		return agents[((Atom)a).getGlobalIndex()].A;
	}

	public Vector gradient(AtomSet a) {
		double rhoSummed = rhoSummed(a);
		double denergy = p.A*p.E/p.Z*(1 + Math.log(rhoSummed/p.Z));
		gradient.Ea1Tv1(denergy, A(a));
		/** 		System.out.println("In EmbeddedAtomMethod.PMany, gradient is " +
        		gradient + 
				", energy is " +
				energy(a));
         		", vector work1 is < " + 
        		ElectronVolt.UNIT.fromSim(work1.x(0)) + ", " +
        		ElectronVolt.UNIT.fromSim(work1.x(1)) + ", " +
        		ElectronVolt.UNIT.fromSim(work1.x(2)) + 
        		">,  and the interatomic distance is " + r)
        System.out.println("In EmbeddedAtomMethod.PMany, " +
        		p.E +
				", rhoSummed is " +
				rhoSummed(a));
				System.exit(0);*/
		return gradient;
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
    
	/**
	 * Returns the radius.
	 * @return double
	 */
	public double getRadius() {
		return radius;
	}

	/**
	 * Sets the radius.
	 * @param radius The radius to set
	 */
	public void setRadius(double radius) {
		this.radius = radius;
	}
	
    private final EmbeddedAtomMethodP2 agentSource;
    private Wrapper[] agents;
}
