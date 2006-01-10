package etomica.eam;
import java.awt.Color;

import etomica.EtomicaElement;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.util.Arrays;

/**
 * This class defines the pairwise term of the Embedded-Atom Method (EAM) interatomic 
 * potential as given by Baskes (1992).  It is also used to calculate the pair-wise
 * terms used in the EmbeddedAtomMethodPMany class.
 * 
 * It was created for use with the EmbeddedAtomMethodPInitial and EmbeddedAtomMethodPMany
 * classes in the EAMMd3D simulation class by A. Schultz and K.R. Schadel July 2005.
 */

public final class EmbeddedAtomMethodP2 extends Potential2SoftSpherical implements EtomicaElement, AtomAgentManager.AgentSource {

	public EmbeddedAtomMethodP2(Simulation sim, ParameterSetEAM p) {
		super(sim.space);
        this.p = p;
        work1 = space.makeVector();
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.speciesRoot);
    }
	
    public void setPhase(Phase phase){
        super.setPhase(phase);
        agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        agents = getAgents(phase);
    }
    
    public Wrapper[] getAgents(Phase phase) {
        agents = (Wrapper[])agentManager[phase.getOrdinal()-1].getAgents();
        return agents;
    }

    /**
     * The following field returns the EAM model's pair-potential term, u, which 
     * defines the interactions between atom i and each of its nearest neighbors, j.  
     * Later, the pair potential is summed over all j for each atom i.  The only 
     * variable in u is the distance between i and j.  The parameters, which are 
     * dependent upon the type of metallic system to be modeled, are provided in the 
     * class ParameterSetEAM.
     * 
     * The pair-wise potential is proportional to the difference between the energy per
     * atom in the reference structure and the embedding energy per atom in the 
     * reference structure.  This difference is then divided by the coordination 
     * number of the reference structure so that it applies to the real system.  (The
     * reference structure for Sn is FCC).
     * 
     * Eref is the energy per atom of the reference crystal structure.  It is 
     * calculated from an equation of state.
     *  
     * Fref is the pairwise contribution to the embedding energy from j per atom in the
     * reference crystal structure. 
     * 
     * Rho, the atomic electron density, also appears in the EAM model's many-body 
     * potential term. The sum over roh for all of atom i's nearest neighbors must be 
     * calculated to determine this part of the potential.  The embedding energy, F, 
     * for each atom i equals A*E/Z*rhoSummed(ln(rhoSummed)-ln(Z)).
     */
    public double u (double r2) {
    	double r = Math.sqrt(r2);
    	double a = p.alpha * ((r/p.Ro) - 1);
    	double rho = Math.exp(-p.beta * ((r/p.Ro) - 1));
    	agents[pair.atom0.getGlobalIndex()].x += rho;
        agents[pair.atom1.getGlobalIndex()].x += rho;
    	double Eref = - p.E * (1 + a) * Math.exp(- a);
    	double Fref = p.A * p.E * p.Zd/p.Z * rho *
		(Math.log(p.Zd/p.Z) - p.beta * ((r/p.Ro) - 1));
    	return (1/p.Z) * (Eref - Fref);
    }
    
    

    /**
     * The derivative du/dr multiplied by r.  Multiplying the derivative by r allows
     * du to represent energy.  
     * 
     * rho is calculated and the allatomAgents array is created again in case the 
     * gradient is required without the energy first being used; unless rho and the 
     * array are defined in the du field, they will not be accessible in such a case.
     */
    public double du(double r2) {
        double r = Math.sqrt(r2);
        double rho = Math.exp(-p.beta * ((r/p.Ro) - 1));
        
        agents[pair.atom0.getGlobalIndex()].x += rho;
        agents[pair.atom1.getGlobalIndex()].x += rho;
        agents[pair.atom0.getGlobalIndex()].A.PEa1Tv1(-rho*p.beta/p.Ro*2,coordinatePair.dr()); // kmb changed from cPair on 8/3/05
        agents[pair.atom1.getGlobalIndex()].A.PEa1Tv1(rho*p.beta/p.Ro*2,coordinatePair.dr());
        return r*(1/p.Z)*(
        		(p.E*p.alpha*p.alpha/p.Ro)*Math.exp(-p.alpha*((r/p.Ro)-1))*((r/p.Ro)-1)
				+ (p.A*p.E*p.Zd*p.beta/p.Z/p.Ro)
				  *Math.exp(-p.beta*((r/p.Ro)-1))
				  *(1+Math.log(p.Zd/p.Z)-p.beta*((r/p.Ro)-1))
			   );
        
    }

   /**
    * The second derivative of the pair energy, times the square of the
    * separation:  r^2 d^2u/dr^2.
    */
    public double d2u(double r2) {
        double r = Math.sqrt(r2);
        return r2*(1/p.Z)*(
        		(p.E*p.alpha*p.alpha/p.Ro/p.Ro)*Math.exp(-p.alpha*((r/p.Ro)-1))
					*(1-p.alpha*((r/p.Ro)-1))
				-(p.A*p.E*p.Zd*p.beta*p.beta/p.Z/p.Ro/p.Ro)
					*Math.exp(-p.beta*((r/p.Ro)-1))
					*(2 + Math.log(p.Zd/p.Z)-p.beta*((r/p.Ro)-1))
				
        );
    }
    
    public double uInt(double rC){
    	return 0; //uInt is not required for this potential.
    }
    
    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(AtomSet pair) {
    	this.pair = (AtomPair)pair;
        return super.energy(pair); 
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        return du(coordinatePair.r2());
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        double r2 = coordinatePair.r2();
        return d2u(r2) + du(r2);
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector gradient(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        double r2 = coordinatePair.r2();
        double r = Math.sqrt(r2);
        this.pair = (AtomPair)pair;
        work1.Ea1Tv1(du(r2)/r2,coordinatePair.dr());
/*        System.out.println("In EmbeddedAtomMethod.P2, du(r2)/r is " +
        		du(r2)/r + 
				", u is " +
				u(r2) +
        		", vector work1 is < " + 
        		ElectronVolt.UNIT.fromSim(work1.x(0)) + ", " +
        		ElectronVolt.UNIT.fromSim(work1.x(1)) + ", " +
        		ElectronVolt.UNIT.fromSim(work1.x(2)) +
        		">,  and the interatomic distance is " + r);*/
			//	System.exit(0); 
        return work1;
    }
    
    /**
     * Same as uInt.
     */
    public double integral(double rC) {
        return uInt(rC);
    }
    
    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    private final Vector work1;
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    private ParameterSetEAM p;
    private double r2Last = -1.0;
    private AtomPair pair;
    protected AtomAgentManager[] agentManager;
    protected Wrapper[] agents;
    private final PhaseAgentManager phaseAgentManager;
    
    public Object makeAgent(Atom atom) {
        if (atom == null) {
            return new Wrapper(null);
        }
        return new Wrapper(atom.type.isLeaf() ? (Vector)((AtomLeaf)atom).coord.position().clone() : null);
    }
    
    public void releaseAgent(Object agent) {}
    
    public static class Wrapper implements java.io.Serializable {
    	public Wrapper(Vector v) {
    		A = v;
    	}
    	public double x; //wraps double x to look like an object
    	public final Vector A; //wraps vector A to look like an object
    }
    
}
