package etomica.meam;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.meam.MEAMP2.Wrapper;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.potential.Potential1;
import etomica.potential.PotentialSoft;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * This class calculates the MEAM potential of each atom i, using the values that 
 * have been summed over all relevant neighbors, j, for each atom i in its bins.
 * 
 * Note that this class does yet include the screening function.  For the time
 * being, we will use a simple nearest-neighbor cutoff function only.
 * 
 * This class was created by A. Schultz and K.R. Schadel July 2005 as part of
 * a pseudo embedded-atom method potential.  In February 2006 it was adapted 
 * to be a part of a modified embedded-atom method potential.
 */

public final class MEAMPMany extends Potential1 implements PotentialSoft {

	public MEAMPMany(Space space, ParameterSetMEAM p, PhaseAgentManager phaseAgentManager) {
		super(space);
        this.p = p;
        this.phaseAgentManager = phaseAgentManager;
        gradRhoi0 = (Vector3D)space.makeVector();
        gradRhoi1 = (Vector3D)space.makeVector();
        gradRhoi2 = (Vector3D)space.makeVector();
        gradRhoi3 = (Vector3D)space.makeVector();
        gradtav1 = (Vector3D)space.makeVector();
        gradtav2 = (Vector3D)space.makeVector();
        gradtav3 = (Vector3D)space.makeVector();
        gradGamma = (Vector3D)space.makeVector();
        gradRhoi = (Vector3D)space.makeVector();
        gradF = (Vector3D)space.makeVector();
        gradEi[0] = (Vector3D)space.makeVector();
    }
	
    public void setPhase(Phase phase) {
        super.setPhase(phase);
        AtomAgentManager[] agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        agents = (Wrapper[])agentManager[phase.getIndex()].getAgents();
    }
    
    protected double rhoi0(Wrapper agent) {
    	return agent.sums[Wrapper.rhoj0Bin];
    }
    
    protected double rhoi1(Wrapper agent) {
		return Math.sqrt( 
				(agent.sums[Wrapper.rhoj1xBin] * agent.sums[Wrapper.rhoj1xBin])
				+ (agent.sums[Wrapper.rhoj1yBin] * agent.sums[Wrapper.rhoj1yBin])
				+ (agent.sums[Wrapper.rhoj1zBin] * agent.sums[Wrapper.rhoj1zBin]));
    }
    
    protected double rhoi2(Wrapper agent) {
		return Math.sqrt(   
			(agent.sums[Wrapper.rhoj2xxBin]* agent.sums[Wrapper.rhoj2xxBin])
			+ (2.0 * agent.sums[Wrapper.rhoj2xyBin] * agent.sums[Wrapper.rhoj2xyBin])
			+ (2.0 * agent.sums[Wrapper.rhoj2xzBin] * agent.sums[Wrapper.rhoj2xzBin])
			+ (agent.sums[Wrapper.rhoj2yyBin] * agent.sums[Wrapper.rhoj2yyBin])
			+ (2.0 * agent.sums[Wrapper.rhoj2yzBin] * agent.sums[Wrapper.rhoj2yzBin])
			+ (agent.sums[Wrapper.rhoj2zzBin] * agent.sums[Wrapper.rhoj2zzBin])
			- ((1.0/3.0) * agent.sums[Wrapper.rhoj2Bin] * agent.sums[Wrapper.rhoj2Bin]));
    }
    
    protected double rhoi3(Wrapper agent) {
		return Math.sqrt(
		     (agent.sums[Wrapper.rhoj3xxxBin] * agent.sums[Wrapper.rhoj3xxxBin])
		+(3.0 * agent.sums[Wrapper.rhoj3xxyBin] * agent.sums[Wrapper.rhoj3xxyBin])
		+(3.0 * agent.sums[Wrapper.rhoj3xxzBin] * agent.sums[Wrapper.rhoj3xxzBin])
		+(3.0 * agent.sums[Wrapper.rhoj3xyyBin] * agent.sums[Wrapper.rhoj3xyyBin])
		+(6.0 * agent.sums[Wrapper.rhoj3xyzBin] * agent.sums[Wrapper.rhoj3xyzBin])
		+(3.0 * agent.sums[Wrapper.rhoj3xzzBin] * agent.sums[Wrapper.rhoj3xzzBin])
		+    (agent.sums[Wrapper.rhoj3yyyBin] * agent.sums[Wrapper.rhoj3yyyBin])
		+(3.0 * agent.sums[Wrapper.rhoj3yyzBin] * agent.sums[Wrapper.rhoj3yyzBin])
		+(3.0 * agent.sums[Wrapper.rhoj3yzzBin] * agent.sums[Wrapper.rhoj3yzzBin])
		+    (agent.sums[Wrapper.rhoj3zzzBin] * agent.sums[Wrapper.rhoj3zzzBin]));
    }
    
    protected double tav1(Wrapper agent) {
    	return agent.sums[Wrapper.t1Rhoj0Bin] / agent.sums[Wrapper.rhoj0Bin];
    }

    protected double tav2(Wrapper agent) {
    	return agent.sums[Wrapper.t2Rhoj0Bin] / agent.sums[Wrapper.rhoj0Bin];
    }
    
    protected double tav3(Wrapper agent) {
    	return agent.sums[Wrapper.t3Rhoj0Bin] / agent.sums[Wrapper.rhoj0Bin];
    }
    
    protected double gamma(Wrapper agent) {
    	double rhoi0 = rhoi0(agent);
		double rhoi1 = rhoi1(agent);
		double rhoi2 = rhoi2(agent);
		double rhoi3 = rhoi3(agent);
		double tav1 = tav1(agent);
		double tav2 = tav2(agent);
		double tav3 = tav3(agent);
    	return (tav1 * (rhoi1/rhoi0) * (rhoi1/rhoi0))
		+ (tav2 * (rhoi2/rhoi0) * (rhoi2/rhoi0))
		+ (tav3 * (rhoi3/rhoi0) * (rhoi3/rhoi0));
    }
    
    protected double rhoi(Wrapper agent) {
    	double rhoi0 = rhoi0(agent);
		double gamma = gamma(agent);
		
		
		
		//The following expression for the background electron density of atom i
		//is appropriate for Sn only.
		return (2.0 * rhoi0) / (1.0 + Math.exp(-gamma));
    }
    
	public double energy(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		
		double rhoi = rhoi(agent);
		
		//The embedding energy of atom i
		double F = p.A * p.Ec * (rhoi/p.rhoScale) * Math.log(rhoi/p.rhoScale);
		
		return F + ((1.0/2.0)*agent.sums[Wrapper.phiBin]);
	}

	public Vector[] gradient(AtomSet a) {
		Wrapper agent = agents[((Atom)a).getGlobalIndex()];
		
		double rhoi0 = rhoi0(agent);
		double rhoi1 = rhoi1(agent);
		double rhoi2 = rhoi2(agent);
		double rhoi3 = rhoi3(agent);
		double tav1 = tav1(agent);
		double tav2 = tav2(agent);
		double tav3 = tav3(agent);
		double gamma = gamma(agent);
		double rhoi = rhoi(agent);
		
		
		gradRhoi0.E(agent.gradientSums[Wrapper.rhoj0Bin]);
		
		gradRhoi1.Ea1Tv1(agent.sums[Wrapper.rhoj1xBin], 
				agent.gradientSums[Wrapper.rhoj1xBin]);
		gradRhoi1.PEa1Tv1(agent.sums[Wrapper.rhoj1yBin], 
				agent.gradientSums[Wrapper.rhoj1yBin]);
		gradRhoi1.PEa1Tv1(agent.sums[Wrapper.rhoj1zBin], 
				agent.gradientSums[Wrapper.rhoj1zBin]);
		gradRhoi1.TE(1.0/rhoi1);
		
		gradRhoi2.Ea1Tv1(agent.sums[Wrapper.rhoj2xxBin], agent.gradientSums[Wrapper.rhoj2xxBin]);
		gradRhoi2.PEa1Tv1(2.0 * agent.sums[Wrapper.rhoj2xyBin], agent.gradientSums[Wrapper.rhoj2xyBin]);
		gradRhoi2.PEa1Tv1(2.0 * agent.sums[Wrapper.rhoj2xzBin], agent.gradientSums[Wrapper.rhoj2xzBin]);
		gradRhoi2.PEa1Tv1(agent.sums[Wrapper.rhoj2yyBin], agent.gradientSums[Wrapper.rhoj2yyBin]);
		gradRhoi2.PEa1Tv1(2.0 * agent.sums[Wrapper.rhoj2yzBin], agent.gradientSums[Wrapper.rhoj2yzBin]);
		gradRhoi2.PEa1Tv1(agent.sums[Wrapper.rhoj2zzBin], agent.gradientSums[Wrapper.rhoj2zzBin]);
		gradRhoi2.PEa1Tv1(-(1.0/3.0) * agent.sums[Wrapper.rhoj2Bin], agent.gradientSums[Wrapper.rhoj2Bin]);
		gradRhoi2.TE(1.0/rhoi2);
		
		gradRhoi3.Ea1Tv1(agent.sums[Wrapper.rhoj3xxxBin], 
						agent.gradientSums[Wrapper.rhoj3xxxBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3xxyBin], 
					agent.gradientSums[Wrapper.rhoj3xxyBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3xxzBin], 
					agent.gradientSums[Wrapper.rhoj3xxzBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3xyyBin], 
					agent.gradientSums[Wrapper.rhoj3xyyBin]);
		gradRhoi3.PEa1Tv1(6.0 * agent.sums[Wrapper.rhoj3xyzBin], 
					agent.gradientSums[Wrapper.rhoj3xyzBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3xzzBin], 
					agent.gradientSums[Wrapper.rhoj3xzzBin]);
		gradRhoi3.PEa1Tv1(agent.sums[Wrapper.rhoj3yyyBin], 
					agent.gradientSums[Wrapper.rhoj3yyyBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3yyzBin], 
					agent.gradientSums[Wrapper.rhoj3yyzBin]);
		gradRhoi3.PEa1Tv1(3.0 * agent.sums[Wrapper.rhoj3yzzBin], 
					agent.gradientSums[Wrapper.rhoj3yzzBin]);
		gradRhoi3.PEa1Tv1(agent.sums[Wrapper.rhoj3zzzBin], 
					agent.gradientSums[Wrapper.rhoj3zzzBin]);
		gradRhoi3.TE(1.0/rhoi3);
		
		gradtav1.Ea1Tv1(-agent.sums[Wrapper.t1Rhoj0Bin]/(rhoi0*rhoi0), agent.gradientSums[Wrapper.rhoj0Bin]);
		gradtav1.PEa1Tv1(1.0/rhoi0, agent.gradientSums[Wrapper.t1Rhoj0Bin]);
		
		gradtav2.Ea1Tv1(-agent.sums[Wrapper.t2Rhoj0Bin]/(rhoi0*rhoi0), agent.gradientSums[Wrapper.rhoj0Bin]);
		gradtav2.PEa1Tv1(1.0/rhoi0, agent.gradientSums[Wrapper.t2Rhoj0Bin]);
		
		gradtav3.Ea1Tv1(-agent.sums[Wrapper.t3Rhoj0Bin]/(rhoi0*rhoi0), agent.gradientSums[Wrapper.rhoj0Bin]);
		gradtav3.PEa1Tv1(1.0/rhoi0, agent.gradientSums[Wrapper.t3Rhoj0Bin]);
		
		gradGamma.Ea1Tv1( (-1.0/rhoi0)*( (tav1*rhoi1*rhoi1) + (tav2*rhoi2*rhoi2)
				+ (tav3*rhoi3*rhoi3) ), gradRhoi0);
		gradGamma.PEa1Tv1(tav1*rhoi1, gradRhoi1);
		gradGamma.PEa1Tv1(tav2*rhoi2, gradRhoi2);
		gradGamma.PEa1Tv1(tav3*rhoi3, gradRhoi3);
		gradGamma.TE(2.0);
		gradGamma.PEa1Tv1(rhoi1*rhoi1, gradtav1);
		gradGamma.PEa1Tv1(rhoi2*rhoi2, gradtav2);
		gradGamma.PEa1Tv1(rhoi3*rhoi3, gradtav3);
		gradGamma.TE(1.0/(rhoi0*rhoi0));
		
		gradRhoi.Ea1Tv1(rhoi0*Math.exp(-gamma)/(1.0 + Math.exp(-gamma)), gradGamma);
		gradRhoi.PE(gradRhoi0);
		gradRhoi.TE(2.0/(1.0+Math.exp(-gamma)));
	
		gradF.Ea1Tv1( (p.A*p.Ec/p.rhoScale)*
				(1.0 + Math.log(rhoi)- Math.log(p.rhoScale)), gradRhoi);
		
		//System.out.println("gradF is " + gradF);
		//System.exit(0);
		
		gradEi[0].E(gradF);
		gradEi[0].PEa1Tv1(0.5, agent.gradientSums[Wrapper.phiBin]);
		
		return gradEi;
		
	}
	
	public double virial(AtomSet atoms) {
	    return 0.0;
    }
	
    private Wrapper[] agents;
    private ParameterSetMEAM p;
    private final Vector3D gradRhoi0;
    private final Vector3D gradRhoi1;
    private final Vector3D gradRhoi2;
    private final Vector3D gradRhoi3;
    private final Vector3D gradtav1;
    private final Vector3D gradtav2;
    private final Vector3D gradtav3;
    private final Vector3D gradGamma;
    private final Vector3D gradRhoi;
    private final Vector3D gradF;
    private final Vector3D[] gradEi = new Vector3D[1];;
	private PhaseAgentManager phaseAgentManager;
}
