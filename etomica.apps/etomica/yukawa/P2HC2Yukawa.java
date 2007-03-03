package etomica.yukawa;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.Potential2SoftSpherical;
import etomica.simulation.Simulation;
import etomica.space.IVectorRandom;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;

/**
 * Hard-core plus two Yukawa fluid (HC2Yukawa): A Lennard-Jones like potential.
 * 
 * ==================================================================================================================================
 * 2001. Pini, Stell, and Wilding. "Liquid-gas phase behaviour of an argon-like fluid modelled by the hard-core two-Yukawa potential"
 * ==================================================================================================================================
 * 	
 * 			| infinity																		r <= sigma
 * U(r) =	| 
 * 			|(A1/r) * exp[-z1 * (r - sigma)] - (A2 * epsilon/r) * exp[-z2 * (r - sigma)]	r > sigma
 * 
 * 	LJ behavior parameters:	A1 = 1.6438 * sigma
 * 							z1 = 14.7 / sigma
 * 							A2 = 2.03 * sigma
 * 							z2 = 2.69 / sigma
 * 
 * @author msellers
 */

public final class P2HC2Yukawa extends Potential2SoftSpherical implements EtomicaElement {


	
	public double d2u(double r2) {
		// TODO Auto-generated method stub
		throw new RuntimeException();
	}

	public double du(double r2) {
		// TODO Auto-generated method stub
		throw new RuntimeException();
	}

	public P2HC2Yukawa(Simulation sim){
		this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
	}
	
	public P2HC2Yukawa(Space space, double sigma, double epsilon){
		super(space);
		
		dr = space.makeVector();
		setSigma(sigma);
		setEpsilon(epsilon);
		setParameters(sigma);
	}
	
	public double getRange(){return Double.POSITIVE_INFINITY;}
	
	public void setPhase(Phase phase) {
		nearestImageTransformer = phase.getBoundary();
	}
	
	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo("Hard-core plus two Yukawa fluid potential");
		return info;
	}
	
	/**
	 * Energy method.  u(double r) form.
	 */
	public double u(double r2){
		double r = Math.sqrt(r2);
		//hard core repulsion
		if (r <= sigma){
			return Double.POSITIVE_INFINITY;
		}
		
		//two-Yukawa tail attraction
		double expRsigma = Math.exp(r - sigma);
		double uterm1 = (A1 / r) * expZ1 * expRsigma;
		double uterm2 = (A2 * epsilon / r) * expZ2 * expRsigma;
		
		return (uterm1 - uterm2);
	}
		
    /**
     * Energy of the pair as given by the u(double) method
     */
    public double energy(AtomSet atoms) {
        AtomPair pair = (AtomPair)atoms;
        dr.Ev1Mv2(((AtomLeaf)pair.atom1).getCoord().getPosition(),((AtomLeaf)pair.atom0).getCoord().getPosition());
        nearestImageTransformer.nearestImage(dr);
        return u(dr.squared());
    }
	
    /**
     * Integral from rC to infinity.
     */
	public double uInt(double rC){
		
		return 0;
	}

	/**
	 * Accessor methods for size and energy parameters.
	 */
	public double getSigma() {return sigma;}
	
	public double getEpsilon() {return epsilon;}
	
	/**
	 * Mutator methods for size and energy parameters.
	 */
	public final void setEpsilon(double eps) {epsilon = eps;}
	
	public final void setSigma(double s) {sigma = s;}
	
	/**
	 * Parameter calculation method.
	 */
	public final void setParameters(double s){
		A1 = 1.6438 * s;
		A2 = 2.03 * s;
		z1 = 14.7 / s;
		z2 = 2.69 / s;
		expZ1 = Math.exp(-z1);
		expZ2 = Math.exp(-z2);
	}
	
    private static final long serialVersionUID = 1L;
	private double sigma;
	private double epsilon;
	private double A1;
	private double A2;
	private double z1;
	private double z2;
	private double expZ1;
	private double expZ2; 
	private final IVectorRandom dr;
	private NearestImageTransformer nearestImageTransformer;
}	
