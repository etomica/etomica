
package etomica.models.water;

import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.*;

public class P2WaterSPC extends Potential2 implements Potential2.Soft {

	public P2WaterSPC(Simulation sim, PotentialTruncation potentialTruncation, Space3D.Boundary boundary) {
		this(sim, potentialTruncation);
		this.boundary = boundary;
	}
	public P2WaterSPC(Simulation sim, PotentialTruncation potentialTruncation) {
		super(Simulation.instance.hamiltonian.potential);
		setSigma(3.1670);
		setEpsilon(Kelvin.UNIT.toSim(78.23));
		this.potentialTruncation = potentialTruncation;
		this.space = (Space3D)sim.space();
		work = (Space3D.Vector)space.makeVector();
		setCharges();
	}   
	public double energy(AtomPair pair){
		double distance  = pair.r2();
		if(potentialTruncation.isZero(distance)) return 0.0;
		else{
			AtomTreeNodeWater node1 = (AtomTreeNodeWater)pair.atom1().node;
			AtomTreeNodeWater node2 = (AtomTreeNodeWater)pair.atom2().node;
			
			Space3D.Vector O1r = (Space3D.Vector)node1.O.coord.position();
			Space3D.Vector O2r = (Space3D.Vector)node2.O.coord.position();
			Space3D.Vector H11r = (Space3D.Vector)node1.H1.coord.position();
			Space3D.Vector H12r = (Space3D.Vector)node1.H2.coord.position();
			Space3D.Vector H21r = (Space3D.Vector)node2.H1.coord.position();
			Space3D.Vector H22r = (Space3D.Vector)node2.H2.coord.position();

            
			double sum = 0.0;
			double r2 = 0.0;
			
			r2 = Space3D.r2(O1r, O2r, boundary, work);
			if(r2<1.6) {
				return Double.POSITIVE_INFINITY;
			}	
			sum += chargeOO/Math.sqrt(r2);
			double s2 = sigma2/(r2);
			double s6 = s2*s2*s2;
			sum += epsilon4*s6*(s6 - 1.0);
			
			r2 = Space3D.r2(O1r, H21r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeOH/Math.sqrt(r2);
			r2 = Space3D.r2(O1r, H22r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeOH/Math.sqrt(r2);
			r2 = Space3D.r2(H11r, O2r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeOH/Math.sqrt(r2);
			r2 = Space3D.r2(H11r, H21r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeHH/Math.sqrt(r2);
			r2 = Space3D.r2(H11r, H22r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeHH/Math.sqrt(r2);
			r2 = Space3D.r2(H12r, O2r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeOH/Math.sqrt(r2);
			r2 = Space3D.r2(H12r, H21r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeHH/Math.sqrt(r2);
			r2 = Space3D.r2(H12r, H22r, boundary, work);
						if(r2<1.6) {
							return Double.POSITIVE_INFINITY;
						}
			sum += chargeHH/Math.sqrt(r2);
			return sum;																					
		}
        
	}//end of energy
    
	public Space.Vector gradient(AtomPair pair){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double hyperVirial(AtomPair pair){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double integral(double rC){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double virial(AtomPair pair){
		throw new etomica.exception.MethodNotImplementedException();
	}
    
	public double getSigma() {return sigma;}
    
	private final void setSigma(double s) {
		sigma = s;
		sigma2 = s*s;
	}
    
	public double getEpsilon() {return epsilon;}
    
	private final void setEpsilon(double eps) {
		epsilon = eps;
		epsilon4 = 4*epsilon;
	}
	private final void setCharges() {
		chargeOO = chargeO * chargeO;
		chargeOH = chargeO * chargeH;
		chargeHH = chargeH * chargeH;
	}
    
	public double sigma , sigma2;
	public double epsilon, epsilon4;
	private PotentialTruncation potentialTruncation;
	private Atom O1, H11, H12, O2, H21, H22;
	private Space3D.Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.41);
	private double chargeO = Electron.UNIT.toSim(-0.82);
	private double chargeOO, chargeOH, chargeHH;
	private Space3D space;
	private Space3D.Vector work;
	/**
	 * Returns the boundary.
	 * @return Space3D.Boundary
	 */
	public Space3D.Boundary getBoundary() {
		return boundary;
	}

	/**
	 * Sets the boundary.
	 * @param boundary The boundary to set
	 */
	public void setBoundary(Space3D.Boundary boundary) {
		this.boundary = boundary;
	}

}