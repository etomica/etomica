
package etomica.models.water;

import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialTruncation;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.*;

public class P2WaterSPC extends Potential2 implements Potential2Soft {

	public P2WaterSPC(Space space, PotentialTruncation potentialTruncation, etomica.space3d.Boundary boundary) {
		this(space, potentialTruncation);
		this.boundary = boundary;
	}
	public P2WaterSPC(Space space, PotentialTruncation potentialTruncation) {
		super(space, potentialTruncation);
		setSigma(3.1670);
		setEpsilon(Kelvin.UNIT.toSim(78.23));
		work = (etomica.space3d.Vector)space.makeVector();
		shift = (etomica.space3d.Vector)space.makeVector();
		setCharges();
	}   
	public double energy(Atom[] pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomTreeNodeWater node1 = (AtomTreeNodeWater)pair[0].node;
		AtomTreeNodeWater node2 = (AtomTreeNodeWater)pair[1].node;
		
		//compute O-O distance to consider truncation	
		etomica.space3d.Vector O1r = (etomica.space3d.Vector)node1.O.coord.position();
		etomica.space3d.Vector O2r = (etomica.space3d.Vector)node2.O.coord.position();

		work.Ev1Mv2(O1r, O2r);
		boundary.nearestImage(work, shift);
		r2 = work.squared();

		if(potentialTruncation.isZero(r2)) return 0.0;

		if(r2<1.6) return Double.POSITIVE_INFINITY;
	
		sum += chargeOO/Math.sqrt(r2);
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		
		etomica.space3d.Vector H11r = (etomica.space3d.Vector)node1.H1.coord.position();
		etomica.space3d.Vector H12r = (etomica.space3d.Vector)node1.H2.coord.position();
		etomica.space3d.Vector H21r = (etomica.space3d.Vector)node2.H1.coord.position();
		etomica.space3d.Vector H22r = (etomica.space3d.Vector)node2.H2.coord.position();
        		
		final boolean zeroShift = shift.isZero();
					
		r2 = (zeroShift) ? O1r.Mv1Squared(H21r) : O1r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);
		
		r2 = (zeroShift) ? O1r.Mv1Squared(H22r) : O1r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(O2r) : H11r.Mv1Pv2Squared(O2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(H21r) : H11r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(H22r) : H11r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(O2r) : H12r.Mv1Pv2Squared(O2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(H21r) : H12r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(H22r) : H12r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		return sum;																					        
	}//end of energy
    
	public etomica.space.Vector gradient(Atom[] pair){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double hyperVirial(Atom[] pair){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double integral(double rC){
		throw new etomica.exception.MethodNotImplementedException();
	}
	public double virial(Atom[] pair){
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
	private etomica.space3d.Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.41);
	private double chargeO = Electron.UNIT.toSim(-0.82);
	private double chargeOO, chargeOH, chargeHH;
	private etomica.space3d.Vector work, shift;
	/**
	 * Returns the boundary.
	 * @return Space3D.Boundary
	 */
	public etomica.space3d.Boundary getBoundary() {
		return boundary;
	}

	/**
	 * Sets the boundary.
	 * @param boundary The boundary to set
	 */
	public void setBoundary(etomica.space3d.Boundary boundary) {
		this.boundary = boundary;
	}

}
