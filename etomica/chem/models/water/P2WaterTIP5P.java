package etomica.chem.models.water;

import etomica.Atom;
import etomica.Space;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialTruncation;
import etomica.space.Vector;
import etomica.space3d.Boundary;
import etomica.space3d.Space3D;
import etomica.units.Electron;
import etomica.units.Kelvin;

public class P2WaterTIP5P extends Potential2 implements Potential2Soft {

	public P2WaterTIP5P(Space space, PotentialTruncation potentialTruncation, Boundary boundary) {
		this(space, potentialTruncation);
		this.boundary = boundary;
	}
	public P2WaterTIP5P(Space space, PotentialTruncation potentialTruncation) {
		super(space, potentialTruncation);
		setSigma(3.12);
		setEpsilon(Kelvin.UNIT.toSim(80.51));
		work = (Vector)space.makeVector();
		shift = (Vector)space.makeVector();
		setCharges();
	}   
	public double energy(Atom[] pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomTreeNodeTIP5PWater node1 = (AtomTreeNodeTIP5PWater)pair[0].node;
		AtomTreeNodeTIP5PWater node2 = (AtomTreeNodeTIP5PWater)pair[1].node;
		
		//compute O-O distance to consider truncation	
		Vector O1r = (Vector)node1.O.coord.position();
		Vector O2r = (Vector)node2.O.coord.position();

		work.Ev1Mv2(O1r, O2r);
		boundary.nearestImage(work, shift);
		r2 = work.squared();

		if(potentialTruncation.isZero(r2)) return 0.0;

		if(r2<1.6) return Double.POSITIVE_INFINITY;
	
		/*sum += chargeOO/Math.sqrt(r2); */
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		final boolean zeroShift = shift.isZero();
			
		
		
		Vector Charge11r = (Vector)node1.Charge1.coord.position();
		Vector Charge12r = (Vector)node1.Charge2.coord.position();
		Vector Charge21r = (Vector)node2.Charge1.coord.position();
		Vector Charge22r = (Vector)node2.Charge2.coord.position();
		Vector H11r = (Vector)node1.H1.coord.position();
		Vector H12r = (Vector)node1.H2.coord.position();
		Vector H21r = (Vector)node2.H1.coord.position();
		Vector H22r = (Vector)node2.H2.coord.position();
        		
		
					
		r2 = (zeroShift) ? H11r.Mv1Squared(H21r) : H11r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);
		
		r2 = (zeroShift) ? H11r.Mv1Squared(H22r) : H11r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(Charge21r) : H11r.Mv1Pv2Squared(Charge21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(Charge22r) : H11r.Mv1Pv2Squared(Charge22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		
		
		
		
		r2 = (zeroShift) ? H12r.Mv1Squared(H21r) : H12r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(H22r) : H12r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(Charge21r) : H12r.Mv1Pv2Squared(Charge21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(Charge22r) : H12r.Mv1Pv2Squared(Charge22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		
		
		
		
		r2 = (zeroShift) ? Charge11r.Mv1Squared(H21r) : Charge11r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);
		
		r2 = (zeroShift) ? Charge11r.Mv1Squared(H22r) : Charge11r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? Charge11r.Mv1Squared(Charge21r) : Charge11r.Mv1Pv2Squared(Charge21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeCharge/Math.sqrt(r2);

		r2 = (zeroShift) ? Charge11r.Mv1Squared(Charge22r) : Charge11r.Mv1Pv2Squared(Charge22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeCharge/Math.sqrt(r2);

				
				
				
				
		r2 = (zeroShift) ? Charge12r.Mv1Squared(H21r) : Charge12r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? Charge12r.Mv1Squared(H22r) : Charge12r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? Charge12r.Mv1Squared(Charge21r) : Charge12r.Mv1Pv2Squared(Charge21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
	    sum += chargeChargeCharge/Math.sqrt(r2);

		r2 = (zeroShift) ? Charge12r.Mv1Squared(Charge22r) : Charge12r.Mv1Pv2Squared(Charge22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeCharge/Math.sqrt(r2);




		return sum;																					        
	}//end of energy
    
	public Vector gradient(Atom[] pair){
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
		chargeChargeCharge = chargeCharge * chargeCharge;
		chargeChargeH = chargeCharge * chargeH;
		chargeHH = chargeH * chargeH;
	}
    
	public double sigma , sigma2;
	public double epsilon, epsilon4;
	private Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.241);
	private double chargeCharge = Electron.UNIT.toSim(-0.241);
	private double chargeChargeCharge, chargeChargeH, chargeHH;
	private Vector work, shift;
	/**
	 * Returns the boundary.
	 * @return Space3D.Boundary
	 */
	public Boundary getBoundary() {
		return boundary;
	}

	/**
	 * Sets the boundary.
	 * @param boundary The boundary to set
	 */
	public void setBoundary(Boundary boundary) {
		this.boundary = boundary;
	}

}
