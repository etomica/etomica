package etomica.chem.models.water;

import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialTruncation;
import etomica.space.Vector;
import etomica.space3d.Boundary;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.*;

public class P2WaterTIP4P extends Potential2 implements Potential2Soft {

	public P2WaterTIP4P(Space space, PotentialTruncation potentialTruncation, Boundary boundary) {
		this(space, potentialTruncation);
		this.boundary = boundary;
	}
	public P2WaterTIP4P(Space space, PotentialTruncation potentialTruncation) {
		super(space, potentialTruncation);
		setSigma(3.15365);
		setEpsilon(Kelvin.UNIT.toSim(77.94));
		this.potentialTruncation = potentialTruncation;
		work = (Vector)space.makeVector();
		shift = (Vector)space.makeVector();
		setCharges();
	}   
	public double energy(Atom[] pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomTreeNodeTIP4PWater node1 = (AtomTreeNodeTIP4PWater)pair[0].node;
		AtomTreeNodeTIP4PWater node2 = (AtomTreeNodeTIP4PWater)pair[1].node;
		
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
		
		
		//work.Ev1Mv2(Charge1r, Charge2r);
		//boundary.nearestImage(work, shift);
		//r2 = work.squared();

		

		
		
		Vector Charge1r = (Vector)node1.Charge.coord.position();
		Vector Charge2r = (Vector)node2.Charge.coord.position();
		
		r2 = (zeroShift) ? Charge1r.Mv1Squared(Charge2r) : Charge1r.Mv1Pv2Squared(Charge2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
	    sum += chargeChargeCharge/Math.sqrt(r2); 
		
		
		
		Vector H11r = (Vector)node1.H1.coord.position();
		Vector H12r = (Vector)node1.H2.coord.position();
		Vector H21r = (Vector)node2.H1.coord.position();
		Vector H22r = (Vector)node2.H2.coord.position();
        		
		
					
		r2 = (zeroShift) ? Charge1r.Mv1Squared(H21r) : Charge1r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);
		
		r2 = (zeroShift) ? Charge1r.Mv1Squared(H22r) : Charge1r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(Charge2r) : H11r.Mv1Pv2Squared(Charge2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(H21r) : H11r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H11r.Mv1Squared(H22r) : H11r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(Charge2r) : H12r.Mv1Pv2Squared(Charge2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeChargeH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(H21r) : H12r.Mv1Pv2Squared(H21r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = (zeroShift) ? H12r.Mv1Squared(H22r) : H12r.Mv1Pv2Squared(H22r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

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
	private PotentialTruncation potentialTruncation;
	private Atom O1, H11, H12, O2, H21, H22, Charge1, Charge2;
	private Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.52);
	private double chargeCharge = Electron.UNIT.toSim(-1.04);
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