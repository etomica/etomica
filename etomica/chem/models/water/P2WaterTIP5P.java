package etomica.chem.models.water;

import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.*;

public class P2WaterTIP5P extends Potential2 implements Potential2.Soft {

	public P2WaterTIP5P(SimulationElement parent, PotentialTruncation potentialTruncation, Space3D.Boundary boundary) {
		this(parent, potentialTruncation);
		this.boundary = boundary;
	}
	public P2WaterTIP5P(SimulationElement parent, PotentialTruncation potentialTruncation) {
		super(parent);
		setSigma(3.12);
		setEpsilon(Kelvin.UNIT.toSim(80.51));
		this.potentialTruncation = potentialTruncation;
		work = (Space3D.Vector)space.makeVector();
		shift = (Space3D.Vector)space.makeVector();
		setCharges();
	}   
	public double energy(AtomPair pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomTreeNodeTIP5PWater node1 = (AtomTreeNodeTIP5PWater)pair.atom1().node;
		AtomTreeNodeTIP5PWater node2 = (AtomTreeNodeTIP5PWater)pair.atom2().node;
		
		//compute O-O distance to consider truncation	
		Space3D.Vector O1r = (Space3D.Vector)node1.O.coord.position();
		Space3D.Vector O2r = (Space3D.Vector)node2.O.coord.position();

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
			
		
		
		Space3D.Vector Charge11r = (Space3D.Vector)node1.Charge1.coord.position();
		Space3D.Vector Charge12r = (Space3D.Vector)node1.Charge2.coord.position();
		Space3D.Vector Charge21r = (Space3D.Vector)node2.Charge1.coord.position();
		Space3D.Vector Charge22r = (Space3D.Vector)node2.Charge2.coord.position();
		Space3D.Vector H11r = (Space3D.Vector)node1.H1.coord.position();
		Space3D.Vector H12r = (Space3D.Vector)node1.H2.coord.position();
		Space3D.Vector H21r = (Space3D.Vector)node2.H1.coord.position();
		Space3D.Vector H22r = (Space3D.Vector)node2.H2.coord.position();
        		
		
					
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
		chargeChargeCharge = chargeCharge * chargeCharge;
		chargeChargeH = chargeCharge * chargeH;
		chargeHH = chargeH * chargeH;
	}
    
	public double sigma , sigma2;
	public double epsilon, epsilon4;
	private PotentialTruncation potentialTruncation;
	private Atom O1, H11, H12, O2, H21, H22, Charge11, Charge12, Charge21, Charge22;
	private Space3D.Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.241);
	private double chargeCharge = Electron.UNIT.toSim(-0.241);
	private double chargeChargeCharge, chargeChargeH, chargeHH;
	private Space3D.Vector work, shift;
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
