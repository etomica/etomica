package etomica.chem.models.water;

import etomica.Atom;
import etomica.AtomSet;
import etomica.Space;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.Electron;
import etomica.units.Kelvin;

public class P2WaterTIP4P extends Potential2 implements Potential2Soft {

	public P2WaterTIP4P(Space space, Boundary boundary) {
		this(space);
		this.boundary = boundary;
	}
	public P2WaterTIP4P(Space space) {
		super(space);
		setSigma(3.15365);
		setEpsilon(Kelvin.UNIT.toSim(77.94));
		work = (Vector3D)space.makeVector();
		shift = (Vector3D)space.makeVector();
		setCharges();
	}   
	public double energy(AtomSet atoms){
		double sum = 0.0;
		double r2 = 0.0;
		
        AtomTreeNodeTIP4PWater node1 = (AtomTreeNodeTIP4PWater)atoms.getAtom(0).node;
		AtomTreeNodeTIP4PWater node2 = (AtomTreeNodeTIP4PWater)atoms.getAtom(1).node;
		
		//compute O-O distance to consider truncation	
		Vector3D O1r = (Vector3D)node1.O.coord.position();
		Vector3D O2r = (Vector3D)node2.O.coord.position();

		shift.Ev1Mv2(O1r, O2r);
        work.E(shift);
		boundary.nearestImage(shift);
        r2 = shift.squared();

		if(r2<1.6) return Double.POSITIVE_INFINITY;

        shift.ME(work);
	
		/*sum += chargeOO/Math.sqrt(r2); */
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		final boolean zeroShift = shift.isZero();
		
		
		//work.Ev1Mv2(Charge1r, Charge2r);
		//boundary.nearestImage(work, shift);
		//r2 = work.squared();

		

		
		
		Vector3D Charge1r = (Vector3D)node1.Charge.coord.position();
		Vector3D Charge2r = (Vector3D)node2.Charge.coord.position();
		
		r2 = (zeroShift) ? Charge1r.Mv1Squared(Charge2r) : Charge1r.Mv1Pv2Squared(Charge2r,shift);
		if(r2<1.6) return Double.POSITIVE_INFINITY;
	    sum += chargeChargeCharge/Math.sqrt(r2); 
		
		
		
		Vector3D H11r = (Vector3D)node1.H1.coord.position();
		Vector3D H12r = (Vector3D)node1.H2.coord.position();
		Vector3D H21r = (Vector3D)node2.H1.coord.position();
		Vector3D H22r = (Vector3D)node2.H2.coord.position();
        		
		
					
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
	private Atom O1, H11, H12, O2, H21, H22, Charge1, Charge2;
	private Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.52);
	private double chargeCharge = Electron.UNIT.toSim(-1.04);
	private double chargeChargeCharge, chargeChargeH, chargeHH;
	private Vector3D work, shift;
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