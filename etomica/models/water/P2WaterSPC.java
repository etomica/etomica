
package etomica.models.water;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.Potential2;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * SPC potential for water.  Requires the molecule node be an
 * AtomTreeNodeWater. 
 */
public class P2WaterSPC extends Potential2 {

	public P2WaterSPC(Space space) {
		super(space);
		setSigma(3.1670);
		setEpsilon(Kelvin.UNIT.toSim(78.23));
		work = space.makeVector();
		shift = space.makeVector();
		setCharges();
	}   

    public void setPhase(Phase phase) {
        boundary = phase.getBoundary();
    }

    public double energy(AtomSet pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomWater3P water1 = (AtomWater3P)pair.getAtom(0);
		AtomWater3P water2 = (AtomWater3P)pair.getAtom(1);
		
		//compute O-O distance to consider truncation	
		IVector O1r = water1.O.getPosition();
		IVector O2r = water2.O.getPosition();

		work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
		r2 = work.squared();

		if(r2<1.6) return Double.POSITIVE_INFINITY;
	
		sum += chargeOO/Math.sqrt(r2);
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		
		IVector H11r = water1.H1.getPosition();
		IVector H12r = water1.H2.getPosition();
		IVector H21r = water2.H1.getPosition();
		IVector H22r = water2.H2.getPosition();
        		
		final boolean zeroShift = shift.squared() < 0.1;
        if (zeroShift) {
            r2 = O1r.Mv1Squared(H21r);
        }
        else {
            shift.PE(O1r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(O1r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);
		
        if (zeroShift) {
            r2 = O1r.Mv1Squared(H22r);
        }
        else {
            shift.PE(O1r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(O1r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H11r.Mv1Squared(O2r);
        }
        else {
            shift.PE(H11r);
            r2 = O2r.Mv1Squared(shift);
            shift.ME(H11r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H11r.Mv1Squared(H21r);
        }
        else {
            shift.PE(H11r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H11r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H11r.Mv1Squared(H22r);
        }
        else {
            shift.PE(H11r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H11r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H12r.Mv1Squared(O2r);
        }
        else {
            shift.PE(H12r);
            r2 = O2r.Mv1Squared(shift);
            shift.ME(H12r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H12r.Mv1Squared(H21r);
        }
        else {
            shift.PE(H12r);
            r2 = H21r.Mv1Squared(shift);
            shift.ME(H12r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

        if (zeroShift) {
            r2 = H12r.Mv1Squared(H22r);
        }
        else {
            shift.PE(H12r);
            r2 = H22r.Mv1Squared(shift);
            shift.ME(H12r);
        }
		if(r2<1.6) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		return sum;																					        
	}//end of energy
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
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
    
    private static final long serialVersionUID = 1L;
	public double sigma , sigma2;
	public double epsilon, epsilon4;
	private etomica.space.Boundary boundary;
	private double chargeH = Electron.UNIT.toSim(0.41);
	private double chargeO = Electron.UNIT.toSim(-0.82);
	private double chargeOO, chargeOH, chargeHH;
	private IVector work, shift;
}
