
package etomica.virial.simulations;

import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.*;
import etomica.models.water.*;

/** 
 * 
 * @author kofke
 *
 * Special version of water model that permits no truncation of the potential
 * and no periodic boundary. Used for cluster calculations.
 */
public class P2WaterSPCE extends Potential2 implements Potential2.Soft {

	public P2WaterSPCE(SimulationElement parent) {
		super(parent);
		setSigma(3.1670);
		setEpsilon(Kelvin.UNIT.toSim(78.21));
		setCharges();
	}   
	public double energy(AtomPair pair){
		double sum = 0.0;
		double r2 = 0.0;
			
		AtomTreeNodeWater node1 = (AtomTreeNodeWater)pair.atom1().node;
		AtomTreeNodeWater node2 = (AtomTreeNodeWater)pair.atom2().node;
		
		//compute O-O distance to consider truncation	
		Space3D.Vector O1r = (Space3D.Vector)node1.O.coord.position();
		Space3D.Vector O2r = (Space3D.Vector)node2.O.coord.position();
		Space3D.Vector H11r = (Space3D.Vector)node1.H1.coord.position();
		Space3D.Vector H12r = (Space3D.Vector)node1.H2.coord.position();
		Space3D.Vector H21r = (Space3D.Vector)node2.H1.coord.position();
		Space3D.Vector H22r = (Space3D.Vector)node2.H2.coord.position();
        
        final double core = 0.1;
        				
		r2 = O1r.Mv1Squared(O2r);
		if(r2<core) return Double.POSITIVE_INFINITY;	
		sum += chargeOO/Math.sqrt(r2);
		double s2 = sigma2/(r2);
		double s6 = s2*s2*s2;
		sum += epsilon4*s6*(s6 - 1.0);
		
		r2 = O1r.Mv1Squared(H21r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);
		
		r2 = O1r.Mv1Squared(H22r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = H11r.Mv1Squared(O2r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = H11r.Mv1Squared(H21r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = H11r.Mv1Squared(H22r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = H12r.Mv1Squared(O2r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeOH/Math.sqrt(r2);

		r2 = H12r.Mv1Squared(H21r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

		r2 = H12r.Mv1Squared(H22r);
		if(r2<core) return Double.POSITIVE_INFINITY;
		sum += chargeHH/Math.sqrt(r2);

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
		chargeOO = chargeO * chargeO;
		chargeOH = chargeO * chargeH;
		chargeHH = chargeH * chargeH;
	}
    
	public double sigma , sigma2;
	public double epsilon, epsilon4;
	private Atom O1, H11, H12, O2, H21, H22;
	private double chargeH = Electron.UNIT.toSim(0.4238);
	private double chargeO = Electron.UNIT.toSim(-0.8476);
	private double chargeOO, chargeOH, chargeHH;
}