
package etomica.models.water;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.exception.MethodNotImplementedException;
import etomica.phase.Phase;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.space.CoordinatePair;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * 
 * @author kofke
 *
 * TIP4P potential for water.  Requires the molecule node to be an
 * AtomTreeNodeWaterTIP4P.
 */
public class P2WaterTIP4P extends Potential2 implements Potential2Soft {

    public P2WaterTIP4P(Space space) {
	    super(space);
	    setSigma(3.1540);
	    setEpsilon(Kelvin.UNIT.toSim(78.02));
	    setCharges();
    }   

    public double energy(AtomSet atoms){
        double sum = 0.0;
        double r2 = 0.0;

        AtomPair pair = (AtomPair)atoms;
        AtomTreeNodeWater4P node1 = (AtomTreeNodeWater4P)pair.atom0.node;
        AtomTreeNodeWater4P node2 = (AtomTreeNodeWater4P)pair.atom1.node;

        Vector O1r = node1.O.coord.position();
        Vector O2r = node2.O.coord.position();
        Vector H11r = node1.H1.coord.position();
        Vector H12r = node1.H2.coord.position();
        Vector H21r = node2.H1.coord.position();
        Vector H22r = node2.H2.coord.position();

        Vector M1r = node1.M.coord.position();
        Vector M2r = node2.M.coord.position();
        
		
        final double core = 0.1;

        //compute O-O distance to consider truncation   
        r2 = O1r.Mv1Squared(O2r);
        
        //System.out.println("O-O distance is " + r2);
        
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOO/Math.sqrt(r2);
        double s2 = sigma2/(r2);
        double s6 = s2*s2*s2;
        sum += epsilon4*s6*(s6 - 1.0);

        // block comment; for TIP4P water O charge is zero; no contribution to sum
/*        r2 = O1r.Mv1Squared(H21r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOH/Math.sqrt(r2);

        r2 = O1r.Mv1Squared(H22r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOH/Math.sqrt(r2);

        r2 = H11r.Mv1Squared(O2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOH/Math.sqrt(r2);

        r2 = H12r.Mv1Squared(O2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOH/Math.sqrt(r2);

        r2 = M1r.Mv1Squared(O2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOM/Math.sqrt(r2);

        System.out.println("sum is " + sum);

        r2 = M2r.Mv1Squared(O1r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeOM/Math.sqrt(r2);

        System.out.println("sum is " + sum);
*/
        r2 = H11r.Mv1Squared(H21r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHH/Math.sqrt(r2);

        r2 = H11r.Mv1Squared(H22r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHH/Math.sqrt(r2);

        r2 = H12r.Mv1Squared(H21r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHH/Math.sqrt(r2);

        r2 = H12r.Mv1Squared(H22r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHH/Math.sqrt(r2);

//        System.out.println("sum of all O-H terms is " + sum);
        
        r2 = M1r.Mv1Squared(H21r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHM/Math.sqrt(r2);
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H22r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHM/Math.sqrt(r2);
        //System.out.println("sum is " + sum);

        r2 = M2r.Mv1Squared(H11r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHM/Math.sqrt(r2);

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H12r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeHM/Math.sqrt(r2);
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(M2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        sum += chargeMM/Math.sqrt(r2);
        //System.out.println("sum is " + sum);

        //System.out.println("energy is " + sum);
        return sum;																					        
    }//end of energy

    public Vector gradient(AtomSet pair){
        throw new MethodNotImplementedException();
    }
    public double hyperVirial(AtomSet pair){
        throw new MethodNotImplementedException();
    }
    public double integral(double rC){
        throw new MethodNotImplementedException();
    }
    public double virial(AtomSet pair){
        throw new MethodNotImplementedException();
    }

    public double getSigma() {return sigma;}

    private final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
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
        chargeMM = chargeM * chargeM;
        chargeOM = chargeO * chargeM;
        chargeHM = chargeH * chargeM;
    }

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.getBoundary());
    }

    public double sigma , sigma2;
    public double epsilon, epsilon4;
    private double chargeH = Electron.UNIT.toSim(0.52);
    private double chargeO = Electron.UNIT.toSim(0);
    private double chargeM = Electron.UNIT.toSim(-1.04);
    private double chargeOO, chargeOH, chargeHH, chargeOM, chargeHM, chargeMM;
    protected CoordinatePair cPair;
}