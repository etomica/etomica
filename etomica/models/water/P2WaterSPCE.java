
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
 * SPC/E potential for water.  Requires the molecule node be an
 * AtomTreeNodeWater. 
 */
public class P2WaterSPCE extends Potential2 implements Potential2Soft {

    public P2WaterSPCE(Space space) {
	    super(space);
	    setSigma(3.1670);
	    setEpsilon(Kelvin.UNIT.toSim(78.21));
	    setCharges();
    }   

    public double energy(AtomSet atoms){
        double sum = 0.0;
        double r2 = 0.0;

        AtomPair pair = (AtomPair)atoms;
        AtomTreeNodeWater3P node1 = (AtomTreeNodeWater3P)pair.atom0.node;
        AtomTreeNodeWater3P node2 = (AtomTreeNodeWater3P)pair.atom1.node;

        Vector O1r = node1.O.coord.position();
        Vector O2r = node2.O.coord.position();
        Vector H11r = node1.H1.coord.position();
        Vector H12r = node1.H2.coord.position();
        Vector H21r = node2.H1.coord.position();
        Vector H22r = node2.H2.coord.position();

        final double core = 0.1;

        //compute O-O distance to consider truncation   
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
    }

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.boundary());
    }

    public double sigma , sigma2;
    public double epsilon, epsilon4;
    private double chargeH = Electron.UNIT.toSim(0.4238);
    private double chargeO = Electron.UNIT.toSim(-0.8476);
    private double chargeOO, chargeOH, chargeHH;
    protected CoordinatePair cPair;
}