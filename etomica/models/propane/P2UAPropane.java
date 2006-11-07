
package etomica.models.propane;

import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.Potential2;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Kelvin;

/** 
 * 
 * @author kofke
 *
 * SPC/E potential for water.  Requires the molecule node be an
 * AtomTreeNodeWater. 
 */
public class P2UAPropane extends Potential2 {

    public P2UAPropane(Space space) {
	    super(space);
	    setSigmaCH3(3.750);
	    setSigmaCH2(3.950);
	    setEpsilonCH3(Kelvin.UNIT.toSim(98.00));
	    setEpsilonCH2(Kelvin.UNIT.toSim(46.00));
//	    setCharges();
    }   

    public double energy(AtomSet atoms){
        double sum = 0.0;
        double r2 = 0.0;

        AtomPair pair = (AtomPair)atoms;
        AtomTreeNodeUAPropane node1 = (AtomTreeNodeUAPropane)pair.atom0.node;
        AtomTreeNodeUAPropane node2 = (AtomTreeNodeUAPropane)pair.atom1.node;

        Vector propane1UA1r = node1.UA1.coord.position();
        Vector propane2UA1r = node2.UA1.coord.position();
        Vector propane1UA2r = node1.UA2.coord.position();
        Vector propane2UA2r = node2.UA2.coord.position();
        Vector propane1UA3r = node1.UA3.coord.position();
        Vector propane2UA3r = node2.UA3.coord.position();

        
        final double core = 0.1;

        //compute CH3-CH3 distance between first UA atoms on each chain to consider truncation   
        r2 = propane1UA1r.Mv1Squared(propane2UA1r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        double s2 = sigmaCH3SQ/(r2);
        double s6 = s2*s2*s2;
        sum += epsilonCH3X4*s6*(6 - 1.0);

        r2 = propane1UA1r.Mv1Squared(propane2UA2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
// Use Lorentz-Berthelot combining rules for unlike CH3-CH2 interaction
        double sigmaIJ = (sigmaCH3+sigmaCH2)/2;
        double epsilonIJ = Math.sqrt(epsilonCH3*epsilonCH2);
        double epsilonIJ4 = 4*epsilonIJ;
        s2 = sigmaIJ/(r2);
        s6 = s2*s2*s2;
        sum += epsilonIJ4*s6*(s6 - 1.0);

        r2 = propane1UA1r.Mv1Squared(propane2UA3r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        s2 = sigmaCH3/(r2);
        s6 = s2*s2*s2;
        sum += epsilonCH3X4*s6*(s6 - 1.0);

        r2 = propane1UA2r.Mv1Squared(propane2UA1r);
        if(r2<core) return Double.POSITIVE_INFINITY;
//      Use Lorentz-Berthelot combining rules for unlike CH3-CH2 interaction
        s2 = sigmaIJ/(r2);
        s6 = s2*s2*s2;
        sum += epsilonIJ4*s6*(s6 - 1.0);
        
        r2 = propane1UA2r.Mv1Squared(propane2UA2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        s2 = sigmaCH2/(r2);
        s6 = s2*s2*s2;
        sum += epsilonCH2X4*s6*(s6 - 1.0);

        r2 = propane1UA2r.Mv1Squared(propane2UA3r);
        if(r2<core) return Double.POSITIVE_INFINITY;
//      Use Lorentz-Berthelot combining rules for unlike CH3-CH2 interaction
        s2 = sigmaIJ/(r2);
        s6 = s2*s2*s2;
        sum += epsilonIJ4*s6*(s6 - 1.0);

        r2 = propane1UA3r.Mv1Squared(propane2UA1r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        s2 = sigmaCH3/(r2);
        s6 = s2*s2*s2;
        sum += epsilonCH3X4*s6*(s6 - 1.0);

        r2 = propane1UA3r.Mv1Squared(propane2UA2r);
        if(r2<core) return Double.POSITIVE_INFINITY;
//      Use Lorentz-Berthelot combining rules for unlike CH3-CH2 interaction
        s2 = sigmaIJ/(r2);
        s6 = s2*s2*s2;
        sum += epsilonIJ4*s6*(s6 - 1.0);
        
        r2 = propane1UA3r.Mv1Squared(propane2UA3r);
        if(r2<core) return Double.POSITIVE_INFINITY;
        s2 = sigmaCH3/(r2);
        s6 = s2*s2*s2;
        sum += epsilonCH3X4*s6*(s6 - 1.0);

/*        if (sum > 10000000) {
        		System.out.println("sum (energy) = " + sum);
        }

        if (sum < 0) {
    		System.out.println("sum (energy) = " + sum);
    }
*/
        
        return sum;																					        
    }//end of energy

    public double getSigmaCH3() {return sigmaCH3;}

    private final void setSigmaCH3(double s) {
        sigmaCH3 = s;
        sigmaCH3SQ = s*s;
    }

    public double getSigmaCH2() {return sigmaCH2;}

    private final void setSigmaCH2(double s) {
        sigmaCH2 = s;
        sigmaCH2SQ = s*s;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public double getEpsilonCH3() {return epsilonCH3;}
    
    private final void setEpsilonCH3(double eps) {
        epsilonCH3 = eps;
        epsilonCH3X4 = 4*epsilonCH3;
    }
    
    public double getEpsilonCH2() {return epsilonCH2;}
    
    private final void setEpsilonCH2(double eps) {
        epsilonCH2 = eps;
        epsilonCH2X4 = 4*epsilonCH2;
    }

    public void setPhase(Phase phase) {
        nearestImageTransformer = phase.getBoundary();
    }

    public double sigmaCH3, sigmaCH3SQ, sigmaCH2, sigmaCH2SQ;
    public double epsilonCH3, epsilonCH2, epsilonCH3X4, epsilonCH2X4;
	protected NearestImageTransformer nearestImageTransformer;

}