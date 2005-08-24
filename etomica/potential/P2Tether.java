package etomica.potential;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.space.CoordinatePairKinetic;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Dimension;
/**
 * Potential that acts like a hard string connecting the centers of two atoms.
 * Meant for use as an intra-molecular interaction.
 * Interaction of atoms is zero if separated by less than the tether length.  Atoms
 * undergo an impulsive attractive collision when attempting to separate by more than the tether distance.
 *
 * @author David Kofke
 */
public class P2Tether extends Potential2HardSpherical {

    public P2Tether(Space space) {
        super(space, new CoordinatePairKinetic(space));
        setTetherLength(0.75*Default.ATOM_SIZE);
        lastCollisionVirialTensor = space.makeTensor();
        dr = space.makeVector();
        cPair = (CoordinatePairKinetic)coordinatePair;
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard string between adjacent atoms, hard sphere for non-adjacents");
        return info;
    }

    /**
     * Accessor method for the tether distance
     */
    public double getTetherLength() {return tetherLength;}
    /**
     * Accessor method for the tether distance
     */
    public void setTetherLength(double t) {
        tetherLength = t;
        tetherLengthSquared = t*t;
    }
    public Dimension getTetherLengthDimension() {return Dimension.LENGTH;}

    /**
     * Implements collision dynamics for pair attempting to separate beyond tether distance
     */
    public final void bump(AtomSet atoms, double falseTime) {
        AtomPair pair = (AtomPair)atoms;
        cPair.reset(pair);
        cPair.resetV();
        dr.E(cPair.dr());
        Vector dv = cPair.dv();
        dr.PEa1Tv1(falseTime,dv);
        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = ((AtomTypeLeaf)pair.atom0.type).rm();
        double rm1 = ((AtomTypeLeaf)pair.atom1.type).rm();
        lastCollisionVirial = 2.0/(rm0 + rm1)*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        ((ICoordinateKinetic)pair.atom0.coord).velocity().PEa1Tv1( rm0,dv);
        ((ICoordinateKinetic)pair.atom1.coord).velocity().PEa1Tv1(-rm1,dv);
        pair.atom0.coord.position().PEa1Tv1(-falseTime*rm0,dv);
        pair.atom1.coord.position().PEa1Tv1( falseTime*rm1,dv);
    }


    public final double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    public final Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }

    /**
     * Returns the tether length. The potential range is in fact infinite, but if
     * the integrator is generating configurations correctly, there will be no atoms
     * interacting beyond the tetherlength distance.  
     */
    public double getRange() {
        return tetherLength;
    }
  
    /**
     * Time at which two atoms will reach the end of their tether, assuming free-flight kinematics
     */
    public final double collisionTime(AtomSet pair, double falseTime) {
        cPair.reset((AtomPair)pair);
        cPair.resetV();
        dr.E(cPair.dr());
        Vector dv = cPair.dv();
        dr.Ea1Tv1(falseTime,dv);
        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        if(Default.FIX_OVERLAP && r2 > tetherLengthSquared && bij > 0) {return 0.0;}  //outside tether, moving apart; collide now
        double discr = bij*bij - v2 * ( r2 - tetherLengthSquared );
        return (-bij + Math.sqrt(discr))/v2 + falseTime;
    }
  
    /**
     * Returns infinity if separation is greater than tether distance, zero otherwise
     */
    public double u(double r2) {
        return (r2 > tetherLengthSquared) ? Double.POSITIVE_INFINITY : 0.0;
    }
    
    public double energyChange() {return 0.0;}
    
    private double tetherLength, tetherLengthSquared;
    private double lastCollisionVirial = 0.0;
    private double lastCollisionVirialr2 = 0.0;
    private final Vector dr;
    private final Tensor lastCollisionVirialTensor;
    private final CoordinatePairKinetic cPair;
   
}//end of P2Tether
  