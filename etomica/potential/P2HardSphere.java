package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.simulation.Simulation;
import etomica.space.ICoordinateKinetic;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.util.Debug;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if spheres overlap, and is zero otherwise.  Collision diameter describes
 * size of spheres.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSphere extends Potential2HardSpherical {
    
    private static final long serialVersionUID = 1L;
    /**
     * Separation at which spheres first overlap
     */
    protected double collisionDiameter;
   
    /**
     * Square of collisionDiameter
     */
    protected double sig2;
    protected double lastCollisionVirial = 0.0;
    protected double lastCollisionVirialr2 = 0.0;
    protected final boolean ignoreOverlap;
    protected final IVector dv;
    protected final Tensor lastCollisionVirialTensor;
    
    public P2HardSphere(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().ignoreOverlap);
    }
    public P2HardSphere(Space space, double d, boolean ignoreOverlap) {
        super(space);
        setCollisionDiameter(d);
        lastCollisionVirialTensor = space.makeTensor();
        dv = space.makeVector();
        this.ignoreOverlap = ignoreOverlap;
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple hard-sphere potential");
        return info;
    }
    
    public double getRange() {
    	return collisionDiameter;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(AtomSet pair, double falseTime) {
        ICoordinateKinetic coord0 = (ICoordinateKinetic)((AtomPair)pair).atom0;
        ICoordinateKinetic coord1 = (ICoordinateKinetic)((AtomPair)pair).atom1;
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        nearestImageTransformer.nearestImage(dr);

        double bij = dr.dot(dv);
        double time = Double.POSITIVE_INFINITY;

        if(bij < 0.0) {
        	if (ignoreOverlap && dr.squared() < sig2) return falseTime;
            double v2 = dv.squared();
            double discriminant = bij*bij - v2 * ( dr.squared() - sig2 );
            if(discriminant > 0) {
                time = (-bij - Math.sqrt(discriminant))/v2;
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || time < 0.0)) {
        	System.out.println("atoms "+pair+" r2 "+dr.squared()+" bij "+bij+" time "+time);
        	if (time < 0.0) throw new RuntimeException("negative collision time for hard spheres");
        }
        return time + falseTime;
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(AtomSet pair, double falseTime) {
        Atom atom0 = ((AtomPair)pair).atom0;
        Atom atom1 = ((AtomPair)pair).atom1;
        ICoordinateKinetic coord0 = (ICoordinateKinetic)atom0;
        ICoordinateKinetic coord1 = (ICoordinateKinetic)atom1;
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        nearestImageTransformer.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double rm0 = ((AtomTypeLeaf)atom0.getType()).rm();
        double rm1 = ((AtomTypeLeaf)atom1.getType()).rm();
        double reducedMass = 2.0/(rm0 + rm1);
        lastCollisionVirial = reducedMass*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        //dv is now change in velocity due to collision
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        coord0.getVelocity().PEa1Tv1( rm0,dv);
        coord1.getVelocity().PEa1Tv1(-rm1,dv);
        coord0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        coord1.getPosition().PEa1Tv1( falseTime*rm1,dv);
    }
    
    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    
    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
//        lastCollisionVirialTensor.DE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;        
    }
    
    /**
     * Accessor method for collision diameter
     */
    public double getCollisionDiameter() {return collisionDiameter;}
    /**
     * Accessor method for collision diameter
     */
    public void setCollisionDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        collisionDiameter = c;
        sig2 = c*c;
    }
    public etomica.units.Dimension getCollisionDiameterDimension() {
        return etomica.units.Length.DIMENSION;
    }
    
    /**
     * Interaction energy of the pair.
     * Zero if separation is greater than collision diameter, infinity otherwise
     */
    public double u(double r2) {
        if (r2 < sig2) {
            return Double.POSITIVE_INFINITY;
        }
        return 0.0;
    }
    
    public double energyChange() {return 0.0;}
    
}//end of P2HardSphere