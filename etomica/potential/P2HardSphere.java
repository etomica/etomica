package etomica.potential;

import etomica.Atom;
import etomica.Debug;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.space.CoordinatePairKinetic;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if spheres overlap, and is zero otherwise.  Collision diameter describes
 * size of spheres.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSphere extends Potential2HardSpherical {
    
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
   protected final Vector dr;
   protected final Tensor lastCollisionVirialTensor;
    
    public P2HardSphere() {
        this(Simulation.getDefault().space, Default.ATOM_SIZE);
    }

    public P2HardSphere(double d) {
        this(Simulation.getDefault().space, d);
    }
    public P2HardSphere(Space space) {
        this(space, Default.ATOM_SIZE);
    }
    public P2HardSphere(Space space, double d) {
        super(space);
        setCollisionDiameter(d);
        lastCollisionVirialTensor = space.makeTensor();
        dr = space.makeVector();
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
    public double collisionTime(Atom[] pair, double falseTime) {
//        dr.E(cPairNbr.reset(pair[0].coord,pair[1].coord));
        Vector dv = ((CoordinatePairKinetic)cPairNbr).resetV(pair[0].coord,pair[1].coord);
        dr.Ev1Pa1Tv2(cPairNbr.reset(),falseTime,dv);
//        dr.E(cPairNbr.dr());
//        Vector dv = ((CoordinatePairKinetic)cPairNbr).dv();
//        dr.PEa1Tv1(falseTime,dv);
        double bij = dr.dot(dv);
        double time = Double.POSITIVE_INFINITY;

        if(bij < 0.0) {
        	if (Default.FIX_OVERLAP && dr.squared() < sig2) return falseTime;
            double v2 = dv.squared();
            double discriminant = bij*bij - v2 * ( dr.squared() - sig2 );
            if(discriminant > 0) {
                time = (-bij - Math.sqrt(discriminant))/v2;
            }
        }
        if (Debug.ON && Debug.DEBUG_NOW && (Debug.allAtoms(pair) || time < 0.0)) {
        	System.out.println("atoms "+pair[0]+" and "+pair[1]+" r2 "+dr.squared()+" bij "+bij+" time "+time);
        	if (time < 0.0) throw new RuntimeException("negative collision time for hard spheres");
        }
        return time + falseTime;
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(Atom[] pair, double falseTime) {
        cPair.reset(pair[0].coord,pair[1].coord);
        ((CoordinatePairKinetic)cPair).resetV();
        dr.E(cPair.dr());
        Vector dv = ((CoordinatePairKinetic)cPair).dv();
        dr.PEa1Tv1(falseTime,dv);
        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double reducedMass = 2.0/(pair[0].type.rm() + pair[1].type.rm());
        lastCollisionVirial = reducedMass*bij;
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        //dv is now change in velocity due to collision
        dv.Ea1Tv1(lastCollisionVirialr2/reducedMass,dr);
        ((ICoordinateKinetic)pair[0].coord).velocity().PE(dv);
        ((ICoordinateKinetic)pair[1].coord).velocity().ME(dv);
        pair[0].coord.position().PEa1Tv1(-falseTime,dv);
        pair[1].coord.position().PEa1Tv1(falseTime,dv);
    }
    
    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    
    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
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
        collisionDiameter = c;
        sig2 = c*c;
    }
    public etomica.units.Dimension getCollisionDiameterDimension() {
        return etomica.units.Dimension.LENGTH;
    }
    
    /**
     * Interaction energy of the pair.
     * Zero if separation is greater than collision diameter, infinity otherwise
     */
    public double u(double r2) {
        return (r2 < sig2) ? Double.MAX_VALUE : 0.0;
    }
    
    public double energyChange() {return 0.0;}
    
}//end of P2HardSphere