package etomica;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if spheres overlap, and is zero otherwise.  Collision diameter describes
 * size of spheres.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSphere extends Potential2 implements PotentialHard {
    
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
   protected final Space.Vector dr;
   protected final Space.Tensor lastCollisionVirialTensor;
    
    public P2HardSphere() {
        this(Simulation.instance.hamiltonian.potential, Default.ATOM_SIZE);
    }

    public P2HardSphere(double d) {
        this(Simulation.instance.hamiltonian.potential, d);
    }
    public P2HardSphere(SimulationElement parent) {
        this(parent, Default.ATOM_SIZE);
    }
    public P2HardSphere(SimulationElement parent, double d) {
        super(parent);
        setCollisionDiameter(d);
        lastCollisionVirialTensor = simulation().space().makeTensor();
        dr = simulation().space().makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple hard-sphere potential");
        return info;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        cPair.resetV();
        double bij = cPair.vDotr();
        if(r2 < sig2) {return (bij > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
        double time = Double.MAX_VALUE;

        if(bij < 0.0) {
          double velocitySquared = cPair.v2();
          double discriminant = bij*bij - velocitySquared * ( r2 - sig2 );
          if(discriminant > 0) {
            time = (-bij - Math.sqrt(discriminant))/velocitySquared;
          }
        }
    	System.out.println("in collisionTime "+pair[0]+" "+pair[1]+" "+time);
        return time;
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        dr.E(cPair.dr());  //used by lastCollisionVirialTensor
		cPair.resetV();
        lastCollisionVirial = 2.0/(pair[0].coord.rm() + pair[1].coord.rm())*cPair.vDotr();
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        cPair.push(lastCollisionVirialr2);
    }
    
    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }
    
    public Space.Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
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
    public double energy(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
    	if(cPair.r2() < sig2) {
//    		System.out.println("uh oh in p2hardsphere");
    	}
        return (cPair.r2() < sig2) ? Double.MAX_VALUE : 0.0;
    }
    
    
}//end of P2HardSphere