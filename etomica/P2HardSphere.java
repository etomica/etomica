package etomica;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if disks overlap, and is zero otherwise.  Collision diameter describes
 * size of disks.
 * Suitable for use in space of any dimension.
 *
 * @author David Kofke
 */
public class P2HardSphere extends Potential2 implements Potential2Hard, EtomicaElement {
    
    public String getVersion() {return "P2HardSphere:01.07.03/"+Potential2.VERSION;}

   /**
    * Separation at which disks first overlap
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
        this(Simulation.instance, Default.ATOM_SIZE);
    }

    public P2HardSphere(double d) {
        this(Simulation.instance, d);
    }
    public P2HardSphere(Simulation sim) {
        this(sim, Default.ATOM_SIZE);
    }
    public P2HardSphere(Simulation sim, double d) {
        super(sim);
        setCollisionDiameter(d);
        lastCollisionVirialTensor = sim.space().makeTensor();
        dr = sim.space().makeVector();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple hard-sphere potential");
        return info;
    }

    /**
     * Time to collision of pair, assuming free-flight kinematics
     */
    public double collisionTime(AtomPair pair) {
        double r2 = pair.r2();
        double bij = pair.vDotr();
        if(r2 < sig2) {return (bij > 0) ? Double.MAX_VALUE : 0.0;}  //inside wall; no collision
        double time = Double.MAX_VALUE;

        if(bij < 0.0) {
          double velocitySquared = pair.v2();
          double discriminant = bij*bij - velocitySquared * ( r2 - sig2 );
          if(discriminant > 0) {
            time = (-bij - Math.sqrt(discriminant))/velocitySquared;
          }
        }
        return time;
    }
    
    /**
     * Implements collision dynamics and updates lastCollisionVirial
     */
    public void bump(AtomPair pair) {
        double r2 = pair.r2();
        dr.E(pair.dr());  //used by lastCollisionVirialTensor
        lastCollisionVirial = 2.0/(pair.atom1().coord.rm() + pair.atom2().coord.rm())*pair.vDotr();
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        pair.cPair.push(lastCollisionVirialr2);
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
    public double energy(AtomPair pair) {
        return (pair.r2() < sig2) ? Double.MAX_VALUE : 0.0;
    }
    
    
}//end of P2HardSphere