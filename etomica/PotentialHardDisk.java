package etomica;

/**
 * Basic hard-(rod/disk/sphere) potential.
 * Energy is infinite if disks overlap, and is zero otherwise.  Collision diameter describes
 * size of disks.
 * Suitable for use in space of any dimension.
 */
public class PotentialHardDisk extends Potential2HardAbstract implements EtomicaElement
{
    public String getVersion() {return "PotentialHardDisk:01.06.17/"+Potential2HardAbstract.VERSION;}

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
    
    public PotentialHardDisk() {
        this(Simulation.instance, Default.ATOM_SIZE);
    }

    public PotentialHardDisk(double d) {
        this(Simulation.instance, d);
    }
    public PotentialHardDisk(Simulation sim) {
        this(sim, Default.ATOM_SIZE);
    }
    public PotentialHardDisk(Simulation sim, double d) {
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
        lastCollisionVirial = 2.0/(pair.atom1().rm() + pair.atom2().rm())*pair.vDotr();
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
    
    /**
     * Correction for trucation of the potential
     * Identically zero for this model
     */
    public double energyLRC(int n1, int n2, double V) {return 0.0;}
    
    /**
     * Returns true if separation of pair is less than collision diameter, false otherwise
     */
    public boolean overlap(AtomPair pair) {return pair.r2() < sig2;}
    
}