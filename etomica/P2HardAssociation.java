package etomica;

/**
 * Purely attractive square-well potential with no repulsive core.
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
public class P2HardAssociation extends Potential2 implements PotentialHard {

    private double wellDiameter, wellDiameterSquared;
    private double epsilon;
    private double r2;
    private double lastCollisionVirial, lastCollisionVirialr2;
    private final Space.Tensor lastCollisionVirialTensor;
    private final Space.Vector dr;
    
    public P2HardAssociation() {
        this(Simulation.getDefault().space, Default.POTENTIAL_CUTOFF_FACTOR, Default.POTENTIAL_WELL);
    }
    public P2HardAssociation(double wellDiameter, double epsilon) {
        this(Simulation.getDefault().space, wellDiameter, epsilon);
    }
    public P2HardAssociation(Space space, double wellDiameter, double epsilon) {
        super(space);
        setEpsilon(epsilon);
        setWellDiameter(wellDiameter);
        dr = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
    }
    
   /**
    * Implements the collision dynamics.  Does not deal with the hard cores, only the wells.  This
    * section is essentially the same as PotentialSquareWell without the hard core section.
    */
    public void bump(Atom[] pair) {
        double eps = 1e-6;
        cPair.reset(pair[0].coord,pair[1].coord);
        r2 = cPair.r2();
        double bij = cPair.vDotr();
        //ke is kinetic energy from the components of velocity
        double reduced_m = 1/(pair[0].coord.rm() + pair[1].coord.rm());
        dr.E(cPair.dr());
        double ke = bij*bij*reduced_m/(2*r2);
        double r2New;
        if (bij > 0.0) {    //Separating
            if (ke < epsilon) {    //Not enough energy to escape the well
                lastCollisionVirial = 2*bij*reduced_m;
                r2New = (1 - eps)*wellDiameterSquared;
            }
            else {  //Escaping the well
                lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
                r2New = (1 + eps)*wellDiameterSquared;
            }
        }
        else {          //Approaching
            lastCollisionVirial = reduced_m*(bij + Math.sqrt(bij*bij + 2.0*r2*epsilon/reduced_m));  //might need to double check
            r2New = (1 - eps)*wellDiameterSquared;
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        cPair.push(lastCollisionVirialr2);
        if(r2New != r2) cPair.setSeparation(r2New);
    }       //end of bump
    
    
    public double lastCollisionVirial() {return lastCollisionVirial;}
    
    public Space.Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }
    
   /**
    * Computes the next time of collision of the given atomPair assuming free flight.  Only computes the next
    * collision of the wells.  Takes into account both separation and convergence.
    */
    public double collisionTime(Atom[] pair) {
        double discr = 0.0;
        cPair.reset(pair[0].coord,pair[1].coord);
		cPair.resetV();
       double bij = cPair.vDotr();
        double r2 = cPair.r2();
        double v2 = cPair.v2();
        double tij = Double.MAX_VALUE;
        
        if (r2 < wellDiameterSquared) {         //check to see if already inside wells
            discr = bij*bij - v2 * (r2 - wellDiameterSquared );
            tij = (-bij + Math.sqrt(discr))/v2;
        }
        else {
            if (bij < 0.) { //Approaching
                discr = bij*bij - v2 * (r2 - wellDiameterSquared );
                if (discr > 0.) {
                    tij = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        return tij;
    }
    
  /**
   * Returns -epsilon if less than well diameter, or zero otherwise.
   */
    public double energy(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
        return (cPair.r2() < wellDiameterSquared) ?  -epsilon : 0.0;
    }
    
   /**
    * Accessor for well-diameter.
    * Since the well-diameter is not a multiplier in this potential as in square well, it is necessary
    * to be able to set this manually if so desired.
    */
    public double getWellDiameter() {return wellDiameter;}
    
   /**
    * Accessor for well-diameter.
    * Allows manual changing of the well diameter since it is not merely a multiple of the core-diameter
    * as in square well.
    */
    
    public void setWellDiameter(double w) {
        wellDiameter = w;
        wellDiameterSquared = w*w;
    }
    
   /**
    * Accessor method for depth of well.
    */
    public double getEpsilon() {return epsilon;}
    
   /**
    * Accessor method for depth of well.
    */
    public void setEpsilon(double s) {
        epsilon = s;
    }
}//end of P2HardAssociation