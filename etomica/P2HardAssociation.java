package etomica;

/**
 * Purely attractive square-well potential with no repulsive core.
 *
 * @author Rob Riggleman
 * @author David Kofke
 */
public class P2HardAssociation extends Potential2Hard implements EtomicaElement {

    public String getVersion() {return "P2HardAssociation:01.07.03/"+Potential2.VERSION;}

    private double wellDiameter, wellDiameterSquared;
    private double epsilon;
    private double r2;
    private double lastCollisionVirial, lastCollisionVirialr2;
    private final Space.Tensor lastCollisionVirialTensor;
    private final Space.Vector dr;
    
    public P2HardAssociation() {
        this(Simulation.instance.hamiltonian.potential, Default.POTENTIAL_CUTOFF_FACTOR, Default.POTENTIAL_WELL);
    }
    public P2HardAssociation(double wellDiameter, double epsilon) {
        this(Simulation.instance.hamiltonian.potential, wellDiameter, epsilon);
    }
    public P2HardAssociation(PotentialGroup parent, double wellDiameter, double epsilon) {
        super(parent);
        setEpsilon(epsilon);
        setWellDiameter(wellDiameter);
        dr = parentSimulation().space().makeVector();
        lastCollisionVirialTensor = parentSimulation().space().makeTensor();
    }
    
   /**
    * Implements the collision dynamics.  Does not deal with the hard cores, only the wells.  This
    * section is essentially the same as PotentialSquareWell without the hard core section.
    */
    public void bump(AtomPair pair) {
        double eps = 1e-6;
        r2 = pair.r2();
        double bij = pair.vDotr();
        //ke is kinetic energy from the components of velocity
        double reduced_m = 1/(pair.atom1().coord.rm() + pair.atom2().coord.rm());
        dr.E(pair.dr());
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
        pair.cPair.push(lastCollisionVirialr2);
        if(r2New != r2) pair.cPair.setSeparation(r2New);
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
    public double collisionTime(AtomPair pair) {
        double discr = 0.0;
        double bij = pair.vDotr();
        double r2 = pair.r2();
        double v2 = pair.v2();
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
    public double energy(AtomPair pair) {
        return (pair.r2() < wellDiameterSquared) ?  -epsilon : 0.0;
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