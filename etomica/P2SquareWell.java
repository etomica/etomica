package etomica;
import etomica.units.Dimension;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 */
public class P2SquareWell extends Potential2 implements PotentialHard {

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected Space.Tensor lastCollisionVirialTensor;
    private double lastEnergyChange;
    protected Space.Vector dr;

    public P2SquareWell() {
        this(Simulation.getDefault().space);
    }
    public P2SquareWell(Space space) {
        this(space,Default.ATOM_SIZE, Default.POTENTIAL_CUTOFF_FACTOR, Default.POTENTIAL_WELL);
    }
    public P2SquareWell(double coreDiameter, double lambda, double epsilon) {
        this(Simulation.getDefault().space, coreDiameter, lambda, epsilon);
    }

    public P2SquareWell(Space space, double coreDiameter, double lambda, double epsilon) {
        super(space);
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        dr = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Simple hard repulsive core with a square-well region of attraction");
        return info;
    }

    /**
     * Returns true if separation of pair is less than core diameter, false otherwise
     */
    public boolean overlap(Atom[] pair) {
        cPair.reset(pair[0].coord,pair[1].coord);
        return cPair.r2() < coreDiameterSquared;
    }

    /**
     * Implements collision dynamics between two square-well atoms.
     * Includes all possibilities involving collision of hard cores, and collision of wells
     * both approaching and diverging
     */
    public void bump(Atom[] pair) {
        cPair.reset(pair[0].coord,pair[1].coord);
        double eps = 1.0e-10;
        double r2 = cPair.r2();
        cPair.resetV();
        double bij = cPair.vDotr();
        double reduced_m = 1.0/(pair[0].coord.rm() + pair[1].coord.rm());
        double r2New;
        if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
            if (Debug.ON && Math.abs(r2 - coreDiameterSquared)/coreDiameterSquared > 1.e-9) {
                throw new RuntimeException("atoms "+pair[0]+" "+pair[1]+" not at the right distance "+r2+" "+coreDiameterSquared);
            }
            lastCollisionVirial = 2.0*reduced_m*bij;
            r2New = r2;
            lastEnergyChange = 0.0;
        }
        else {    // Well collision
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(bij > 0.0) {         // Separating
                if(ke < epsilon) {     // Not enough kinetic energy to escape
                    if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                        throw new RuntimeException("atoms "+pair[0]+" "+pair[1]+" not at the right distance "+r2+" "+wellDiameterSquared);
                    }
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    r2New = (1-eps)*wellDiameterSquared;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                        throw new RuntimeException("atoms "+pair[0]+" "+pair[1]+" not at the right distance "+r2+" "+wellDiameterSquared);
                    }
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*epsilon/reduced_m));
                    r2New = (1+eps)*wellDiameterSquared;
                    lastEnergyChange = epsilon;
                }
            }
            else {                  // Approaching
                if (Debug.ON && Math.abs(r2 - wellDiameterSquared)/wellDiameterSquared > 1.e-9) {
                    throw new RuntimeException("atoms "+pair[0]+" "+pair[1]+" not at the right distance "+r2+" "+wellDiameterSquared);
                }
                lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij+2.0*r2*epsilon/reduced_m));
                r2New = (1-eps)*wellDiameterSquared;
                lastEnergyChange = -epsilon;
            }
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        cPair.push(lastCollisionVirialr2);
        if(r2New != r2) cPair.setSeparation(r2New);
    }//end of bump method

    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public Space.Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.E(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(Atom[] pair) {
        cPair.reset(pair[0].coord,pair[1].coord);
        double discr = 0.0;

        cPair.resetV();
        double bij = cPair.vDotr();
        double r2 = cPair.r2();
        double v2 = cPair.v2();

        double tij = Double.MAX_VALUE;

        if(r2 < wellDiameterSquared) {  // Already inside wells

            if(bij < 0.0) {    // Check for hard-core collision
                if(Default.FIX_OVERLAP && r2 < coreDiameterSquared) {   // Inside core; collision now
                    return 0.0;
                }

                discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                if(discr > 0) {  // Hard cores collide next
                    tij = (-bij - Math.sqrt(discr))/v2;
                }
                else {           // Moving toward each other, but wells collide next
                    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                    tij = (-bij + Math.sqrt(discr))/v2;
                }
            }
            else {           // Moving away from each other, wells collide next
                discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                tij = (-bij + Math.sqrt(discr))/v2;
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                if(discr > 0) {
                    tij = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        return tij;
    }

    /**
     * Returns infinity if overlapping, -epsilon if otherwise less than well diameter, or zero if neither.
     */
    public double energy(Atom[] pair) {
        cPair.reset(pair[0].coord,pair[1].coord);
        double r2 = cPair.r2();
        return ( r2 < wellDiameterSquared) ? 
            ((r2 < coreDiameterSquared) ? Double.MAX_VALUE : -epsilon) : 0.0;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return coreDiameter;}
    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getCoreDiameterDimension() {return Dimension.LENGTH;}

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {return lambda;}
    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when 
     * this is changed
     */
    public void setLambda(double lam) {
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getLambdaDimension() {return Dimension.NULL;}

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for depth of well
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Dimension.ENERGY;}

}
  
