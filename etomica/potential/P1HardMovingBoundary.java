package etomica.potential;

import etomica.Atom;
import etomica.Debug;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Space;
import etomica.Integrator.IntervalEvent;
import etomica.integrator.IntegratorHard;
import etomica.space.Boundary;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Potential that places hard repulsive walls that move and 
 * accelerate subject to an external force field (pressure).
 */
 
public class P1HardMovingBoundary extends Potential1 implements PotentialHard, Integrator.IntervalListener {
    
    private double collisionRadius = 0.0;
    private final int D;
    private final int wallD;
    private double wallPosition;
    private double wallVelocity;
    private double wallMass;
    private double force;
    private double pressure;
    private double lastVirial;
    private final Boundary pistonBoundary;
    private Atom atom;
    
    /**
     * Constructor for a hard moving (and accelerating) boundary.
     * @param space
     * @param wallDimension dimension which the wall is perpendicular to
     */
    public P1HardMovingBoundary(Space space, Boundary boundary, int wallDimension, double mass) {
        super(space);
        D = space.D();
        wallD = wallDimension;
        wallPosition = 0.0;
        wallMass = mass;
        force = 0.0;
        pistonBoundary = boundary;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard moving wall");
        return info;
    }
    
    public int getPriority() {return 100;}
    
    public void setWallPosition(double p) {
        wallPosition = p;
    }
    public double getWallPosition() {
        return wallPosition;
    }
    public double getWallVelocity() {
        return wallVelocity;
    }
    public void setWallVelocity(double v) {
        wallVelocity = v;
    }
    public void setForce(double f) {
        force = f;
        pressure = -1.0;
    }
    public double getForce() {
        return force;
    }
    public void setPressure(double p) {
        pressure = p;
    }
    
    /**
     * @return Returns the mass.
     */
    public double getMass() {
        return wallMass;
    }
    /**
     * @param mass The mass to set.
     */
    public void setMass(double mass) {
        wallMass = mass;
    }
    
    public double energy(Atom[] a) {
        double dx = a[0].coord.position().x(wallD) - wallPosition;
        if (dx*dx < collisionRadius*collisionRadius) {
            return Double.POSITIVE_INFINITY;
        }
        return 0.0;
    }
     
    public double energyChange() {return 0.0;}
    
    public double collisionTime(Atom[] a, double falseTime) {
    	atom = a[0];
        double dr = atom.coord.position().x(wallD) - wallPosition;
        double dv = ((ICoordinateKinetic)atom.coord).velocity().x(wallD) - wallVelocity;
        dr += dv*falseTime;
        double t = Double.POSITIVE_INFINITY;
        if (dv * dr < 0.0) {
            if (dv > 0.0) {
                dr = dr - collisionRadius;
            }
            else {
                dr = -dr + collisionRadius;
            }
            t = dr / dv;
        }
        //approaching wall
        if (Default.FIX_OVERLAP && t<0.0) t = 0.0;
        return t + falseTime;
    }
                
    public void bump(Atom[] a, double falseTime) {
    	atom = a[0];
        double r = atom.coord.position().x(wallD);
        Vector v = ((ICoordinateKinetic)atom.coord).velocity();
        lastVirial = 2.0/(1/wallMass + atom.type.rm())*(wallPosition-r)*(wallVelocity-v.x(wallD));
        if (pressure >= 0.0) {
            double area = 1.0;
            if (pressure > 0.0) {
                final Vector dimensions = pistonBoundary.dimensions();
                for (int i=0; i<D; i++) {
                    if (i != wallD) {
                        area *= dimensions.x(i);
                    }
                }
            }
            force = pressure/area;
        }
        double trueWallVelocity = wallVelocity + falseTime*force/wallMass;
        double dv = 2.0/(1/wallMass + atom.type.rm())*(trueWallVelocity-v.x(wallD));
        lastVirial = -Math.abs(dv)*collisionRadius;
        if (Debug.ON) {
            double trueWallPosition = wallPosition + wallVelocity*falseTime + 0.5*falseTime*(force/wallMass)*(force/wallMass);
            if (Math.abs(trueWallPosition-(r+v.x(wallD)*falseTime)) - collisionRadius > 1.e-9*collisionRadius) {
                System.out.println("bork at "+falseTime+" ! "+atom+" "+r+" "+v);
                System.out.println("wall bork! "+wallPosition+" "+wallVelocity+" "+force);
                System.out.println("dr bork! "+((r+v.x(wallD)*falseTime)-trueWallPosition)+" "+collisionRadius);
                throw new RuntimeException("bork!");
            }
        }
        v.setX(wallD,v.x(wallD)+dv*atom.type.rm());
        atom.coord.position().setX(wallD,r-dv*atom.type.rm()*falseTime);
        wallVelocity -= dv/wallMass;
        wallPosition += dv/wallMass*falseTime;
    }
    
    public double lastCollisionVirial() {return lastVirial;}
    
    /**
     * not yet implemented.
     */
    public Tensor lastCollisionVirialTensor() {return null;}
    
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public void setCollisionRadius(double d) {collisionRadius = d;}
    /**
     * Distance from the center of the sphere to the boundary at collision.
     */
    public double getCollisionRadius() {return collisionRadius;}
    /**
     * Indicates collision radius has dimensions of Length.
     */
    public etomica.units.Dimension getCollisionRadiusDimension() {return etomica.units.Dimension.LENGTH;}

    
    public void intervalAction(IntervalEvent evt) {
        if (pressure >= 0.0) {
            double area = 1.0;
            if (pressure > 0.0) {
                final Vector dimensions = pistonBoundary.dimensions();
                for (int i=0; i<D; i++) {
                    if (i != wallD) {
                        area *= dimensions.x(i);
                    }
                }
            }
            force = pressure/area;
        }
        double t = ((IntegratorHard)evt.getSource()).getTimeStep();
        double a = force/wallMass;
        wallPosition += wallVelocity * t + 0.5*t*a*a;
        wallVelocity += t*a;
        System.out.println("pressure => velocity "+(t*a)+" "+wallVelocity);
    }
 /*   
    public static void main(String[] args) {
        
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        sim.species.setNMolecules(10);
        P1HardBoundary potential = new P1HardBoundary();
        potential.set(sim.phase.speciesMaster);
        potential.setCollisionRadius(0.5*Default.ATOM_SIZE);
  //      sim.phase.setBoundary(sim.space().makeBoundary(Space2D.Boundary.NONE));
        sim.elementCoordinator.go();
        
        Simulation.makeAndDisplayFrame(sim);
    }//end of main
*/
}//end of P1HardBoundary
   
