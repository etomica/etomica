package etomica.potential;

import etomica.Atom;
import etomica.Debug;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Space;
import etomica.space.Boundary;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Potential that places hard repulsive walls that move and 
 * accelerate subject to an external force field (pressure).
 */
 
public class P1HardMovingBoundary extends Potential1 implements PotentialHard {
    
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
    
    public double collisionTime(Atom[] atoms, double falseTime) {
    	atom = atoms[0];
        double dr = atom.coord.position().x(wallD) - wallPosition;
        double dv = ((ICoordinateKinetic)atom.coord).velocity().x(wallD) - wallVelocity;
        dr += dv*falseTime;
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
        double a = -force/wallMass;   // atom acceleration - wall acceleration
        dv += a*falseTime;
        dr += 0.5*a*falseTime*falseTime;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(atoms)) {
            System.out.println(dr+" "+dv+" "+falseTime+" "+atom);
            System.out.println(((ICoordinateKinetic)atom.coord).velocity().x(wallD));
            System.out.println(atom.coord.position().x(wallD));
        }
        double t = Double.POSITIVE_INFINITY;
        double discr = -1.0;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(atoms)) {
            System.out.println(atom);
        }
        if (dr*dv < 0.0 || dr*a < 0.0) {
            // either moving toward or accelerating toward each other
            if (Math.abs(dr) < collisionRadius && dr*dv < 0.0) {
                throw new RuntimeException("overlap "+atom+" "+dr+" "+dv+" "+a);
            }
            double drc;
            if (dr>0.0) {
                drc = dr - collisionRadius;
            }
            else {
                drc = dr + collisionRadius;
            }
            discr = dv*dv - 2.0*a*drc;
            if (discr >= 0.0) {
                discr = Math.sqrt(discr);
                if (dr*dv > 0.0 && dr*a < 0.0) {
                    t = -dv/a + discr/Math.abs(a);
                }
                else if (dr*dv < 0.0 && dr*a > 0.0) {
                    t = -dv/a - discr/Math.abs(a);
                }
                else if (dr*dv < 0.0 && dr*a < 0.0) {
                    t = -dv/a + discr/Math.abs(a);
                }
                else {
                    throw new RuntimeException("oops");
                }
            }
        }
/*        if (dv * dr < 0.0) {
            if (dv > 0.0) {
                dr = dr - collisionRadius;
            }
            else {
                dr = -dr + collisionRadius;
            }
            t = dr / dv;
        }*/
        if (Default.FIX_OVERLAP && t<0.0) t = 0.0;
        if (Debug.ON && (t<0.0 || Debug.DEBUG_NOW && Debug.anyAtom(atoms))) {
            System.out.println(a+" "+dr+" "+dv+" "+discr+" "+t+" "+atom);
            if (t<0) throw new RuntimeException("foo");
        }
        return t + falseTime;
    }
                
    public void bump(Atom[] a, double falseTime) {
    	atom = a[0];
        double r = atom.coord.position().x(wallD);
        Vector v = ((ICoordinateKinetic)atom.coord).velocity();
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
        double dp = 2.0/(1/wallMass + atom.type.rm())*(trueWallVelocity-v.x(wallD));
        lastVirial = -Math.abs(dp)*collisionRadius;
        if (Debug.ON) {
            double trueWallPosition = wallPosition + wallVelocity*falseTime + 0.5*falseTime*falseTime*(force/wallMass);
//            if (Math.abs(Math.abs(trueWallPosition-(r+v.x(wallD)*falseTime)) - collisionRadius) > 1.e-9*collisionRadius) {
                System.out.println("bork at "+falseTime+" ! "+atom+" "+(r+v.x(wallD)*falseTime)+" "+v.x(wallD));
                System.out.println("wall bork! "+trueWallPosition+" "+trueWallVelocity+" "+force);
                System.out.println("dr bork! "+((r+v.x(wallD)*falseTime)-trueWallPosition)+" "+collisionRadius);
                System.out.println(atom.coord.position().x(wallD));
//                throw new RuntimeException("bork!");
//            }
        }
        v.setX(wallD,v.x(wallD)+dp*atom.type.rm());
        atom.coord.position().setX(wallD,r-dp*atom.type.rm()*falseTime);
        wallVelocity -= dp/wallMass;
        wallPosition += dp/wallMass*falseTime;
        
        double dr = atom.coord.position().x(wallD) - wallPosition;
        double dv = ((ICoordinateKinetic)atom.coord).velocity().x(wallD) - wallVelocity;
        dr += dv*falseTime;
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
        double acc = -force/wallMass;   // atom acceleration - wall acceleration
        dv += acc*falseTime;
        dr += 0.5*acc*falseTime*falseTime;
        System.out.println("*****"+atom+" "+dr+" "+dv+" "+acc);
        if (dv < 0) {
            throw new RuntimeException("dv must be positive");
        }
        
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

    
    public void advanceAcrossTimeStep(double tStep) {
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
        double a = force/wallMass;
        wallPosition += wallVelocity * tStep + 0.5*tStep*tStep*a;
        wallVelocity += tStep*a;
        System.out.println("pressure => velocity "+(tStep*a)+" "+wallVelocity+" "+wallPosition+" "+tStep);
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
   
