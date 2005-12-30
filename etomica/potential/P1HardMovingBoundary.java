package etomica.potential;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.graphics.Drawable;
import etomica.space.Boundary;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Length;
import etomica.util.Debug;

/**
 * Potential that places hard repulsive walls that move and 
 * accelerate subject to an external force field (pressure).
 */
 
public class P1HardMovingBoundary extends Potential1 implements PotentialHard, Drawable {
    
    /**
     * Constructor for a hard moving (and accelerating) boundary.
     * @param space
     * @param wallDimension dimension which the wall is perpendicular to
     */
    public P1HardMovingBoundary(Space space, Boundary boundary, int wallDimension, double mass,
            boolean ignoreOverlap) {
        super(space);
        D = space.D();
        wallD = wallDimension;
        wallPosition = 0.0;
        setMass(mass);
        force = 0.0;
        pistonBoundary = boundary;
        virialSum = 0.0;
        this.ignoreOverlap = ignoreOverlap;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard moving wall");
        return info;
    }
    
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
    public double getPressure() {
        return pressure;
    }
    
    public void setStationary(boolean b) {
        if(b) {
            wallMass = Double.POSITIVE_INFINITY;
            wallVelocity = 0.0;
        } else {
            wallMass = setWallMass;
        }
    }
    
    public boolean isStationary() {
        return Double.isInfinite(wallMass);
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
        setWallMass = mass;
    }
    
    public double energy(AtomSet a) {
        double dx = ((Atom)a).coord.position().x(wallD) - wallPosition;
        if (dx*dx < collisionRadius*collisionRadius) {
            return Double.POSITIVE_INFINITY;
        }
        return 0.0;
    }
     
    public double energyChange() {return 0.0;}
    
    public double collisionTime(AtomSet atoms, double falseTime) {
        double dr = ((Atom)atoms).coord.position().x(wallD) - wallPosition;
        double dv = ((ICoordinateKinetic)((Atom)atoms).coord).velocity().x(wallD) - wallVelocity;
        dr += dv*falseTime;
        if (pressure >= 0.0) {
            double area = 1.0;
            if (pressure > 0.0) {
                final Vector dimensions = pistonBoundary.getDimensions();
                for (int i=0; i<D; i++) {
                    if (i != wallD) {
                        area *= (dimensions.x(i)-collisionRadius*2.0);
                    }
                }
            }
            force = pressure*area;
        }
        double a = -force/wallMass;   // atom acceleration - wall acceleration
        dv += a*falseTime;
        dr += 0.5*a*falseTime*falseTime;
        if (Debug.ON && Debug.DEBUG_NOW && Debug.anyAtom(atoms)) {
            System.out.println(dr+" "+dv+" "+falseTime+" "+atoms);
            System.out.println(((ICoordinateKinetic)((Atom)atoms).coord).velocity().x(wallD));
            System.out.println(((Atom)atoms).coord.position().x(wallD));
        }
        double t = Double.POSITIVE_INFINITY;
        double discr = -1.0;
        if (dr*dv < 0.0 || dr*a < 0.0) {
            // either moving toward or accelerating toward each other
            if ((Debug.ON || ignoreOverlap) && Math.abs(dr) < collisionRadius && dr*dv < 0.0) {
                if (ignoreOverlap) return falseTime;
                throw new RuntimeException("overlap "+atoms+" "+dr+" "+dv+" "+a);
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
                if (dr*a < 0.0) {
                    t = -dv/a + discr/Math.abs(a);
                }
                else if (a == 0.0) {
                    if (dr*dv < 0.0) t = -drc/dv;
                }
                else if (dr*dv < 0.0 && dr*a > 0.0) {
                    t = -dv/a - discr/Math.abs(a);
                } else {
                    throw new RuntimeException("oops");
                }
            }
        }
        if (ignoreOverlap && t<0.0) t = 0.0;
        if (Debug.ON && (t<0.0 || Debug.DEBUG_NOW && Debug.anyAtom(atoms))) {
            System.out.println(atoms+" "+a+" "+dr+" "+dv+" "+discr+" "+t+" "+(t+falseTime)+" "+(((Atom)atoms).coord.position().x(wallD)+((ICoordinateKinetic)((Atom)atoms).coord).velocity().x(wallD)*(t+falseTime))+" "+(wallPosition+wallVelocity*(t+falseTime)-0.5*a*(t+falseTime)*(t+falseTime)));
            if (t<0) throw new RuntimeException("foo");
        }
        return t + falseTime;
    }
                
    public void bump(AtomSet a, double falseTime) {
        double r = ((Atom)a).coord.position().x(wallD);
        Vector v = ((ICoordinateKinetic)((Atom)a).coord).velocity();
        if (pressure >= 0.0) {
            double area = 1.0;
            if (pressure > 0.0) {
                final Vector dimensions = pistonBoundary.getDimensions();
                for (int i=0; i<D; i++) {
                    if (i != wallD) {
                        area *= (dimensions.x(i)-collisionRadius*2.0);
                    }
                }
            }
            force = pressure*area;
        }
        double trueWallVelocity = wallVelocity + falseTime*force/wallMass;
        if (Debug.ON) {
            double trueWallPosition = wallPosition + wallVelocity*falseTime + 0.5*falseTime*falseTime*(force/wallMass);
            if (Math.abs(Math.abs(trueWallPosition-(r+v.x(wallD)*falseTime)) - collisionRadius) > 1.e-7*collisionRadius) {
                System.out.println("bork at "+falseTime+" ! "+a+" "+(r+v.x(wallD)*falseTime)+" "+v.x(wallD));
                System.out.println("wall bork! "+trueWallPosition+" "+trueWallVelocity+" "+force);
                System.out.println("dr bork! "+((r+v.x(wallD)*falseTime)-trueWallPosition)+" "+collisionRadius);
                System.out.println(((Atom)a).coord.position().x(wallD));
                throw new RuntimeException("bork!");
            }
        }
        double dp = 2.0/(1/wallMass + ((AtomTypeLeaf)((Atom)a).type).rm())*(trueWallVelocity-v.x(wallD));
        virialSum += dp;
        v.setX(wallD,v.x(wallD)+dp*((AtomTypeLeaf)((Atom)a).type).rm());
        ((Atom)a).coord.position().setX(wallD,r-dp*((AtomTypeLeaf)((Atom)a).type).rm()*falseTime);
        wallVelocity -= dp/wallMass;
        wallPosition += dp/wallMass*falseTime;
        
    }
    
    public double lastWallVirial() {
        double area = 1.0;
        final Vector dimensions = pistonBoundary.getDimensions();
        for (int i=0; i<D; i++) {
            if (i != wallD) {
                area *= (dimensions.x(i)-collisionRadius*2.0);
            }
        }
        double s = virialSum / area;
        virialSum = 0.0;
        return s;
    }
        
    
    public double lastCollisionVirial() {return 0;}
    
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
    public etomica.units.Dimension getCollisionRadiusDimension() {return etomica.units.Length.DIMENSION;}

    
    public void advanceAcrossTimeStep(double tStep) {
        if (pressure >= 0.0) {
            double area = 1.0;
            final Vector dimensions = pistonBoundary.getDimensions();
            for (int i=0; i<D; i++) {
                if (i != wallD) {
                    area *= (dimensions.x(i)-collisionRadius*2.0);
                }
            }
            force = pressure*area;
        }
        double a = force/wallMass;
        wallPosition += wallVelocity * tStep + 0.5*tStep*tStep*a;
        wallVelocity += tStep*a;
//        System.out.println("pressure => velocity "+a+" "+wallVelocity+" "+wallPosition+" "+tStep);
    }
    
    public void setThickness(double t) {
        thickness = t;
    }
    
    public void draw(java.awt.Graphics g, int[] origin, double toPixel) {
        g.setColor(java.awt.Color.gray);
        int xP = origin[0] + (wallD==0 ? (int)((wallPosition-thickness)*toPixel) : 0);
        int yP = origin[1] + (wallD==1 ? (int)((wallPosition-thickness)*toPixel) : 0);
        int t = Math.max(1,(int)(thickness*toPixel));
        int wP = wallD==0 ? t : (int)(toPixel*pistonBoundary.getDimensions().x(0));
        int hP = wallD==1 ? t : (int)(toPixel*pistonBoundary.getDimensions().x(1));
        g.fillRect(xP,yP,wP,hP);
    }
    
    private double collisionRadius = 0.0;
    private final int D;
    private final int wallD;
    private double wallPosition;
    private double wallVelocity;
    private double wallMass;
    private double setWallMass;
    private double force;
    private double pressure;
    private final Boundary pistonBoundary;
    private double thickness = 0.0;
    private double virialSum;
    private boolean ignoreOverlap;
    
}//end of P1HardBoundary
   
