package etomica.potential;

import etomica.Atom;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Simulation;
import etomica.Space;
import etomica.Integrator.Forcible;
import etomica.atom.AtomTypeWall;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Interaction between a hard sphere and a hard stationary wall.
 * Wall may be horizontal or vertical, as specified by its AtomType object.
 * Vertical wall is at constant x, horizontal wall is at constant y.
 * Wall is taken to be of infinite length.
 *
 * Designed for 2-D, but may work in other dimensional spaces.
 *
 * @author David Kofke
 */
 
public class P2HardSphereWall extends Potential2 implements PotentialHard {
    
    protected double collisionDiameter, collisionRadius;
    private double lastCollisionVirial;
    
    private boolean isothermal = false;
    private double temperature = Default.TEMPERATURE;

    public P2HardSphereWall() {this(Simulation.getDefault().space, Default.ATOM_SIZE);}
    
    public P2HardSphereWall(double d) {
        this(Simulation.getDefault().space, d);
    }
    
    public P2HardSphereWall(Space space) {
        this(space, Default.ATOM_SIZE);
    }
    
    public P2HardSphereWall(Space space, double d) {
        super(space);
        setCollisionDiameter(d);
        ZERO = space.makeTensor();//temporary
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard repulsive potential between a sphere and a wall");
        return info;
    }

    /**
     * Returns infinity if overlap is true, zero otherwise
     */
    public double energy(Atom[] pair) {return overlap(pair) ? Double.MAX_VALUE : 0.0;}
    
    public double energyChange() {return 0.0;}

    /**
     * True if perpendicular distance between wall and sphere is less than
     * collision radius (diameter/2), false otherwise
     */
    public boolean overlap(Atom[] pair) {
        cPair.reset(pair[0].coord,pair[1].coord);
        Atom sphere;
        Atom wall;
        if(pair[1].type instanceof AtomTypeWall) {
           sphere = pair[0];
           wall = pair[1];
        }
        else {
           sphere = pair[1];
           wall = pair[0];
        }
        
        double time = Double.MAX_VALUE;
        double dr, dv;
        int i;
        
        AtomTypeWall wallType = (AtomTypeWall)wall.type;
        if(wallType.isVertical()) {i = 0;}
        else {i = 1;}
        
        dr = wall.coord.position().x(i) - sphere.coord.position().x(i);
        return (Math.abs(dr) < collisionRadius);
    }
  
    /**
     * Time to collision of sphere and wall, considering that one or both may
     * be under the influence of a constant force.
     */
    public double collisionTime(Atom[] pair) {
    	cPair.reset(pair[0].coord,pair[1].coord);
//    	System.out.println("p2hardspherewall "+pair.toString());
        Atom sphere;
        Atom wall;
        if(pair[1].type instanceof AtomTypeWall) {
           sphere = pair[0];
           wall = pair[1];
        }
        else {
           sphere = pair[1];
           wall = pair[0];
        }
        AtomTypeWall wallType = (AtomTypeWall)wall.type;
                
        int i = (((AtomTypeWall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
        double dr = cPair.dr(i);  //dr = atom2 - atom1
		cPair.resetV();
        double dv = cPair.dv(i);
        if(pair[0] == wall) {  //make sure dr = wall - sphere
            dr *= -1;
            dv *= -1;
        }
        double a = 0.0;
        if(wall.ia instanceof Integrator.Forcible  && !wall.coord.isStationary()) {
            a = ((Integrator.Forcible)wall.ia).force().x(i) * wall.type.rm();
        }
        if(sphere.ia instanceof Integrator.Forcible  && !sphere.coord.isStationary()) {
            a -= ((Integrator.Forcible)sphere.ia).force().x(i) * sphere.type.rm();
        }
        //wall or sphere has non-zero force
        double time = 0.0;
        if(a != 0.0) {//collision time with acceleration
            time = Double.MAX_VALUE;
            if(Math.abs(dr) <= collisionRadius) {   //inside wall; collision now if approaching, never otherwise
                if(dr*dv < 0.0) return 0.0;  //approaching; collide now
                else {  //separating; move just outside contact and continue to compute collision time
                    moveToContact(i,pair);
                    dr = wall.coord.position().x(i)-sphere.coord.position().x(i);
                }
            }  
            dr += (dr > 0.0) ? -collisionRadius : +collisionRadius;
            double discrim = dv*dv - 2*a*dr;
            if(discrim > 0) {
                boolean adr = (a*dr > 0); //accelerating away from each other
                boolean adv = (a*dv > 0); //moving away from each other
                int aSign = (a > 0) ? +1 : -1;
                if(adr && adv) {time = Double.MAX_VALUE;}  //moving and accelerating away; no collision
                else if(adr) {time = (-dv - aSign*Math.sqrt(discrim))/a;} //moving toward, accelerating away (- root is positive)
                // else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));}
                // series expansion for small acceleration; causes missed
                // collisions--may need to handle dr*dv>0 too?
                else {time = (-dv + aSign*Math.sqrt(discrim))/a;} //moving away, accelerating toward (- root is negative)
            }
        }
        else { //force free
            if(dr*dv > 0.0) {return Double.MAX_VALUE;}    //Separating, no collision
            double adr = Math.abs(dr);
            if(Default.FIX_OVERLAP && adr < collisionRadius) {return 0.0;}            // Inside core and approaching; collision now
            else {return (adr-collisionRadius)/Math.abs(dv);}  //Outside core and approaching; collision at core
        }
        return time;
    }//end of collisionTime

    /**
     * Hard-sphere/piston collision dynamics
     */
    public void bump(Atom[] pair) 
    {
    	cPair.reset(pair[0].coord,pair[1].coord);
        Atom sphere;
        Atom wall;
        if(pair[1].type instanceof AtomTypeWall) {
           sphere = pair[0];
           wall = pair[1];
        }
        else {
           sphere = pair[1];
           wall = pair[0];
        }
        AtomTypeWall wallType = (AtomTypeWall)wall.type;
    
        int i = (((AtomTypeWall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
        cPair.resetV();
        double pOld = sphere.coord.momentum().x(i);
        
        if(wall.coord.isStationary()) {
         /*   if(isothermal) {//specific to 2D
          //      double oldp2 = sphere.momentum().squared();
          //      double newp2 = sphere.mass()*temperature*parentSimulation().space().D();
          //      sphere.momentum().TE(Math.sqrt(newp2/oldp2));
          //      sphere.momentum().TE(i,-1.0);
                double px = MaxwellBoltzmann.randomMomentumComponent(temperature,sphere.coord.mass());
                double py = MaxwellBoltzmann.randomMomentumComponent(temperature,sphere.coord.mass());
                //enforce reflection from wall; new momentum must have opposite sign to old momentum
                if(i==0 && px*sphere.coord.momentum().x(i) > 0) px = -px;
                else if(i==1 && py*sphere.coord.momentum().x(i) > 0) py = -py;
                sphere.coord.momentum().setx(0,px);
                sphere.coord.momentum().setx(1,py);
            }
            else {
                sphere.coord.momentum().TE(i,-1.0);
    //            wallType.pAccumulator += 2*sphere.momentum().x(i);
            }*/
                sphere.coord.momentum().TE(i,-1.0);
 //               wallType.pAccumulator += 2.*sphere.coord.momentum(i);
            
        }
        else {
          double dv = wall.coord.momentum().x(i)*wall.coord.rm()-sphere.coord.momentum(i)*sphere.coord.rm();
          double dp = -2.0/(wall.coord.rm() + sphere.coord.rm())*dv;
 //         wall.coord.momentum().PE(i,+dp);  
          sphere.coord.momentum().PE(i,-dp);
  //        wallType.pAccumulator -= dp;
        }
        
        if(wall.coord.isStationary()) {//randomize x-y components of sphere momentum
            if(isothermal) {
                double oldp2 = sphere.coord.momentum().squared();
                double newp2 = sphere.coord.mass()*temperature*space.D();
                sphere.coord.momentum().TE(Math.sqrt(newp2/oldp2));
            }

            double px = sphere.coord.momentum(0);
            double py = sphere.coord.momentum(1);
            double xSign = px/Math.abs(px);
            double ySign = py/Math.abs(py);
            double p2Tot = px*px + py*py;
            px = xSign*Math.sqrt(Simulation.random.nextDouble()*p2Tot);
            py = ySign*Math.sqrt(p2Tot - px*px);
            sphere.coord.momentum().setX(0,px);
            sphere.coord.momentum().setX(1,py);
        }
        double deltaP = sphere.coord.momentum(i) - pOld;
        if(!wall.coord.isStationary()) wall.coord.momentum().PE(i,-deltaP);
        
        wallType.pAccumulator += deltaP;
        moveToContact(i, pair);
        
    }
    
    /**
     * Separates sphere and wall so that they are just outside point of contact.
     */
    // assume that cPair is already set up since this method is private
    private void moveToContact(int i, Atom[] pair) {
        double dr = cPair.dr(i);  //dr = atom2 - atom1
        double mult = (dr > 0.0) ? -1 : +1;
        double delta = Math.abs(dr) - collisionRadius;
        if(delta < 0.0) {   //inside wall; set apart to contact point
 //           double mult = (dr > 0.0) ? -2.0 : +2.0;
            Vector r1 = pair[0].coord.position();
            Vector r2 = pair[1].coord.position();
            if(pair[1].coord.isStationary()) 
                r1.setX(i,r2.x(i)+mult*(collisionRadius+1e-6));
            else if(pair[0].coord.isStationary())
                r2.setX(i,r1.x(i)-mult*(collisionRadius+1e-6));
            else {
                double mid = (r1.x(i)+r2.x(i))*0.5;
                r1.setX(i,mid+0.5*mult*(collisionRadius+1e-6));
                r2.setX(i,mid-0.5*mult*(collisionRadius+1e-6));
     //           pair.atom2.r.PE(i,delta);
     //           pair.atom1.r.PE(i,-delta);
            }
     //       delta = Math.abs(pair.atom2.r.x(i)-pair.atom1.r.x(i)) - collisionRadius;
     //       System.out.println(delta);
        }
    }    
    
    /**
     * Always returns zero (not yet implemented)
     */
    public double lastCollisionVirial() {return 0.0;}
    final Tensor ZERO; //temporary
    public Tensor lastCollisionVirialTensor() {return ZERO;}


    /**
     * Accessor method for collision diameter.
     * Wall and sphere begin to overlap when center of sphere is one half this distance from the wall,
     * as measured perpendicularly from the wall.
     */
    public double getCollisionDiameter() {return collisionDiameter;}
    /**
     * Accessor method for collision diameter.
     * Wall and sphere begin to overlap when center of sphere is one half this distance from the wall,
     * as measured perpendicularly from the wall.
     */
    public final void setCollisionDiameter(double c) {
        collisionDiameter = c;
        collisionRadius = 0.5*c;
    }
    
    public void setIsothermal(boolean b) {isothermal = b;}
    public boolean isIsothermal() {return isothermal;}
    public void setTemperature(double t) {temperature = t;}
    public double getTemperature() {return temperature;}
    public etomica.units.Dimension getTemperatureDimension() {
        return etomica.units.Dimension.TEMPERATURE;
    }
}
