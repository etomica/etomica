package etomica;
import etomica.statmech.MaxwellBoltzmann;

/**
 * Interaction between a hard sphere and a hard stationary wall.
 * Wall may be horizontal or vertical, as specified by its AtomType object.
 * Vertical wall is at constant x, horizontal wall is at constant y.
 * Wall is taken to be of infinite length.
 *
 * Designed for 2-D, but may work in other dimensional spaces.
 */
 
 //not sure of the difference between this and PotentialHardDiskPiston

public class PotentialHardDiskWall extends Potential implements Potential.Hard, EtomicaElement {
    
    public final String getVersion() {return "PotentialHardDiskWall:01.02.15.0/"+Potential.VERSION;}

    protected double collisionDiameter, collisionRadius;
    
    private boolean isothermal = false;
    private double temperature = Default.TEMPERATURE;

    public PotentialHardDiskWall() {this(Simulation.instance, Default.ATOM_SIZE);}
    
    public PotentialHardDiskWall(double d) {
        this(Simulation.instance, d);
    }
    
    public PotentialHardDiskWall(Simulation sim, double d) {
        super(sim);
        setCollisionDiameter(d);
        ZERO = sim.space.makeTensor();//temporary
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard repulsive potential between a sphere and a wall");
        return info;
    }

 /**
  * Returns infinity if overlap is true, zero otherwise
  */
  public double energy(AtomPair pair) {return overlap(pair) ? Double.MAX_VALUE : 0.0;}
  /**
   * Long-range correction to potential.  Always returns zero for this model.
   */
  public double energyLRC(int n1, int n2, double V) {return 0.0;}

  /**
   * True if perpendicular distance between wall and disk is less than collision radius (diameter/2), false otherwise
   */
  public boolean overlap(AtomPair pair) {
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
        }
        
        double time = Double.MAX_VALUE;
        double dr, dv;
        int i;
        
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
        if(wallType.isVertical()) {i = 0;}
        else {i = 1;}
        
        dr = wall.position(i) - disk.position(i);
        return (Math.abs(dr) < collisionRadius);
    }
  
    /**
     * Time to collision of disk and wall, considering that one or both may be under the influence of a constant force.
     */
    public double collisionTime(AtomPair pair) {
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
        }
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
                
        int i = (((AtomType.Wall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate
        double dr = pair.dr(i);  //dr = atom2 - atom1
        double dv = pair.dv(i);
        if(pair.atom1() == wall) {  //make sure dr = wall - disk
            dr *= -1;
            dv *= -1;
        }
        double a = 0.0;
        if(wall.ia instanceof Integrator.Agent.Forcible) {
            a = ((Integrator.Agent.Forcible)wall.ia).force().component(i) * wall.rm();
        }
        if(disk.ia instanceof Integrator.Agent.Forcible) {
            a -= ((Integrator.Agent.Forcible)disk.ia).force().component(i) * disk.rm();
        }
        //wall or disk has non-zero force
        double time = 0.0;
        if(a != 0.0) {//collision time with acceleration
            time = Double.MAX_VALUE;
            if(Math.abs(dr) <= collisionRadius) {   //inside wall; collision now if approaching, never otherwise
                if(dr*dv < 0.0) return 0.0;  //approaching; collide now
                else {  //separating; move just outside contact and continue to compute collision time
                    moveToContact(i,pair);
                    dr = wall.r.component(i)-disk.r.component(i);
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
            //    else if(-a*dr/(dv*dv) < 1.e-7) {if(dr*dv<0) time = -dr/dv*(1+0.5*dr*a/(dv*dv));} //series expansion for small acceleration; causes missed collisions--may need to handle dr*dv>0 too?
                else {time = (-dv + aSign*Math.sqrt(discrim))/a;} //moving away, accelerating toward (- root is negative)
            }
        }
        else { //force free
            if(dr*dv > 0.0) {return Double.MAX_VALUE;}    //Separating, no collision
            double adr = Math.abs(dr);
            if(adr < collisionRadius) {return 0.0;}            // Inside core and approaching; collision now
            else {return (adr-collisionRadius)/Math.abs(dv);}  //Outside core and approaching; collision at core
        }
        return time;
    }//end of collisionTime

    /**
     * Hard-disk/piston collision dynamics
     */
    public void bump(AtomPair pair) 
    {
        Atom disk;
        Atom wall;
        if(pair.atom2().type instanceof AtomType.Wall) {
           disk = pair.atom1();
           wall = pair.atom2();
        }
        else {
           disk = pair.atom2;
           wall = pair.atom1;
        }
        AtomType.Wall wallType = (AtomType.Wall)wall.type;
    
        int i = (((AtomType.Wall)wall.type).isHorizontal()) ? 1 : 0;  //indicates if collision affects x or y coordinate

        if(wall.isStationary()) {
            if(isothermal) {//specific to 2D
          //      double oldp2 = disk.momentum().squared();
          //      double newp2 = disk.mass()*temperature*parentSimulation().space().D();
          //      disk.momentum().TE(Math.sqrt(newp2/oldp2));
          //      disk.momentum().TE(i,-1.0);
                double px = MaxwellBoltzmann.randomMomentumComponent(temperature,disk.mass());
                double py = MaxwellBoltzmann.randomMomentumComponent(temperature,disk.mass());
                //enforce reflection from wall; new momentum must have opposite sign to old momentum
                if(i==0 && px*disk.momentum().component(i) > 0) px = -px;
                else if(i==1 && py*disk.momentum().component(i) > 0) py = -py;
                disk.momentum().setComponent(0,px);
                disk.momentum().setComponent(1,py);
            }
            else disk.momentum().TE(i,-1.0);
        }
        else {
          double dv = wall.momentum(i)*wall.rm()-disk.momentum(i)*disk.rm();
          double dp = -2.0/(wall.rm() + disk.rm())*dv;
          wall.momentum().PE(i,+dp);  
          disk.momentum().PE(i,-dp);  
        }
        
        moveToContact(i, pair);
        
    }
    
    /**
     * Separates disk and wall so that they are just outside point of contact.
     */
    private void moveToContact(int i, AtomPair pair) {
        double dr = pair.dr(i);  //dr = atom2 - atom1
        double delta = Math.abs(dr) - collisionRadius;
        if(delta < 0.0) {   //inside wall; set apart to contact point
            double mult = (dr > 0.0) ? -2.0 : +2.0;
            pair.atom2.r.PE(i,mult*delta);
            pair.atom1.r.PE(i,-mult*delta);
            double diff = pair.atom2.r.component(i)-pair.atom1.r.component(i);
        }
    }    
    
    /**
     * Always returns zero (not yet implemented)
     */
  public double lastCollisionVirial() {return 0.0;}
  final Space.Tensor ZERO; //temporary
  public Space.Tensor lastCollisionVirialTensor() {return ZERO;}


    /**
     * Accessor method for collision diameter.
     * Wall and disk begin to overlap when center of disk is one half this distance from the wall,
     * as measured perpendicularly from the wall.
     */
    public double getCollisionDiameter() {return collisionDiameter;}
    /**
     * Accessor method for collision diameter.
     * Wall and disk begin to overlap when center of disk is one half this distance from the wall,
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