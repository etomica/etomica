//includes main method
package etomica.potential;

import etomica.Atom;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.space.ICoordinateKinetic;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 * Potential that places hard repulsive walls coinciding with the
 * boundary of the phase, which is assumed to be rectangular in shape.
 *
 * @author David Kofke
 */
 
 //extends PotentialAbstract instead of Potential1HardAbstract because potential depends
 //on property (boundary.dimensions) of phase, and this is more readily available to
 //the Agent than to the parent potential.  
 //perhaps Potential1HardAbstract should be redesigned to permit easier access to features
 //of the phase by the parent potential, since this is probably a common situation for
 //one-body potentials
 
public class P1HardBoundary extends Potential1 implements PotentialHard {
    
    private double collisionRadius = 0.0;
    private boolean isothermal = false;
    private double temperature;
    private Atom atom;
    private double lastVirial;
    
    public P1HardBoundary() {
        this(Simulation.getDefault().space);
    }
    
    public P1HardBoundary(Space space) {
        super(space);
        temperature = Default.TEMPERATURE;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard repulsive potential at the phase boundaries");
        return info;
    }
    
    public double energy(Atom[] a) {
        atom = a[0];
        double e = 0.0;
        Vector dimensions = boundary.dimensions();
        double collisionRadiusSquared = collisionRadius*collisionRadius;
        double rx = atom.coord.position().x(0);
        double ry = atom.coord.position().x(1);
        double dx0_2 = rx*rx;
        double dy0_2 = ry*ry;
        double dx1_2 = (rx - dimensions.x(0))*(rx - dimensions.x(0));
        double dy1_2 = (ry - dimensions.x(1))*(ry - dimensions.x(1));
        if((dx0_2 < collisionRadiusSquared)||(rx < 0.0)||(rx > dimensions.x(0) )||(dx1_2 < collisionRadiusSquared)||
            (dy0_2 < collisionRadiusSquared)||(ry < 0.0)||(ry > dimensions.x(1) )||(dy1_2 < collisionRadiusSquared)){
            e = Double.POSITIVE_INFINITY;
        }else{e = 0.0;}
        return e;
    }
     
    public double energyChange() {return 0.0;}
    
    public double collisionTime(Atom[] a, double falseTime) {
    	atom = a[0];
        Vector r = atom.coord.position();
        Vector v = ((ICoordinateKinetic)atom.coord).velocity();
        r.PEa1Tv1(falseTime,v);
        Vector dimensions = boundary.dimensions();
        double tmin = Double.POSITIVE_INFINITY;
        for(int i=r.length()-1; i>=0; i--) {
            double vx = v.x(i);
            if(vx == 0.0) continue;
            double rx = r.x(i);
            double dx = dimensions.x(i);
            double t = (vx > 0.0) ? (dx - rx - collisionRadius)/vx : (-rx + collisionRadius)/vx;
            if(t < tmin) tmin = t;
        }
        if (Default.FIX_OVERLAP && tmin<0.0) tmin = 0.0;
        return tmin + falseTime;
    }
                
//    public void bump(IntegratorHard.Agent agent) {
//        Atom a = agent.atom();
    public void bump(Atom[] a, double falseTime) {
    	atom = a[0];
        Vector r = atom.coord.position();
        Vector v = ((ICoordinateKinetic)atom.coord).velocity();
        r.PEa1Tv1(falseTime,v);
        Vector dimensions = boundary.dimensions();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for(int i=r.length()-1; i>=0; i--) {
            double rx = r.x(i);
            double vx = v.x(i);
            double dx = dimensions.x(i);
            double del = (vx > 0.0) ? Math.abs(dx - rx - collisionRadius) : Math.abs(-rx + collisionRadius);
            if(del < delmin) {
                delmin = del;
                imin = i;
            }
        }
        lastVirial = atom.type.getMass()*2.0*v.x(imin)*collisionRadius;
        v.setX(imin,-v.x(imin));
        // dv = 2*NewVelocity
        double newP = atom.coord.position().x(imin) - falseTime*v.x(imin)*2.0;
        atom.coord.position().setX(imin,newP);
    }//end of bump
    
    public void setIsothermal(boolean b) {isothermal = b;}
    public boolean isIsothremal() {return isothermal;}
    
    public void setTemperature(double t) {temperature = t;}
    public double getTemperature() {return temperature;}
    public etomica.units.Dimension getTemperatureDimension() {return etomica.units.Dimension.TEMPERATURE;}
        
    /**
     * not yet implemented
     */
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
   
