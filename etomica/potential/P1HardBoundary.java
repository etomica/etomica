//includes main method
package etomica.potential;

import etomica.Atom;
import etomica.Default;
import etomica.EtomicaInfo;
import etomica.Simulation;
import etomica.Space;
import etomica.Space.Tensor;
import etomica.Space.Vector;

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
    private final int D;
    private double[] lastVirial = new double[2];
    private Atom atom;
    
    public P1HardBoundary() {
        this(Simulation.getDefault().space);
    }
    
    public P1HardBoundary(Space space) {
        super(space);
        temperature = Default.TEMPERATURE;
        D = space.D();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard repulsive potential at the phase boundaries");
        return info;
    }
    
    public double energy(Atom[] a) {
        atom = a[0];
        double e = 0.0;
        Space.Vector dimensions = atom.node.parentPhase().dimensions();
        double collisionRadiusSquared = collisionRadius*collisionRadius;
        double rx = atom.coord.position(0);
        double ry = atom.coord.position(1);
        double dx0_2 = rx*rx;
        double dy0_2 = ry*ry;
        double dx1_2 = (rx - dimensions.x(0))*(rx - dimensions.x(0));
        double dy1_2 = (ry - dimensions.x(1))*(ry - dimensions.x(1));
        if((dx0_2 < collisionRadiusSquared)||(rx < 0.0)||(rx > dimensions.x(0) )||(dx1_2 < collisionRadiusSquared)||
            (dy0_2 < collisionRadiusSquared)||(ry < 0.0)||(ry > dimensions.x(1) )||(dy1_2 < collisionRadiusSquared)){
            e = Double.MAX_VALUE;
        }else{e = 0.0;}
        return e;
    }
     
    public double energyChange() {return 0.0;}
    
    public double collisionTime(Atom[] a, double falseTime) {
    	atom = a[0];
        Space.Vector r = atom.coord.truePosition(falseTime);
        Space.Vector p = atom.coord.momentum();
        Space.Vector dimensions = atom.node.parentPhase().dimensions();
        double tmin = Double.POSITIVE_INFINITY;
        for(int i=r.length()-1; i>=0; i--) {
            double px = p.x(i);
            if(px == 0.0) continue;
            double rx = r.x(i);
            double dx = dimensions.x(i);
            double t = (px > 0.0) ? (dx - rx - collisionRadius)/px : (-rx + collisionRadius)/px;
            if(t < tmin) tmin = t;
        }
        if (Default.FIX_OVERLAP && tmin<0.0) tmin = 0.0;
        return atom.coord.mass()*tmin + falseTime;
    }
                
//    public void bump(IntegratorHard.Agent agent) {
//        Atom a = agent.atom();
    public void bump(Atom[] a, double falseTime) {
    	atom = a[0];
        Space.Vector r = atom.coord.truePosition(falseTime);
        Space.Vector p = atom.coord.momentum();
        Space.Vector dimensions = atom.node.parentPhase().dimensions();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for(int i=r.length()-1; i>=0; i--) {
            double rx = r.x(i);
            double px = p.x(i);
            double dx = dimensions.x(i);
            double del = (px > 0.0) ? Math.abs(dx - rx - collisionRadius) : Math.abs(-rx + collisionRadius);
            if(del < delmin) {
                delmin = del;
                imin = i;
            }
        }
//        pAccumulator += 2*Math.abs(p.x(imin));
        lastVirial[0] = 0.0;
        lastVirial[1] = 0.0;
        if(imin == 0){
            if(p.x(0) > 0.0){
                lastVirial[0] = 2.0*p.x(0)*(dimensions.x(0)-r.x(0));
            }else{
                lastVirial[1] = 2.0*p.x(0)*(0.0-r.x(0));
            }
        }
        Space.Vector newP = Space.makeVector(D);
        newP.E(p);
        newP.setX(imin,-p.x(imin)); //multiply momentum component by -1
        if(isothermal) {
            //XXX Evil.  If anything, atoms should get velocity from Maxwell-Boltzmann
            //distribution such that the atom is moving away from the wall
            newP.TE(Math.sqrt(D*temperature*atom.coord.mass()/newP.squared()));
        }
        atom.coord.trueAccelerateTo(newP,falseTime);
    }//end of bump
    
    public void setIsothermal(boolean b) {isothermal = b;}
    public boolean isIsothremal() {return isothermal;}
    
    public void setTemperature(double t) {temperature = t;}
    public double getTemperature() {return temperature;}
    public etomica.units.Dimension getTemperatureDimension() {return etomica.units.Dimension.TEMPERATURE;}
        
    public double lastCollisionVirial(int i) {return lastVirial[i];}

    /**
     * not yet implemented
     */
    public double lastCollisionVirial() {return lastVirial[0];}
    //public double lastCollisionVirial() {return Double.NaN;}
    
    /**
     * not yet implemented.
     */
    public Space.Tensor lastCollisionVirialTensor() {return null;}
    
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
   
