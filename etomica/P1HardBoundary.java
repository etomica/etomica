//includes main method
package etomica;

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
 
public class P1HardBoundary extends Potential1 implements Potential1.Hard {
    
    public String getVersion() {return "P1HardBoundary:01.06.29/"+Potential1.VERSION;}

    private double collisionRadius = 0.0;
    private boolean isothermal = false;
    private double temperature;
    private final int D;
    private double []lastVirial = new double[2];
    
    public P1HardBoundary() {
        this(Simulation.instance.hamiltonian.potential);
    }
    
    public P1HardBoundary(PotentialGroup parent) {
        super(parent);
        temperature = Default.TEMPERATURE;
        D = parent.parentSimulation().space.D();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Hard repulsive potential at the phase boundaries");
        return info;
    }
    
    public double energy(Atom a) {
        double e = 0.0;
        Space.Vector dimensions = a.node.parentPhase().dimensions();
        double collisionRadiusSquared = collisionRadius*collisionRadius;
        double rx = a.coord.position(0);
        double ry = a.coord.position(1);
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
     
    public double collisionTime(Atom a) {
        Space.Vector r = a.coord.position();
        Space.Vector p = a.coord.momentum();
        Space.Vector dimensions = a.node.parentPhase().dimensions();
        double tmin = Double.MAX_VALUE;
        for(int i=r.length(); i>=0; i--) {
            double px = p.x(i);
            if(px == 0.0) continue;
            double rx = r.x(i);
            double dx = dimensions.x(i);
            double t = (px > 0.0) ? (dx - rx - collisionRadius)/px : (-rx + collisionRadius)/px;
            if(t < tmin) tmin = t;
        }
        return a.coord.mass()*tmin;
    }
                
//    public void bump(IntegratorHardAbstract.Agent agent) {
//        Atom a = agent.atom();
    public void bump(AtomPair pair) {
        bump(pair.atom1());
    }
    public void bump(Atom a) {
        Space.Vector r = a.coord.position();
        Space.Vector p = a.coord.momentum();
        Space.Vector dimensions = a.node.parentPhase().dimensions();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for(int i=r.length(); i>=0; i--) {
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
        p.setX(imin,-p.x(imin)); //multiply momentum component by -1
        if(isothermal) {p.TE(Math.sqrt(D*temperature*a.coord.mass()/p.squared()));}
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
   
