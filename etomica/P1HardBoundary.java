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
 
public class P1HardBoundary extends Potential1 implements Potential1Hard {
    
    public String getVersion() {return "P1HardBoundary:01.06.29/"+Potential1.VERSION;}

    private double collisionRadius = 0.0;
    
    P1HardBoundary(Simulation sim) {
        super(sim);
    }
    
    public double energy(Atom a) {return 0.0;}
     
    public double collisionTime(Atom a) {
        Space.Vector r = a.coordinate().position();
        Space.Vector p = a.coordinate().momentum();
        Space.Vector dimensions = a.parentPhase().dimensions();
        double tmin = Double.MAX_VALUE;
        for(int i=r.length(); i>=0; i--) {
            double rx = r.component(i);
            double px = p.component(i);
            double dx = dimensions.component(i);
            double t = (px > 0.0) ? (dx - rx - collisionRadius)/px : (-rx + collisionRadius)/px;
            if(t < tmin) tmin = t;
        }
        return a.mass()*tmin;
    }
                
//    public void bump(IntegratorHardAbstract.Agent agent) {
//        Atom a = agent.atom();
    public void bump(AtomPair pair) {
        bump(pair.atom1());
    }
    public void bump(Atom a) {
        Space.Vector r = a.coordinate().position();
        Space.Vector p = a.coordinate().momentum();
        Space.Vector dimensions = a.parentPhase().dimensions();
        double delmin = Double.MAX_VALUE;
        int imin = 0;
        //figure out which component is colliding
        for(int i=r.length(); i>=0; i--) {
            double rx = r.component(i);
            double px = p.component(i);
            double dx = dimensions.component(i);
            double del = (px > 0.0) ? Math.abs(dx - rx - collisionRadius) : Math.abs(-rx + collisionRadius);
            if(del < delmin) {
                delmin = del;
                imin = i;
            }
        }
//        pAccumulator += 2*Math.abs(p.component(imin));
        p.setComponent(imin,-p.component(imin)); //multiply momentum component by -1
    }//end of bump

    /**
     * not yet implemented
     */
    public double lastCollisionVirial() {return Double.NaN;}
    
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

    
    public static void main(String[] args) {
        
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        sim.species.setNMolecules(10);
        P1HardBoundary potential = new P1HardBoundary(sim);
        potential.setCollisionRadius(0.5*Default.ATOM_SIZE);
  //      sim.phase.setBoundary(sim.space().makeBoundary(Space2D.Boundary.NONE));
        sim.elementCoordinator.go();
        Potential1.Agent potentialAgent = (Potential1.Agent)potential.getAgent(sim.phase);
        potentialAgent.setIterator(sim.phase.iteratorFactory().makeAtomIterator());
        
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,500);
        f.getContentPane().add(sim.panel());
        f.pack();
        f.show();
        f.addWindowListener(Simulation.WINDOW_CLOSER);
    }//end of main

}//end of P1HardBoundary
   
