//includes main method
package etomica;

public class P1HardBoundary extends PotentialAbstract {

    private double collisionRadius = 0.0;
    
    P1HardBoundary(Simulation sim) {
        super(sim);
    }
    
    public PotentialAgent makeAgent(Phase p) {
        return new Agent(this, p);
    }
        
    public void setCollisionRadius(double d) {collisionRadius = d;}
    public double getCollisionRadius() {return collisionRadius;}
    public etomica.units.Dimension getCollisionRadiusDimension() {return etomica.units.Dimension.LENGTH;}

    //P1HardBoundary.Agent
    public class Agent extends PotentialAgent implements PotentialAgent.Hard {

        protected AtomIterator iterator;
        private Space.Vector dimensions;
        private double pAccumulator;
        /**
         * @param potential The parent potential making this agent
         * @param phase The phase in which this agent will be placed
         */
        public Agent(PotentialAbstract potential, Phase phase) {
            super(potential, phase);
            dimensions = parentPhase.boundary().dimensions();
            phase.boundaryMonitor.addObserver( new java.util.Observer() {
                public void update(java.util.Observable obs, Object obj) {
                    dimensions = (Space2D.Vector)parentPhase.boundary().dimensions();}
            });
        }
        
        protected void makeIterator() {
            //use default iterator only if potential was not previously turned off
            //in the phase by setting its iterator to NULL
           // if(!(iterator == AtomIterator.NULL)) 
                iterator = parentPhase.iteratorFactory().makeAtomIteratorUp();
        }
            
        public void setIterator(AtomIterator iterator) {
            this.iterator = iterator;
        }
        public AtomIterator iterator() {return iterator;}
    
        public double energy(IteratorDirective id) {return 0.0;}

        public void findCollisions(IteratorDirective id, 
                                    final IntegratorHardAbstract.CollisionHandler collisionHandler) {
            //collisions with one-body hard potentials are considered "uplist" of
            //all atoms, so no collisions with potential arise when doing a downlist iteration
            if(id.direction() == IteratorDirective.DOWN) return;
            
            collisionHandler.setPotential(this);
            iterator.reset(id);
            while(iterator.hasNext()) {
                Atom atom = iterator.next();
                double time = collisionTime(atom);
               //if(time >= 0) etc.  (maybe this would be more efficient)
                if(time < Double.MAX_VALUE) collisionHandler.addCollision(atom, time);
            }//end while
        }//end of findCollisions
    
        public double collisionTime(Atom a) {
            Space.Vector r = a.coordinate().position();
            Space.Vector p = a.coordinate().momentum();
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
                
        public void bump(IntegratorHardAbstract.Agent agent) {
            Atom a = agent.atom();
            Space.Vector r = a.coordinate().position();
            Space.Vector p = a.coordinate().momentum();
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
            pAccumulator += 2*Math.abs(p.component(imin));
            p.setComponent(imin,-p.component(imin)); //multiply momentum component by -1
        }//end of bump
        
    }//end of Agent
    
    public static void main(String[] args) {
        
        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
        Simulation.instance = sim;
        sim.species.setNMolecules(10);
        new P1HardBoundary(sim);
  //      sim.phase.setBoundary(sim.space().makeBoundary(Space2D.Boundary.NONE));
        sim.elementCoordinator.go();
        
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,500);
        f.getContentPane().add(sim.panel());
        f.pack();
        f.show();
        f.addWindowListener(Simulation.WINDOW_CLOSER);
    }//end of main

}//end of P1HardBoundary
   
