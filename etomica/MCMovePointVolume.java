package etomica;

import java.util.Random;
import etomica.units.Dimension;
import etomica.lattice.*;

public class MCMovePointVolume extends MCMove {
    
    private final Random rand = new Random();
    protected double pressure;
    protected PhaseAction.Inflate inflate;
    private MyLattice lattice;

    public MCMovePointVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        inflate = new PhaseAction.Inflate(phase);
    }
    
    public void thisTrial() {
        double hOld, hNew, vOld, vNew;
        vOld = phase.volume();
        hOld = phase.energy.potential() + pressure*vOld;
        double vScale = (2.*Math.random()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/(double)phase.parentSimulation().space().D());
        inflate.actionPerformed(phase,rScale);
        hNew = phase.energy.potential() + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature+(phase.moleculeCount+1)*vScale)
                < Math.random()) 
            {  //reject
              inflate.retractAction();
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
    
    public static void main(String args[]) {
        
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);
        
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        Species species = new SpeciesDisks(sim);
        Potential2 potential = new P2SimpleWrapper(sim,new PotentialHardDisk(sim));
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom();
        MCMove mcMoveVolume = new MCMoveVolume();
        Controller controller = new Controller(sim);
        Phase phase = new Phase(sim);
        Display displayPhase = new DisplayPhase(sim);
        DeviceSlider slider = new DeviceSlider(mcMoveVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);
        
        mcMoveVolume.setPhase(phase);
        integrator.add(mcMoveAtom);
        integrator.add(mcMoveVolume);

		Simulation.instance.elementCoordinator.go();
		                                    
        f.add(Simulation.instance.panel()); //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
    
    /**
     * Array of 2D triangular cells formed by bisecting the cells of a square lattice.
     */
    private class MyLattice implements AbstractLattice {
        
        SquareLattice squareLattice;
        SiteIterator.List iterator = new SiteIterator.List();
        /**
         * Size is the number of sites in each dimension; considered non-periodically
         */
        MyLattice(int size) {
            squareLattice = new SquareLattice(size, new MySiteFactory(), 1.0/(double)(size-1));
            int[] index = new int[2];
            MySite[][] sites = new MySite[size][size];
            for(index[0]=0; index[0]<size; index[0]++) { //make array of sites for convenience of access
                for(index[1]=0; index[1]<size; index[1]++) {
                    sites[index[0]][index[1]] = (MySite)squareLattice.site(index);
                }
            }
            //construct two triangular cells for each site
            for(int i=0; i<size-1; i++) {    //loops don't do last row/column of lattice
                for(int j=0; j<size-1; j++) {
                    MySite site = sites[i][j];
                    site.c1 = new MyCell(site, sites[i][j+1], sites[i+1][j+1]);
                    site.c2 = new MyCell(site, sites[i+1][j], sites[i+1][j+1]);
                    deform(site);
                } 
            }
        }//end of MyLattice constructor
        
        public MyCell getDeformedCell(Space.Vector r) {
            //find and return cell containing r in deformed lattice
            return null;
        }
        public MyCell getOriginalCell(Space.Vector r) {
            //find and return cell containing r in undeformed lattice
            return null;
        }
        
        private void deform(MySite site) {
            //deformation algorithm
        }
        
        //AbstractLattice methods
        public int D() {return 2;}
        public int siteCount() {return 2*squareLattice.siteCount();}
        public Site site(AbstractLattice.Coordinate coord) {return null;} //not used
        public Site randomSite() {return iterator.randomSite();}
        public SiteIterator iterator() {return iterator;}     //iterator for all sites in lattice
    }
    
    private class MySite extends Site {
        MyCell c1, c2;
        Space.Vector originalPosition, deformedPosition;
        MySite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
            super(parent, iterator, coord);
            originalPosition = ((BravaisLattice.Coordinate)coordinate()).position();
            deformedPosition = new Space2D.Vector();
            deformedPosition.E(originalPosition);
        }
    }
    private class MySiteFactory implements SiteFactory {
       public Site makeSite(AbstractLattice parent, 
                            SiteIterator.Neighbor iterator, 
                            AbstractLattice.Coordinate coord) {
           return new MySite(parent, iterator, coord);
       }
    }
    private class MyCell {
        public MySite[] vertex = new MySite[3];
        MyCell(MySite s0, MySite s1, MySite s2) {
            vertex[0] = s0;
            vertex[1] = s1;
            vertex[2] = s2;
        }
    }
}