package etomica;

import java.util.Random;
import etomica.units.Dimension;
import java.awt.*;

/** 
 * Monte Carlo move that changes the system volume by distorting the space, causing
 * more volume to be placed at one point selected at random.  Space is distorted from
 * this point until reaching the boundaries, where it becomes a uniform translation.
 *
 * @author Nandou Lu
 * @author David Kofke
 */

public class MCMovePointVolume extends MCMove {
    
    private final Random rand = new Random();
    protected double pressure;
    private ActionPointVolume action = new ActionPointVolume();
    private Space.Vector r0;

    public MCMovePointVolume() {
        super();
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.50);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        action.setPhase(phase);
    }
        
    public void thisTrial() {
        double hOld, hNew, vOld, vNew;
        vOld = phase.volume();
        hOld = phase.energy.potential() + pressure*vOld;
        r0 = phase.randomPosition();
        double step = rand.nextDouble()*stepSizeMax;
        if(Math.random() < 0.5) { //transform square to distorted
            action.actionPerformed(phase, true, step, r0);
        }
        else { //transform distorted back to square
            action.actionPerformed(phase, false, step, r0);
        }
        vNew = phase.volume();
        hNew = phase.energy.potential() + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||  //not correct yet
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature + action.lastLnJacobian())
                < Math.random()) 
            {  //reject
              action.undo();
            }
        nAccept++;   //accept
    }
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
        
 /*   public static void main(String args[]) {
        
        javax.swing.JFrame f = new javax.swing.JFrame();   //create a window
        f.setSize(600,350);
        
        Default.ATOM_SIZE = 1.0;
        
        Simulation sim = new Simulation(new Space2D());
        sim.setUnitSystem(new etomica.units.UnitSystem.LJ());
        Simulation.instance = sim;
        Species species = new SpeciesDisks(sim);
        Potential2 potential = new P2SimpleWrapper(sim,new PotentialHardDisk(sim));
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom();
        MCMovePointVolume mcMovePointVolume = new MCMovePointVolume();
        MCMoveVolume mcMoveVolume = new MCMoveVolume();
        Controller controller = new Controller(sim);
        final Phase phase = new Phase(sim);
        final DisplayPhase displayPhase = new DisplayPhase(sim);
        displayPhase.setScale(0.7);
 //       DeviceSlider slider = new DeviceSlider(mcMovePointVolume, "pressure");
        DeviceSlider slider = new DeviceSlider(mcMoveVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);
        
        Meter dMeter = new MeterDensity();
        DisplayBox box = new DisplayBox();
        box.setMeter(dMeter);
        box.setWhichValue(MeterAbstract.MOST_RECENT);
        
        mcMovePointVolume.setPhase(phase);
        integrator.add(mcMoveAtom);
   //     integrator.add(mcMovePointVolume);
        integrator.add(mcMoveVolume);

        species.setNMolecules(40);
        
        DisplayPlot plot = new DisplayPlot();
//        plot.setDataSource(dMeter);
        dMeter.setHistorying(true);
                
		Simulation.instance.elementCoordinator.go();
		
        f.getContentPane().add(Simulation.instance.panel()); //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }*/
}