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
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private AtomIterator affectedAtomIterator;

    public MCMovePointVolume(IntegratorMC parentIntegrator) {
        super(parentIntegrator);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.50);
        setPressure(Default.PRESSURE);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        action.setPhase(phase);
        affectedAtomIterator = phase.makeAtomIterator();
    }
        
    public boolean thisTrial() {
        double hOld, hNew, vOld, vNew;
        vOld = phase.volume();
        hOld = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum() + pressure*vOld;
        r0 = phase.randomPosition();
        double step = rand.nextDouble()*stepSizeMax;
        if(Math.random() < 0.5) { //transform square to distorted
            action.actionPerformed(phase, true, step, r0);
        }
        else { //transform distorted back to square
            action.actionPerformed(phase, false, step, r0);
        }
        vNew = phase.volume();
        hNew = potential.set(phase).calculate(iteratorDirective, energy.reset()).sum() + pressure*vNew;
        if(hNew >= Double.MAX_VALUE ||  //not correct yet
             Math.exp(-(hNew-hOld)/parentIntegrator.temperature + action.lastLnJacobian())
                < Math.random()) {  //reject
              action.undo();
              return false;
            }
        return true;   //accept
    }
    
    public AtomIterator affectedAtoms() {return affectedAtomIterator;}
    
    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Dimension.PRESSURE;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,(double)lp));}
        
    public static void main(String args[]) {
        
        Default.ATOM_SIZE = 1.0;
        
        Simulation sim = new Simulation(new Space2D());
        sim.setUnitSystem(new etomica.units.UnitSystem.LJ());
        Simulation.instance = sim;
        Species species = new SpeciesSpheres(sim);
        Potential2 potential = new P2HardSphere();
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom(integrator);
        MCMovePointVolume mcMovePointVolume = new MCMovePointVolume(integrator);
        MCMoveVolume mcMoveVolume = new MCMoveVolume(integrator);
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

        species.setNMolecules(40);
        
        DisplayPlot plot = new DisplayPlot();
//        plot.setDataSource(dMeter);
        dMeter.setHistorying(true);
                
		Simulation.instance.elementCoordinator.go();
		Simulation.makeAndDisplayFrame(Simulation.instance);
    }
}