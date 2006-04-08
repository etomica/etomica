package etomica.integrator.mcmove;

import etomica.action.PhaseInflate;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.atom.iterator.AtomIteratorNull;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Pressure;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */

public class MCMoveVolume extends MCMoveStep {
    
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    protected final PhaseInflate inflate;
    private final int D;

    private transient double uOld, hOld, vNew, vScale;
    private transient double uNew = Double.NaN;

    public MCMoveVolume(Simulation sim) {
        this(sim.potentialMaster, sim.space, sim.getDefaults().pressure);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolume(PotentialMaster potentialMaster, Space space, double pressure) {
        super(potentialMaster, new MCMoveStepTracker(), 1);
        this.D = space.D();
        inflate = new PhaseInflate(potentialMaster.getSpace());
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        setName("MCMoveVolume");
        
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
        inflate.setPhase(p);
    }
    
    public boolean doTrial() {
        double vOld = phase.volume();
        uOld = energyMeter.getDataAsScalar();
        hOld = uOld + pressure*vOld;
        vScale = (2.*Simulation.random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        inflate.setScale(rScale);
        inflate.actionPerformed();
        uNew = Double.NaN;
        return true;
    }//end of doTrial
    
    public double getA() {
        return Math.exp((phase.moleculeCount()+1)*vScale);
    }
    
    public double getB() {
        uNew = energyMeter.getDataAsScalar();
        double hNew = uNew + pressure*vNew;
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
        inflate.undo();
    }

    public double energyChange(Phase p) {return (p == phase) ? uNew - uOld : 0.0;}
    
    public AtomIterator affectedAtoms(Phase p) {
        if(p != phase) return AtomIteratorNull.INSTANCE;
        return new AtomIteratorAllMolecules(phase);
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public final double pressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}
    public final void setLogPressure(int lp) {setPressure(Math.pow(10.,lp));}
    
    /**
     * main method to test and demonstrate this class
     */
/*    public static void main(String args[]) {
                    
        Simulation sim = new Simulation(new Space2D());
        Simulation.instance = sim;
        Species species = new SpeciesSpheresMono(sim);
//        species.setNMolecules(2);
        P2HardSphere potential = new P2HardSphere();
//        Potential potential = new P2LennardJones();
        IntegratorMC integrator = new IntegratorMC(sim);
        MCMove mcMoveAtom = new MCMoveAtom(integrator);
        MCMove mcMoveVolume = new MCMoveVolume(integrator);
        Controller controller = new Controller(sim);
        Phase phase = new Phase(sim);
        Display displayPhase = new DisplayPhase(sim);
        DeviceSlider slider = new DeviceSlider(mcMoveVolume, "pressure");
        slider.setMinimum(0);
        slider.setMaximum(100);

		Meter energy = new MeterPotentialEnergy();
		energy.setPhase(phase);
		energy.setHistorying(true);
		energy.setActive(true);		
		energy.getHistory().setNValues(500);		
		DisplayPlot plot = new DisplayPlot();
		plot.setLabel("Energy");
		plot.setDataSources(energy.getHistory());

	    Simulation.instance.elementCoordinator.go();
        
  //      phase.setDensity(0.1);
    		                                    
        Simulation.makeAndDisplayFrame(sim);
    }//end of main  
    */
}