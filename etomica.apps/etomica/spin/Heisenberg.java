package etomica.simulation.prototypes;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPhaseSpin2D;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.site.PotentialMasterSite;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Space2D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.spin.ConfigurationAligned;
import etomica.spin.MCMoveSpinFlip;
import etomica.spin.MeterSpin;
import etomica.spin.NeighborCellManagerFixed;
import etomica.spin.P1MagneticField;
import etomica.spin.P2Spin;
import etomica.units.systems.LJ;


/**
 * Simulation of a simple 2D Ising model.  Prototype
 * for simulation of a more general magentic system.
 *
 * @author David Kofke
 *
 */

/*
 * History
 * Created on May 22, 2005 by kofke
 */
public class Heisenberg extends Simulation {

    public Heisenberg() {
        this(Space2D.getInstance());
    }
    
    /**
     * 
     */
    public Heisenberg(Space space) {
        super(space, false, new PotentialMasterSite(space));
        defaults.makeLJDefaults();
        phase = new Phase(this);
        int nCells = 60;
        int numAtoms = space.powerD(nCells);
        phase.setCellManager(new NeighborCellManagerFixed(phase,nCells));
        spins = new SpeciesSpheresMono(this);
        spins.setNMolecules(numAtoms);
        phase.makeMolecules();
        new ConfigurationAligned(space).initializeCoordinates(phase);
        
        potential = new P2Spin(space);
        field = new P1MagneticField(space);
        integrator = new IntegratorMC(this);
        mcmove = new MCMoveSpinFlip(potentialMaster);
        integrator.getMoveManager().addMCMove(mcmove);
        
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setDoSleep(false);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);

        AtomType type = spins.getFactory().getType();
        potentialMaster.addPotential(field, new AtomType[] {type});
        potentialMaster.addPotential(potential, new AtomType[] {type, type});
        
        integrator.setPhase(phase);
        ((PotentialMasterSite)potentialMaster).updateTypeList(phase);
        phase.getCellManager().assignCellAll();
        
        meter = new MeterSpin(space);
        meter.setPhase(phase);
        dAcc = new AccumulatorAverage(this);
        pump = new DataPump(meter, dAcc);
        adapter = new IntervalActionAdapter(pump,integrator);
        adapter.setActionInterval(10);

    }

    public Phase phase;
    public Species spins;
    public P2Spin potential;
    public P1MagneticField field;
    private IntegratorMC integrator;
    public MCMoveSpinFlip mcmove;
    public MeterSpin meter;
    public DataPump pump;
    public IntervalActionAdapter adapter;
    public AccumulatorAverage dAcc;
    
    public static void main(String[] args) {
        Heisenberg sim = new Heisenberg(Space2D.getInstance());
        sim.register(sim.integrator);
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DisplayPhase displayPhase = simGraphic.getDisplayPhase(sim.phase);
        simGraphic.remove(displayPhase);
        displayPhase.setPhaseCanvas(new DisplayPhaseSpin2D(displayPhase));
        simGraphic.add(displayPhase);
        DeviceSlider temperatureSlider = new DeviceSlider(sim.getController(), sim.integrator,"temperature");
        temperatureSlider.setMinimum(0.5);
        temperatureSlider.setMaximum(10.0);
        temperatureSlider.setShowBorder(true);
        LJ lj = new LJ();
        temperatureSlider.setUnit(lj.temperature());
        simGraphic.add(temperatureSlider);
        temperatureSlider.setValue(sim.integrator.getTemperature());
        DeviceSlider fieldSlider = new DeviceSlider(sim.getController(), sim.field, "h");
        fieldSlider.setMinimum(-5.);
        fieldSlider.setMaximum(+5.);
        fieldSlider.setNMajor(5);
        fieldSlider.setValue(0.0);
        fieldSlider.setShowBorder(true);
        fieldSlider.setLabel("Magnetic field");
        simGraphic.add(fieldSlider);
        
        DisplayBoxesCAE boxes = new DisplayBoxesCAE();
        boxes.setAccumulator(sim.dAcc);
        boxes.setLabel("Magnetization");
        boxes.setLabelType(DisplayBox.BORDER);
        simGraphic.add(boxes);

        simGraphic.makeAndDisplayFrame();
    }
}
