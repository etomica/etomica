package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.DisplayTextBoxesCAE;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMd3D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public MeterPotentialEnergy energy;
    public AccumulatorAverageCollapsing avgEnergy;
    public DataPump pump;


    public LjMd3D() {
        super(Space3D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster(space);
        double sigma = 1.0;
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integrator.setTimeStep(0.02);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getMoleculeType().getChildTypes()[0]).setDiameter(sigma);
        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, 50);
        potential = new P2LennardJones(space, sigma, 1.0);
        AtomTypeLeaf leafType = species.getLeafType();

        potentialMaster.addPotential(potential,new AtomType[]{leafType,leafType});
        
        integrator.setBox(box);
        BoxImposePbc imposepbc = new BoxImposePbc(space);
        imposepbc.setBox(box);
        integrator.addIntervalAction(imposepbc);
		
        ConfigurationLattice configuration = new ConfigurationLattice(new LatticeCubicFcc(), space);
        configuration.initializeCoordinates(box);
        energy = new MeterPotentialEnergy(potentialMaster);
        energy.setBox(box);
        avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        pump = new DataPump(energy, avgEnergy);
        integrator.addIntervalAction(pump);
        integrator.setActionInterval(pump, 10);
    }

    public static class Applet extends javax.swing.JApplet {

        public void init() {
            final String APP_NAME = "LjMd3D";
            LjMd3D sim= new LjMd3D();
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, APP_NAME, 3, sim.space);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            simGraphic.getController().getDataStreamPumps().add(sim.pump);

            DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
            display.setAccumulator(sim.avgEnergy);
            simGraphic.add(display);
            getContentPane().add(simGraphic.getPanel());
        }
    }

    public static void main(String[] args) {
    	final String APP_NAME = "LjMd3D";
    	final LjMd3D sim = new LjMd3D();
    	final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 3, sim.space);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.getController().getDataStreamPumps().add(sim.pump);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        DisplayTextBoxesCAE display = new DisplayTextBoxesCAE();
        display.setAccumulator(sim.avgEnergy);
        simGraphic.add(display);
    }

}