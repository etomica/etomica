package etomica.simulations;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2;
import etomica.space.BoundaryNone;
import etomica.space2d.Space2D;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public Potential2 potential;
    public P1HardBoundary wallPotential;
    public ActivityIntegrate ai;
    public double lambda;

    public PistonCylinder() {
        super(new Space2D());
        lambda = 1.5;
        
	    species = new SpeciesSpheresMono(this);
//	    speciesPC = new SpeciesPistonCylinder(this);
//        speciesPC.setLength(20.);
        
	    phase = new Phase(space);
//        phase.setConfiguration(new ConfigurationFile(space,"pc"));
        phase.setBoundary(new BoundaryNone(space));
        phase.speciesMaster.addSpecies(species);
	    
	    potential = new P2SquareWell(space,Default.ATOM_SIZE,lambda,Default.POTENTIAL_WELL);
//        potential = new P2HardSphere(space,Default.ATOM_SIZE);
	    potentialMaster.setSpecies(potential,new Species[]{species,species});
	    
        wallPotential = new P1HardBoundary(space);
        wallPotential.setCollisionRadius(Default.ATOM_SIZE*0.5); //potential.getCoreDiameter()*0.5);
        potentialMaster.setSpecies(wallPotential,new Species[]{species});
        
	    
//        PotentialGroup potentialDiskPC = new PotentialGroup(2);
//	    potentialDiskPC.setSpecies(new Species[] {species, speciesPC});
	    
//	    potentialHardDiskWall = new P2HardSphereWall(potentialDiskPC, Default.ATOM_SIZE);
//	    potentialHardDiskWall.setIterator(new ApiGeneral(Simulation.instance.space,
//	                new AtomIteratorSinglet(), new AtomIteratorSequential()));
	    
	    integrator = new IntegratorHard(potentialMaster);
        integrator.addPhase(phase);
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
	    
    }

    public static class Applet extends javax.swing.JApplet {
        public DisplayPhase display;
        public DataSourceCountSteps meterCycles;
//        public DisplayBox displayCycles;
        public MeterTemperature thermometer;

        public void init() {
            PistonCylinder pc = new PistonCylinder();
            Default.DO_SLEEP = true;
            SimulationGraphic sg = new SimulationGraphic(pc);
            getContentPane().add(sg.panel());
            display = new DisplayPhase(pc.phase);
            
            meterCycles = new DataSourceCountSteps(pc.integrator);
//            displayCycles = new DisplayBox(meterCycles);
            
            sg.panel().setBackground(java.awt.Color.yellow);

//            pc.phase.setBoundary(speciesPC.new Boundary(phase)); //have piston-cylinder system define boundary of phase
            
            //part unique to this class
            thermometer = new MeterTemperature();
//            DisplayBox tBox = new DisplayBox();
//            tBox.setMeter(thermometer);
//            tBox.setUnit(new PrefixedUnit(Kelvin.UNIT));
            display.setAlign(1,DisplayPhase.BOTTOM);
            
//            DisplayTimer timer = new DisplayTimer(pc.integrator);
//            timer.setUpdateInterval(10);
//          Simulation.instance.elementCoordinator.go();

    /*        DeviceSlider pressureSlider = new DeviceSlider((SpeciesPistonCylinder.PistonPressureField)phase.firstField(),"pressure");
            Simulation.instance.add(pressureSlider.graphic(null));
            pressureSlider.setUnit(new Unit(Bar.UNIT));
            pressureSlider.setMinimum(50);
            pressureSlider.setMaximum(1000);
            */
        }
}
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        PistonCylinder sim = new PistonCylinder();
        sim.ai.setMaxSteps(10000);
        sim.getController().actionPerformed();
//		sim.makeAndDisplayFrame();
   }//end of main
    
}
