package etomica.simulations;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Species;
import etomica.SpeciesPistonCylinder;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.DeviceSlider;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayTimer;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorHardPiston;
import etomica.potential.P1HardBoundary;
import etomica.potential.P1HardMovingBoundary;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2;
import etomica.space.BoundaryNone;
import etomica.space2d.Space2D;
import etomica.units.Bar;
import etomica.units.Kelvin;
import etomica.units.PrefixedUnit;
import etomica.units.Unit;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends Simulation {
    
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public Potential2 potential;
    public P1HardBoundary wallPotential;
    public P1HardMovingBoundary pistonPotential;
    public ActivityIntegrate ai;
    public double lambda;

    public PistonCylinder() {
        super(new Space2D());
        lambda = 1.5;
        
        Default.ATOM_SIZE=2.0;
	    species = new SpeciesSpheresMono(this);
	    species.setNMolecules(120);
        Default.BOX_SIZE=50;
	    phase = new Phase(space);
        phase.setBoundary(new BoundaryNone(space));
        phase.speciesMaster.addSpecies(species);
	    
	    potential = new P2SquareWell(space,Default.ATOM_SIZE,lambda,Default.POTENTIAL_WELL);
//        potential = new P2HardSphere(space,Default.ATOM_SIZE);
	    potentialMaster.setSpecies(potential,new Species[]{species,species});
	    
        wallPotential = new P1HardBoundary(space);
        wallPotential.setCollisionRadius(Default.ATOM_SIZE*0.5); //potential.getCoreDiameter()*0.5);
        potentialMaster.setSpecies(wallPotential,new Species[]{species});
        
        pistonPotential = new P1HardMovingBoundary(space,phase.boundary(),1,Default.ATOM_MASS*1000);
        pistonPotential.setCollisionRadius(Default.ATOM_SIZE*0.5);
        pistonPotential.setWallPosition(-Default.ATOM_SIZE);
        pistonPotential.setWallVelocity(0.5);
        pistonPotential.setPressure(200000.0);
        potentialMaster.setSpecies(pistonPotential,new Species[]{species});
	    
	    integrator = new IntegratorHardPiston(potentialMaster,pistonPotential);
        integrator.addPhase(phase);
        integrator.setIsothermal(true);
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
	    
    }

    public static class Applet extends javax.swing.JApplet {
        public DisplayPhase display;
        public DataSourceCountSteps meterCycles;
        public DisplayBox displayCycles;
        public MeterTemperature thermometer;

        public void init() {
            PistonCylinder pc = new PistonCylinder();
            pc.ai.setDoSleep(true);
            pc.ai.setSleepPeriod(20);

            SimulationGraphic sg = new SimulationGraphic(pc);
            getContentPane().add(sg.panel());
            display = new DisplayPhase(pc.phase);
            
            meterCycles = new DataSourceCountSteps(pc.integrator);
            displayCycles = new DisplayBox(meterCycles);
            pc.integrator.addIntervalListener(displayCycles);
            
            sg.panel().add(displayCycles.graphic());
            sg.panel().setBackground(java.awt.Color.yellow);

//            pc.phase.setBoundary(speciesPC.new Boundary(phase)); //have piston-cylinder system define boundary of phase
            
            //part unique to this class
            thermometer = new MeterTemperature();
            thermometer.setPhase(new Phase[]{pc.phase});
            DisplayBox tBox = new DisplayBox();
            pc.integrator.addIntervalListener(tBox);
            tBox.setDataSource(thermometer);
            tBox.setUnit(new PrefixedUnit(Kelvin.UNIT));
            sg.panel().add(tBox.graphic());
            display.setAlign(1,DisplayPhase.BOTTOM);
            
            DisplayTimer timer = new DisplayTimer(pc.integrator);
            pc.integrator.addIntervalListener(timer);
            timer.setUpdateInterval(10);
            sg.panel().add(timer.graphic());

/*            DeviceSlider pressureSlider = new DeviceSlider((SpeciesPistonCylinder.PistonPressureField)phase.firstField(),"pressure");
            sg.panel().add(pressureSlider.graphic());
            pressureSlider.setUnit(new PrefixedUnit(Bar.UNIT));
            pressureSlider.setMinimum(50);
            pressureSlider.setMaximum(1000);*/
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
