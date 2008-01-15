package etomica.simulation.prototypes;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 2D
 */
 
public class LjMd2D extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Box box;
    public P2LennardJones potential;
    public Controller controller;
    public DisplayBox display;
    public DisplayPlot plot;
    public MeterEnergy energy;

    public LjMd2D() {
        super(Space2D.getInstance(), false);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        box = new Box(this);
        addBox(box);
        box.setNMolecules(species, 50);
        new ConfigurationLattice(new LatticeOrthorhombicHexagonal()).initializeCoordinates(box);
        potential = new P2LennardJones(space);
        potentialMaster.addPotential(potential,new AtomType[]{species.getLeafType(),species.getLeafType()});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        integrator.setBox(box);
		
		energy = new MeterEnergy(potentialMaster);
//		energy.setHistorying(true);
//		
//		energy.getHistory().setHistoryLength(500);
//		
//		plot = new DisplayPlot(this);
//		plot.setLabel("Energy");
//		plot.setDataSources(energy.getHistory());
		
    }
    
}