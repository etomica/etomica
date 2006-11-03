package etomica.simulation.prototypes;

import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;

/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class HSMD3D extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    /**
     * The Phase holding the atoms. 
     */
    public final Phase phase;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorHard integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    public final P2HardSphere potential;
    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public HSMD3D() {
        this(new Default());
    }
    
    public HSMD3D(Default defaults) {
        this(Space3D.getInstance(), defaults);
    }
    
    //we use a second, private constructor to permit the space to
    //appear twice in the call to the superclass constructor; alternatively
    //we could have passed Space3D.getInstance() twice
    private HSMD3D(Space space, Default defaults) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(space, true, new PotentialMasterList(space, 1.6), Default.BIT_LENGTH, defaults);

        int numAtoms = 256;
        double neighborRangeFac = 1.6;
        defaults.makeLJDefaults();
        defaults.atomSize = 1.0;
        defaults.boxSize = 14.4573*Math.pow((numAtoms/2020.0),1.0/3.0);
        ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*defaults.atomSize);

        integrator = new IntegratorHard(this);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);
        this.register(integrator);

        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager();
        nbrManager.setRange(defaults.atomSize*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        
        species = new SpeciesSpheresMono(this);
//        Crystal crystal = new LatticeCubicFcc(space);
//        ConfigurationLattice configuration = new ConfigurationLattice(space, crystal);
//        phase.setConfiguration(configuration);
        potential = new P2HardSphere(this);
//        this.potentialMaster.setSpecies(potential,new Species[]{species,species});

        potentialMaster.addPotential(potential,new Species[]{species,species});

        phase = new Phase(this);
        phase.getAgent(species).setNMolecules(numAtoms);
        new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(phase);
        integrator.setPhase(phase);
 //       integrator.addIntervalListener(new PhaseImposePbc(phase));
        
        //ColorSchemeByType.setColor(speciesSpheres0, java.awt.Color.blue);

 //       MeterPressureHard meterPressure = new MeterPressureHard(integrator);
 //       DataAccumulator accumulatorManager = new DataAccumulator(meterPressure);
        // 	DisplayBox box = new DisplayBox();
        // 	box.setDatumSource(meterPressure);
 //       phase.setDensity(0.7);
    } //end of constructor

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default defaults = new Default();
        defaults.doSleep = false;
        defaults.ignoreOverlap = true;
        etomica.simulation.prototypes.HSMD3D sim = new etomica.simulation.prototypes.HSMD3D(defaults);
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpeciesAgent(sim.phase.getAgent(sim.species));
        simGraphic.add(nSelector);
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
        simGraphic.panel().setBackground(java.awt.Color.yellow);
    }//end of main

}//end of class
