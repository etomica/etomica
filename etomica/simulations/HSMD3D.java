package etomica.simulations;

import etomica.AtomType;
import etomica.ConfigurationLattice;
import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.activity.ActivityIntegrate;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.potential.P2HardSphere;
import etomica.space3d.Space3D;

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
        this(Space3D.getInstance());
    }
    
    //we use a second, private constructor to permit the space to
    //appear twice in the call to the superclass constructor; alternatively
    //we could have passed Space3D.getInstance() twice
    private HSMD3D(Space space) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(space, true, new PotentialMasterNbr(space));

        int numAtoms = 256;
        double neighborRangeFac = 1.6;
        Default.makeLJDefaults();
        Default.ATOM_SIZE = 1.0;
        Default.BOX_SIZE = 14.4573*Math.pow((numAtoms/2020.0),1.0/3.0);
        int nCells = (int)(Default.BOX_SIZE/neighborRangeFac);
        System.out.println("nCells: "+nCells);
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);

        integrator = new IntegratorHard(potentialMaster);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);
        this.register(integrator);

        NeighborListManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        nbrManager.setRange(Default.ATOM_SIZE*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(numAtoms);
//        Crystal crystal = new LatticeCubicFcc(space);
//        ConfigurationLattice configuration = new ConfigurationLattice(space, crystal);
//        phase.setConfiguration(configuration);
        potential = new P2HardSphere(space);
//        this.potentialMaster.setSpecies(potential,new Species[]{species,species});

        NeighborCriterion criterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        potential.setCriterion(criterion);
        potentialMaster.setSpecies(potential,new Species[]{species,species});

        nbrManager.addCriterion(criterion,new AtomType[]{species.getFactory().getType()});

        phase = new Phase(this);
        new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(phase);
        integrator.addPhase(phase);
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
        Default.DO_SLEEP = false;
        Default.FIX_OVERLAP = true;
        etomica.simulations.HSMD3D sim = new etomica.simulations.HSMD3D();
        SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DeviceNSelector nSelector = new DeviceNSelector(sim,sim.phase.getAgent(sim.species));
        simGraphic.add(nSelector);
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType.setColor(sim.species.getFactory().getType(), java.awt.Color.red);
        simGraphic.panel().setBackground(java.awt.Color.yellow);
    }//end of main

}//end of class
