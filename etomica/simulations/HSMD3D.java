// Source file generated by Etomica

package etomica.simulations;

import etomica.DataManager;
import etomica.Default;
import etomica.IntegratorHard;
import etomica.MeterPressureHard;
import etomica.P2HardSphere;
import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheresMono;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;

public class HSMD3D extends Simulation {

    public Phase phase;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public P2HardSphere potential;
    
    public HSMD3D() {
        this(new etomica.Space3D());
    }
    private HSMD3D(Space space) {
        super(space, new PotentialMaster(space));
//        super(space, new PotentialMasterNbr(space));
        Default.makeLJDefaults();
        integrator = new IntegratorHard(potentialMaster);
//        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        integrator.setTimeStep(0.01);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(108);
        phase = new Phase(space);
//        Crystal crystal = new LatticeCubicFcc(space);
//        ConfigurationLattice configuration = new ConfigurationLattice(space, crystal);
//        phase.setConfiguration(configuration);
        potential = new P2HardSphere(space);
        this.potentialMaster.setSpecies(potential,new Species[]{species,species});
        
//      elementCoordinator.go();
        //explicit implementation of elementCoordinator activities
        phase.speciesMaster.addSpecies(species);
        integrator.addPhase(phase);
        integrator.addIntervalListener(new PhaseImposePbc(phase));
        
        //ColorSchemeByType.setColor(speciesSpheres0, java.awt.Color.blue);

        MeterPressureHard meterPressure = new MeterPressureHard();
        DataManager accumulatorManager = new DataManager(meterPressure);
        // 	DisplayBox box = new DisplayBox();
        // 	box.setDatumSource(meterPressure);
 //       phase.setDensity(0.5);
    } //end of constructor

}//end of class
