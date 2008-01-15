package etomica.modules.vle;

import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotential;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

public class VLESimLJ extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final Species species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    public final PotentialMasterCell potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final IPotential potential;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double density;
    
    public VLESimLJ() {
        super(Space3D.getInstance());
        int initNumMolecules = 2000;
        sigma = 1.0;
        temperature = 1.0;
        epsilon = 1.0;
        density = 0.15;
        double cutoff = 3*sigma;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize));
        addBox(boxLiquid);
        boxVapor = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize));
        addBox(boxVapor);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc());
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);
        
        potentialMaster = new PotentialMasterCell(this, cutoff);
        P2LennardJones p2LJ = new P2LennardJones(space, sigma, epsilon);
        potential = new P2SoftSphericalTruncated(p2LJ, cutoff);
        potentialMaster.addPotential(potential, new Species[]{species, species});
        potentialMaster.setCellRange(3);
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        integratorVapor = new IntegratorMC(potentialMaster, random, temperature);
        integratorVapor.setBox(boxVapor);
        atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        integratorGEMC = new IntegratorManagerMC(random);
        integratorGEMC.setGlobalMoveInterval(2);
        integratorGEMC.addIntegrator(integratorLiquid);
        integratorGEMC.addIntegrator(integratorVapor);
        MCMoveVolumeExchangeVLE volumeExchange = new MCMoveVolumeExchangeVLE(
                potentialMaster, random, integratorLiquid,integratorVapor);
        ((MCMoveStepTracker)volumeExchange.getTracker()).setNoisyAdjustment(true);
        MCMoveMoleculeExchangeVLE moleculeExchange = new MCMoveMoleculeExchangeVLE(
                potentialMaster, random, integratorLiquid,integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);
        
        activityIntegrate = new ActivityIntegrate(integratorGEMC);
        getController().addAction(activityIntegrate);

        integratorLiquid.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(boxLiquid).makeMCMoveListener());
        integratorVapor.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(boxVapor).makeMCMoveListener());
        potentialMaster.getNbrCellManager(boxLiquid).assignCellAll();
        potentialMaster.getNbrCellManager(boxVapor).assignCellAll();
    }
}
