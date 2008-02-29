package etomica.modules.vle;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.IPotential;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedBox;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

public class VLESim2 extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final ISpecies species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    public final PotentialMaster potentialMaster;
    public final ActivityIntegrate activityIntegrate;
    public final IPotential potential;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double density;
    
    public VLESim2() {
        super(Space3D.getInstance());
        int initNumMolecules = 500;
        sigma = 1.0;
        temperature = 0.8;
        epsilon = 1.0;
        density = 0.1;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxLiquid);
        boxVapor = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxVapor);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);
        
        potentialMaster = new PotentialMaster(space);
        P2LennardJones p2LJ = new P2LennardJones(space, sigma, epsilon);
        potential = new P2SoftSphericalTruncatedBox(p2LJ);
        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        BoxImposePbc pbc = new BoxImposePbc(boxLiquid, space);
        integratorLiquid.addIntervalAction(pbc);
        integratorLiquid.setActionInterval(pbc, 100);
        
        integratorVapor = new IntegratorMC(potentialMaster, random, temperature);
        integratorVapor.setBox(boxVapor);
        atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        pbc = new BoxImposePbc(boxVapor, space);
        integratorVapor.addIntervalAction(pbc);
        integratorVapor.setActionInterval(pbc, 100);
        
        integratorGEMC = new IntegratorManagerMC(random);
        integratorGEMC.setGlobalMoveInterval(2);
        integratorGEMC.addIntegrator(integratorLiquid);
        integratorGEMC.addIntegrator(integratorVapor);
        MCMoveVolumeExchangeVLE volumeExchange = new MCMoveVolumeExchangeVLE(
                potentialMaster, random, integratorLiquid,integratorVapor);
        volumeExchange.setStepSize(0.05);
//        ((MCMoveStepTracker)volumeExchange.getTracker()).setNoisyAdjustment(true);
        MCMoveMoleculeExchangeVLE moleculeExchange = new MCMoveMoleculeExchangeVLE(
                potentialMaster, random, integratorLiquid,integratorVapor);
        integratorGEMC.getMoveManager().addMCMove(volumeExchange);
        integratorGEMC.getMoveManager().addMCMove(moleculeExchange);
        integratorGEMC.getMoveManager().setFrequency(volumeExchange, 0.01);
        
        activityIntegrate = new ActivityIntegrate(integratorGEMC);
        getController().addAction(activityIntegrate);
    }
}
