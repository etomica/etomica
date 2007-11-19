package etomica.modules.vle;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2LJQ;
import etomica.potential.P2SoftTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Debye;
import etomica.units.Kelvin;

public class VLESim extends Simulation {

    public final Box boxLiquid, boxVapor;
    public final SpeciesSpheresRotating species;
    public final IntegratorMC integratorLiquid, integratorVapor;
    public final IntegratorManagerMC integratorGEMC;
    public final ActivityIntegrate activityIntegrate;
    protected double sigma;
    protected double temperature;
    protected double epsilon;
    protected double moment;
    protected double density;
    protected final P2LJQ p2LJQ;
    protected final P2SoftTruncated p2Truncated;
    
    public VLESim() {
        super(Space3D.getInstance());
        int initNumMolecules = 200;
        sigma = 3;
        temperature = Kelvin.UNIT.toSim(250);
        epsilon = Kelvin.UNIT.toSim(150);
        moment = Debye.UNIT.toSim(5);
        moment *= moment;
        density = 0.004;

        double initBoxSize = Math.pow(initNumMolecules/density, (1.0/3.0));
        
        species = new SpeciesSpheresRotating(this);
        getSpeciesManager().addSpecies(species);
        ((AtomTypeSphere)species.getMoleculeType()).setDiameter(sigma);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize));
        addBox(boxLiquid);
        boxVapor = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize));
        addBox(boxVapor);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc());
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);
        
        PotentialMaster potentialMaster = new PotentialMaster(space);
        p2LJQ = new P2LJQ(space, sigma, epsilon, moment);
        p2LJQ.setTemperature(temperature);
        p2Truncated = new P2SoftTruncated(p2LJQ, 4.0*sigma);
//        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getMoleculeType(), species.getMoleculeType()});
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        MCMoveRotate rotateMove = new MCMoveRotate(potentialMaster, random);
        integratorLiquid.getMoveManager().addMCMove(rotateMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        BoxImposePbc pbc = new BoxImposePbc(boxLiquid);
        integratorLiquid.addIntervalAction(pbc);
        integratorLiquid.setActionInterval(pbc, 100);
        
        integratorVapor = new IntegratorMC(potentialMaster, random, temperature);
        integratorVapor.setBox(boxVapor);
        atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
        rotateMove = new MCMoveRotate(potentialMaster, random);
        integratorVapor.getMoveManager().addMCMove(rotateMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        pbc = new BoxImposePbc(boxVapor);
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
    
    public void setSigma(double newSigma) {
        sigma = newSigma;
        p2LJQ.setSigma(sigma);
        p2Truncated.setTruncationRadius(4.0*sigma);
        ((AtomTypeSphere)species.getMoleculeType()).setDiameter(sigma);
    }

    public double getSigma() {
        return sigma;
    }
    
    public void setEpsilon(double newEpsilon) {
        p2LJQ.setEpsilon(newEpsilon);
    }
    
    public double getEpsilon() {
        return p2LJQ.getEpsilon();
    }
    
    public void setMoment(double newQ) {
        p2LJQ.setQuadrupolarMomentSquare(newQ*newQ);
    }
    
    public double getMoment() {
        return Math.sqrt(p2LJQ.getQuadrupolarMomentSquare());
    }
}
