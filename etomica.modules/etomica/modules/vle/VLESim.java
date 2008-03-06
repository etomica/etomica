package etomica.modules.vle;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorManagerMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LJQ;
import etomica.potential.P2SoftMoleculeMonatomicTruncated;
import etomica.potential.P2SoftTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Debye;
import etomica.units.Kelvin;

public class VLESim extends Simulation {

    public final IBox boxLiquid, boxVapor;
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
        boolean doNBR = false;
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
        ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);

        boxLiquid = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxLiquid);
        boxVapor = new Box(new BoundaryRectangularPeriodic(space, random, initBoxSize), space);
        addBox(boxVapor);
        boxLiquid.setNMolecules(species, initNumMolecules);
        boxVapor.setNMolecules(species, initNumMolecules);
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(), space);
        config.initializeCoordinates(boxLiquid);
        config.initializeCoordinates(boxVapor);

        double range = 15.0;
        PotentialMaster potentialMaster = new PotentialMaster(space);
        if (doNBR) {
            potentialMaster = new PotentialMasterCell(this, range, space);
            ((PotentialMasterCell)potentialMaster).setCellRange(2);
        }
        p2LJQ = new P2LJQ(space, sigma, epsilon, moment);
        p2LJQ.setTemperature(temperature);
        p2Truncated = new P2SoftTruncated(p2LJQ, range);
//        ((P2SoftSphericalTruncatedBox)potential).setTruncationFactor(0.35);
        if (doNBR) {
            potentialMaster.addPotential(p2Truncated, new AtomType[]{species.getLeafType(), species.getLeafType()});
        }
        else {
            potentialMaster.addPotential(new P2SoftMoleculeMonatomicTruncated(space, p2Truncated), new ISpecies[]{species,species});
        }
        
        integratorLiquid = new IntegratorMC(potentialMaster, random, temperature);
        integratorLiquid.setBox(boxLiquid);
        MCMoveAtom atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorLiquid.getMoveManager().addMCMove(atomMove);
        MCMoveRotate rotateMove = new MCMoveRotate(potentialMaster, random);
        integratorLiquid.getMoveManager().addMCMove(rotateMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);
        
        integratorVapor = new IntegratorMC(potentialMaster, random, temperature);
        integratorVapor.setBox(boxVapor);
        atomMove = new MCMoveAtom(potentialMaster, random, 0.5, 5.0, true);
        integratorVapor.getMoveManager().addMCMove(atomMove);
        rotateMove = new MCMoveRotate(potentialMaster, random);
        integratorVapor.getMoveManager().addMCMove(rotateMove);
//        ((MCMoveStepTracker)atomMove.getTracker()).setNoisyAdjustment(true);

        if (!doNBR) {
            BoxImposePbc pbc = new BoxImposePbc(boxLiquid, space);
            integratorLiquid.addIntervalAction(pbc);
            integratorLiquid.setActionInterval(pbc, 100);
            pbc = new BoxImposePbc(boxVapor, space);
            integratorVapor.addIntervalAction(pbc);
            integratorVapor.setActionInterval(pbc, 100);
        }
        
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

        if (doNBR) {
            ((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxLiquid)).assignCellAll();
            ((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxVapor)).assignCellAll();
            integratorLiquid.getMoveEventManager().addListener(((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxLiquid)).makeMCMoveListener());
            integratorVapor.getMoveEventManager().addListener(((NeighborCellManager)((PotentialMasterCell)potentialMaster).getCellAgentManager().getAgent(boxVapor)).makeMCMoveListener());
        }
    }
    
    public void setSigma(double newSigma) {
        sigma = newSigma;
        p2LJQ.setSigma(sigma);
        p2Truncated.setTruncationRadius(4.0*sigma);
        ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);
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
