package etomica.tests;

import etomica.Default;
import etomica.Phase;
import etomica.Simulation;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.config.ConfigurationFile;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionMolecular;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterNbr;
import etomica.potential.P1BondedHardSpheres;
import etomica.potential.P2HardBond;
import etomica.potential.P2SquareWell;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheres;

/**
 * Simple square-well chain simulation.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
 
public class TestSWChain extends Simulation {
    
    public IntegratorHard integrator;
    public Phase phase;

    public TestSWChain(Space space, int numMolecules) {
        super(space, true, new PotentialMasterNbr(space));
        int chainLength = 10;
        int numAtoms = numMolecules * chainLength;
        double sqwLambda = 1.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        Default.makeLJDefaults();
        double timeStep = 0.005;
        double simTime = 100000.0/numAtoms;
        int nSteps = (int)(simTime / timeStep);

        // makes eta = 0.35
        Default.BOX_SIZE = 14.4094*Math.pow((numAtoms/2000.0),1.0/3.0);
        integrator = new IntegratorHard(potentialMaster);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        NeighborListManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        integrator.addListener(nbrManager);
        nbrManager.setRange(Default.ATOM_SIZE*sqwLambda*neighborRangeFac);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(nSteps);
        int nCells = (int)(2*Default.BOX_SIZE/(neighborRangeFac*sqwLambda*Default.ATOM_SIZE));
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);

        P2SquareWell potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,0.5*Default.POTENTIAL_WELL);
        NeighborCriterion nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());

        SpeciesSpheres species = new SpeciesSpheres(this,chainLength);
        species.setNMolecules(numMolecules);
        P1BondedHardSpheres potentialChainIntra = new P1BondedHardSpheres(space);
        ((P2HardBond)potentialChainIntra.bonded).setBondLength(Default.ATOM_SIZE);
        ((P2HardBond)potentialChainIntra.bonded).setBondDelta(bondFactor);
        CriterionBondedSimple criterion = new CriterionBondedSimple(nbrCriterion);
        criterion.setBonded(false);
        potential.setCriterion(criterion);
        potentialChainIntra.setNonbonded(potential);
        criterion = new CriterionBondedSimple(NeighborCriterion.ALL);
        criterion.setBonded(true);
        potentialChainIntra.bonded.setCriterion(criterion);
        potentialMaster.setSpecies(potentialChainIntra, new Species[] {species});
        ((ConformationLinear)species.getFactory().getConformation()).setBondLength(Default.ATOM_SIZE);

        potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,0.5*Default.POTENTIAL_WELL);
        nbrCriterion = new CriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionMolecular criterionMolecular = new CriterionMolecular(nbrCriterion);
        criterionMolecular.setIntraMolecular(false);
        potential.setCriterion(criterionMolecular);
        
        AtomTypeSphere sphereType = (AtomTypeSphere)((AtomFactoryHomo)species.moleculeFactory()).childFactory().getType();
        potentialMaster.addPotential(potential,new AtomType[]{sphereType,sphereType});
        ((PotentialMasterNbr)potentialMaster).getNeighborManager().addCriterion(nbrCriterion,new AtomType[]{sphereType});

        phase = new Phase(this);

        integrator.addPhase(phase);
        new ConfigurationFile(space,"SWChain"+Integer.toString(numMolecules)).initializeCoordinates(phase);
    }
    
    public static void main(String[] args) {
        int numMolecules = 500;
        if (args.length > 0) {
            numMolecules = Integer.valueOf(args[0]).intValue();
        }
        TestSWChain sim = new TestSWChain(Space3D.getInstance(), numMolecules);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space,sim.integrator); 
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        energyMeter.setPhase(sim.phase);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage();
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        pMeter.setPhase(sim.phase);
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        double avgPE = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(AccumulatorAverage.AVERAGE.index)).x;
        avgPE /= numMolecules;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(AccumulatorAverage.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numMolecules;
        System.out.println("Cv/k="+Cv);
        
        if (Double.isNaN(Z) || Math.abs(Z-4.5) > 1.5) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+19.32) > 0.12) {
            System.exit(1);
        }
        // actual value ~2
        if (Double.isNaN(Cv) || Cv < 0.5 || Cv > 4.5) {
            System.exit(1);
        }
    }
}