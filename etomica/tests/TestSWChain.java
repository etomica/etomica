package etomica.tests;

import etomica.ConfigurationFile;
import etomica.ConformationLinear;
import etomica.DataSink;
import etomica.DataSource;
import etomica.Default;
import etomica.IntegratorPotentialEnergy;
import etomica.Phase;
import etomica.Simulation;
import etomica.Space;
import etomica.Species;
import etomica.SpeciesSpheres;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.iterator.ApiIntergroup;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceCOM;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionAll;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbratom.CriterionBondedSimple;
import etomica.nbratom.CriterionMolecular;
import etomica.nbratom.NeighborManager;
import etomica.nbratom.PotentialMasterNbr;
import etomica.potential.P1BondedHardSpheres;
import etomica.potential.P2HardBond;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialGroup;
import etomica.space3d.Space3D;

/**
 * Simple square-well chain simulation.
 * Initial configurations at http://gordon.eng.buffalo.edu/etomica/tests/
 */
 
public class TestSWChain extends Simulation {
    
    public IntegratorHard integrator;
    public Phase phase;

    public TestSWChain(Space space, int numMolecules) {
        super(space, new PotentialMasterNbr(space));
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
        NeighborManager nbrManager = ((PotentialMasterNbr)potentialMaster).getNeighborManager();
        integrator.addIntervalListener(nbrManager);
        nbrManager.setRange(Default.ATOM_SIZE*sqwLambda*neighborRangeFac);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(nSteps);
        int nCells = (int)(2*Default.BOX_SIZE/(neighborRangeFac*sqwLambda*Default.ATOM_SIZE));
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);
        ((PotentialMasterNbr)potentialMaster).setAtomPositionDefinition(new DataSourceCOM(space));

        P2SquareWell potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,0.5*Default.POTENTIAL_WELL);
        NeighborCriterion nbrCriterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());

        SpeciesSpheres species = new SpeciesSpheres(this,potentialMaster.sequencerFactory(),chainLength);
        species.setNMolecules(numMolecules);
        P1BondedHardSpheres potentialChainIntra = new P1BondedHardSpheres(space);
        ((P2HardBond)potentialChainIntra.bonded).setBondLength(Default.ATOM_SIZE);
        ((P2HardBond)potentialChainIntra.bonded).setBondDelta(bondFactor);
        CriterionBondedSimple criterion = new CriterionBondedSimple(nbrCriterion);
        criterion.setBonded(false);
        potential.setCriterion(criterion);
        potentialChainIntra.setNonbonded(potential);
        criterion = new CriterionBondedSimple(new NeighborCriterionAll());
        criterion.setBonded(true);
        potentialChainIntra.bonded.setCriterion(criterion);
        potentialMaster.setSpecies(potentialChainIntra, new Species[] {species});
        ((ConformationLinear)species.getFactory().getConformation()).setBondLength(Default.ATOM_SIZE);

        
        potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,0.5*Default.POTENTIAL_WELL);
        nbrCriterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        CriterionMolecular criterionMolecular = new CriterionMolecular(nbrCriterion);
        criterionMolecular.setIntraMolecular(false);
        potential.setCriterion(criterionMolecular);
        
        PotentialGroup p2Inter = new PotentialGroup(2,space);
        p2Inter.addPotential(potential,new ApiIntergroup());
        ((PotentialMasterNbr)potentialMaster).setSpecies(p2Inter,new Species[]{species,species});
        ((PotentialMasterNbr)potentialMaster).getNeighborManager().addCriterion(criterion);
        ((AtomFactoryHomo)species.moleculeFactory()).childFactory().getType().getNbrManagerAgent().addCriterion(criterionMolecular);

        phase = new Phase(this);

        integrator.addPhase(phase);
        phase.setConfiguration(new ConfigurationFile(space,"SWChain"+Integer.toString(numMolecules)));
    }
    
    public static void main(String[] args) {
        int numMolecules = 500;
        if (args.length > 0) {
            numMolecules = Integer.valueOf(args[0]).intValue();
        }
        TestSWChain sim = new TestSWChain(new Space3D(), numMolecules);

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator); 
        DataSource energyMeter = new IntegratorPotentialEnergy(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage();
        DataPump energyManager = new DataPump(energyMeter,new DataSink[]{energyAccumulator});
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar(sim.phase)*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.temperature());
        double[] data = energyAccumulator.getData();
        double PE = data[AccumulatorAverage.AVERAGE.index]/numMolecules;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+PE);
        double temp = sim.integrator.temperature();
        double Cv = data[AccumulatorAverage.STANDARD_DEVIATION.index]/(temp*temp*numMolecules);
        System.out.println("Cv/k="+Cv);
        
        if (Math.abs(Z-4.5) > 1.5) {
            System.exit(1);
        }
        if (Math.abs(PE+19.32) > 0.12) {
            System.exit(1);
        }
        if (Math.abs(Cv-0.05) > 0.05) {
            System.exit(1);
        }
    }
}