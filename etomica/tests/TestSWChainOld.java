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
import etomica.atom.iterator.AtomsetIteratorFiltered;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.DataSourceCOM;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionSimple;
import etomica.nbr.PotentialMasterNbr;
import etomica.potential.P1BondedHardSpheres;
import etomica.potential.P2HardBond;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialGroup;
import etomica.space3d.Space3D;

/**
 * Simple square-well chain simulation.
 */
 
public class TestSWChainOld extends Simulation {
    
    public IntegratorHard integrator;
    public Phase phase;

    public TestSWChainOld(Space space, int numMolecules) {
        super(space, new PotentialMasterNbr(space));
        int chainLength = 4;
        double sqwLambda = 1.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        Default.makeLJDefaults();
        double moleculeRange = Default.ATOM_SIZE*((chainLength-1)*(1+bondFactor) + sqwLambda*neighborRangeFac); 

        // makes eta = 0.35
        Default.BOX_SIZE = 14.4094*Math.pow((numMolecules*chainLength/2000.0),1.0/3.0);
        integrator = new IntegratorHard(potentialMaster);
        integrator.setTimeStep(0.006);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        integrator.addIntervalListener(((PotentialMasterNbr)potentialMaster).getNeighborManager());
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(500000/numMolecules);
        int nCells = (int)(2*Default.BOX_SIZE/moleculeRange);
        ((PotentialMasterNbr)potentialMaster).setNCells(nCells);
        ((PotentialMasterNbr)potentialMaster).setAtomPositionDefinition(new DataSourceCOM(space));

        P2SquareWell potential = new P2SquareWell(space,Default.ATOM_SIZE,sqwLambda,0.5*Default.POTENTIAL_WELL);

        SpeciesSpheres species = new SpeciesSpheres(this,potentialMaster.sequencerFactory(),chainLength);
        species.setNMolecules(numMolecules);
        P1BondedHardSpheres potentialChainIntra = new P1BondedHardSpheres(space);
        ((P2HardBond)potentialChainIntra.bonded).setBondLength(Default.ATOM_SIZE);
        ((P2HardBond)potentialChainIntra.bonded).setBondDelta(bondFactor);
        potentialChainIntra.setNonbonded(potential);
        potentialMaster.setSpecies(potentialChainIntra, new Species[] {species});
        ((ConformationLinear)species.getFactory().getConformation()).setBondLength(Default.ATOM_SIZE);

        PotentialGroup p2Inter = new PotentialGroup(2,space);
        NeighborCriterion criterion = new NeighborCriterionSimple(space,potential.getRange(),neighborRangeFac*potential.getRange());
        AtomsetIteratorFiltered interIterator = new AtomsetIteratorFiltered(new ApiIntergroup(),criterion);
        p2Inter.addPotential(potential,interIterator);
        ((PotentialMasterNbr)potentialMaster).setSpecies(p2Inter,new Species[]{species,species},moleculeRange);
        ((PotentialMasterNbr)potentialMaster).getNeighborManager().addCriterion(criterion);
        ((AtomFactoryHomo)species.moleculeFactory()).childFactory().getType().getNbrManagerAgent().addCriterion(criterion);

        phase = new Phase(this);

        phase.setConfiguration(null);
        integrator.addPhase(phase);
        phase.setConfiguration(new ConfigurationFile(space,"chain"+Integer.toString(numMolecules)));
    }
    
    public static void main(String[] args) {
        int numMolecules = 1000;
        if (args.length > 0) {
            numMolecules = Integer.valueOf(args[0]).intValue();
        }
        TestSWChainOld sim = new TestSWChainOld(new Space3D(), numMolecules);

        MeterPressureHard pMeter = new MeterPressureHard(sim.integrator); 
        DataSource energyMeter = new IntegratorPotentialEnergy(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage();
        DataPump energyManager = new DataPump(energyMeter,new DataSink[]{energyAccumulator});
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar(sim.phase)*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.temperature());
        double[] data = energyAccumulator.getData();
        double PE = data[1]/numMolecules;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+PE);
        double t2 = sim.integrator.temperature();
        t2 *= t2;
        double Cv = data[3]/t2/numMolecules;
        System.out.println("Cv/k="+Cv);
        
        if (Math.abs(Z-3.23) > 0.4) {
            System.exit(1);
        }
        if (Math.abs(PE+8.07) > 0.03) {
            System.exit(1);
        }
        if (Math.abs(Cv-0.026) > 0.01) {
            System.exit(1);
        }
    }
}