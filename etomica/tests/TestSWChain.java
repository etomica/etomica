package etomica.tests;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeSphere;
import etomica.atom.iterator.ApiBuilder;
import etomica.config.ConfigurationFile;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressureHard;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntervalActionAdapter;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardBond;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
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

    public TestSWChain() {
        this(500);
    }
    
    public TestSWChain(int numMolecules) {
        super(Space3D.getInstance(), true, new PotentialMasterList(Space3D.getInstance()));
        int chainLength = 10;
        int numAtoms = numMolecules * chainLength;
        double sqwLambda = 1.5;
        double neighborRangeFac = 1.2;
        double bondFactor = 0.15;
        defaults.makeLJDefaults();
        double timeStep = 0.005;
        double simTime = 100000.0/numAtoms;
        int nSteps = (int)(simTime / timeStep);

        // makes eta = 0.35
        defaults.boxSize = 14.4094*Math.pow((numAtoms/2000.0),1.0/3.0);
        integrator = new IntegratorHard(this);
        integrator.setTimeStep(timeStep);
        integrator.setIsothermal(true);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager();
        integrator.addListener(nbrManager);
        nbrManager.setRange(defaults.atomSize*sqwLambda*neighborRangeFac);
        getController().addAction(activityIntegrate);
        activityIntegrate.setMaxSteps(nSteps);
        ((PotentialMasterList)potentialMaster).setCellRange(2);
        ((PotentialMasterList)potentialMaster).setRange(neighborRangeFac*sqwLambda*defaults.atomSize);

        SpeciesSpheres species = new SpeciesSpheres(this,chainLength);
        species.setNMolecules(numMolecules);
        P2HardBond bonded = new P2HardBond(this);
        PotentialGroup potentialChainIntra = potentialMaster.makePotentialGroup(1);
        potentialChainIntra.addPotential(bonded, ApiBuilder.makeAdjacentPairIterator());

        bonded.setBondLength(defaults.atomSize);
        bonded.setBondDelta(bondFactor);

        potentialMaster.addPotential(potentialChainIntra, new Species[] {species});
        ((ConformationLinear)species.getFactory().getConformation()).setBondLength(defaults.atomSize);

        P2SquareWell potential = new P2SquareWell(space,defaults.atomSize,sqwLambda,0.5*defaults.potentialWell,false);

        AtomTypeSphere sphereType = (AtomTypeSphere)((AtomFactoryHomo)species.moleculeFactory()).getChildFactory().getType();
        potentialMaster.addPotential(potential,new AtomType[]{sphereType,sphereType});
        CriterionInterMolecular sqwCriterion = (CriterionInterMolecular)((PotentialMasterList)potentialMaster).getCriterion(potential);
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        sqwCriterion.setIntraMolecularCriterion(nonBondedCriterion);

        phase = new Phase(this);

        integrator.setPhase(phase);
        new ConfigurationFile(space,"SWChain"+Integer.toString(numMolecules)).initializeCoordinates(phase);
    }
    
    public static void main(String[] args) {
        int numMolecules = 500;
        if (args.length > 0) {
            numMolecules = Integer.valueOf(args[0]).intValue();
        }
        TestSWChain sim = new TestSWChain(numMolecules);

        MeterPressureHard pMeter = new MeterPressureHard(sim.space);
        pMeter.setIntegrator(sim.integrator);
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        energyMeter.setPhase(sim.phase);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        new IntervalActionAdapter(energyManager, sim.integrator);
        
        sim.getController().actionPerformed();
        
        double Z = pMeter.getDataAsScalar()*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        double avgPE = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.AVERAGE.index)).x;
        avgPE /= numMolecules;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.STANDARD_DEVIATION.index)).x;
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