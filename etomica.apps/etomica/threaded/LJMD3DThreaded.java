package etomica.threaded;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureTensorFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Simple Lennard-Jones Monte Carlo simulation in 3D.
 * Initial configurations at http://rheneas.eng.buffalo.edu/etomica/tests/
 */
 
public class LJMD3DThreaded extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesSpheresMono species;
    public Phase phase;
    public P2LennardJones p2lj;
    public P2SoftSphericalTruncated[] potential;
    public PotentialThreaded potentialThreaded;
    public Controller controller;
    public ActivityIntegrate activityIntegrate;
   
    
    public LJMD3DThreaded() {
        this(500, 1);
    }

    public LJMD3DThreaded(int numAtoms, int numThreads) {
        super(Space3D.getInstance(), true, new PotentialMasterListThreaded(Space3D.getInstance()));
        defaults.makeLJDefaults();
        // need optimization of fac and time step
        double neighborFac = 1.35;

        // THREAD ZONE
        integrator = new IntegratorVelocityVerletThreaded(this, numThreads);
        integrator.setTemperature(1.0);
        integrator.setIsothermal(true);
        integrator.setTimeStep(0.001);
        activityIntegrate = new ActivityIntegrate(this,integrator);
        //activityIntegrate.setMaxSteps(500000);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
        phase = new Phase(this);
        species.getAgent(phase).setNMolecules(numAtoms);
        phase.setDensity(0.65);
        integrator.addListener(((PotentialMasterList)potentialMaster).getNeighborManager(phase));
        
       
        p2lj = new P2LennardJones(this);
        
        double truncationRadius = 2.5*p2lj.getSigma();
        if(truncationRadius > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
        }
        
        
        // -----------THREAD ZONE--------------\\
        ((PotentialMasterListThreaded)potentialMaster).setNumThreads(numThreads);
     
        potential = new P2SoftSphericalTruncated[numThreads];

        for(int i=0; i<numThreads; i++){
            potential[i] = new P2SoftSphericalTruncated(p2lj, truncationRadius);
        }
        
        potentialThreaded = new PotentialThreaded(space, potential);
        
        ((PotentialMasterListThreaded)potentialMaster).setCellRange(1);
        ((PotentialMasterListThreaded)potentialMaster).setRange(neighborFac * truncationRadius);
        ((PotentialMasterListThreaded)potentialMaster).getNeighborManager(phase).setQuiet(true);
        potentialMaster.addPotential(potentialThreaded, new Species[] {species, species});
        //--------------------------------------\\
        
//        new ConfigurationFile(space,"LJMC3D"+Integer.toString(numAtoms)).initializeCoordinates(phase);
        new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(phase);
        integrator.setPhase(phase);
//        WriteConfiguration writeConfig = new WriteConfiguration("LJMC3D"+Integer.toString(numAtoms),phase,1);
//        integrator.addListener(writeConfig);
    }

       
    public static void main(String[] args) {
        int numAtoms = 1000;
        int numThreads = 1;
        
        if (args.length > 0) {
            numAtoms = Integer.valueOf(args[0]).intValue();
        }
        LJMD3DThreaded sim = new LJMD3DThreaded(numAtoms, numThreads);
        sim.getController().actionPerformed();
        SimulationGraphic simgraphic = new SimulationGraphic(sim);
        simgraphic.makeAndDisplayFrame();
        /**
        sim.getDefaults().blockSize = 10;
        MeterPressureTensorFromIntegrator pMeter = new MeterPressureTensorFromIntegrator();
        pMeter.setIntegrator(sim.integrator);
        AccumulatorAverage pAccumulator = new AccumulatorAverage(sim);
        DataPump pPump = new DataPump(pMeter,pAccumulator);
        IntervalActionAdapter iaa = new IntervalActionAdapter(pPump,sim.integrator);
        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        energyMeter.setIncludeLrc(false);
        energyMeter.setPhase(sim.phase);
        AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        iaa = new IntervalActionAdapter(energyManager, sim.integrator);
        iaa.setActionInterval(50);

        MeterKineticEnergy KEMeter = new MeterKineticEnergy();
        KEMeter.setPhase(sim.phase);
        AccumulatorAverage KEAccumulator = new AccumulatorAverage(sim);
        DataPump KEPump = new DataPump(KEMeter, KEAccumulator);
        KEAccumulator.setBlockSize(50);
        IntervalActionAdapter iaaKE = new IntervalActionAdapter(KEPump, sim.integrator);
        iaaKE.setActionInterval(20);
        
        sim.getController().actionPerformed();
        
//        System.exit(0);
        
        double Z = ((DataTensor)((DataGroup)pAccumulator.getData()).getData(StatType.AVERAGE.index)).x.trace()/3*sim.phase.volume()/(sim.phase.moleculeCount()*sim.integrator.getTemperature());
        double avgPE = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.AVERAGE.index)).x;
        avgPE /= numAtoms;
        System.out.println("Z="+Z);
        System.out.println("PE/epsilon="+avgPE);
        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.STANDARD_DEVIATION.index)).x;
        Cv /= temp;
        Cv *= Cv/numAtoms;
        System.out.println("Cv/k="+Cv);
        double avgKE = ((DataDouble)((DataGroup)KEAccumulator.getData()).getData(StatType.AVERAGE.index)).x;
        avgKE /= numAtoms;
        System.out.println("KE="+avgKE);
        
        
        if (Double.isNaN(Z) || Math.abs(Z-0.11) > 0.15) {
            System.exit(1);
        }
        if (Double.isNaN(avgPE) || Math.abs(avgPE+4.355) > 0.04) {
            System.exit(1);
        }
        if (Double.isNaN(Cv) || Math.abs(Cv-0.6) > 0.4) {
            System.exit(1);
        }
        **/
    }

}
