package etomica.yukawa;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.phase.Phase;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * A Yukawa Monte-Carlo simulation in 3D
 */

public class TestYukawaMC3D extends Simulation{
	
	public IntegratorMC integrator;
	public MCMoveAtom mcMoveAtom;
	public SpeciesSpheresMono species;
	public Phase phase;
	public P2Yukawa potential;
	public Controller controller;
	
	public TestYukawaMC3D(){
		this(500);
	}
	
	public TestYukawaMC3D(int numAtoms){
		super(Space3D.getInstance(), false, new PotentialMasterCell(Space3D.getInstance()));
		
		defaults.makeLJDefaults();
		
		
		integrator = new IntegratorMC(this);
		mcMoveAtom = new MCMoveAtom(this);
		mcMoveAtom.setAtomSource(new AtomSourceRandomLeaf());
		mcMoveAtom.setStepSize(0.2*defaults.atomSize);
		integrator.getMoveManager().addMCMove(mcMoveAtom);
		integrator.setEquilibrating(false);
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
		getController().addAction(activityIntegrate);
		species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
		phase.getAgent(species).setNMolecules(numAtoms);
		phase = new Phase(this);
		phase.setDensity(0.65);
		potential = new P2Yukawa(this);
		double truncationRadius = 3.0*potential.getKappa();
		if(truncationRadius > 0.5*phase.getBoundary().getDimensions().x(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*phase.getBoundary().getDimensions().x(0));
		}
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
		((PotentialMasterCell)potentialMaster).setCellRange(3);
		((PotentialMasterCell)potentialMaster).setRange(potentialTruncated.getRange());
		potentialMaster.addPotential(potentialTruncated, new Species[] {species, species});
			
		integrator.getMoveEventManager().addListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).makeMCMoveListener());
		
		new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(phase);
		integrator.setPhase(phase);
		
		((PotentialMasterCell)potentialMaster).getNbrCellManager(phase).assignCellAll();
		
	}
	
	public static void main(String[] args){
		
		int numAtoms = 500;
		if (args.length > 0){
			numAtoms = Integer.valueOf(args[0]).intValue();
		}
		TestYukawaMC3D sim = new TestYukawaMC3D(numAtoms);
		
		sim.getDefaults().blockSize = 10;
		MeterPressure pMeter = new MeterPressure(sim.getSpace());
		pMeter.setIntegrator(sim.integrator);
		AccumulatorAverage pAccumulator = new AccumulatorAverage(sim);
		DataPump pPump = new DataPump(pMeter,pAccumulator);
		IntervalActionAdapter iaa = new IntervalActionAdapter(pPump,sim.integrator);
		iaa.setActionInterval(2*numAtoms);
		MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
		energyMeter.setPhase(sim.phase);
		AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
		DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
		energyAccumulator.setBlockSize(50);
		new IntervalActionAdapter(energyManager, sim.integrator);

		SimulationGraphic simGraphic = new SimulationGraphic(sim);
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpeciesAgent(sim.phase.getAgent(sim.species));
        simGraphic.add(nSelector);
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
		
	}
}
