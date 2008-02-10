package etomica.yukawa;
import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

/**
 * A Yukawa Monte-Carlo simulation in 3D
 */

public class TestYukawaMC3D extends Simulation{
	
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "Test Yukawa MC3D";
    public IntegratorMC integrator;
	public MCMoveAtom mcMoveAtom;
	public SpeciesSpheresMono species;
	public Box box;
	public P2Yukawa potential;
	public Controller controller;
	
	public TestYukawaMC3D(){
		this(500);
	}
	
	public TestYukawaMC3D(int numAtoms){
		super(Space3D.getInstance(), false);
        PotentialMasterCell potentialMaster = new PotentialMasterCell(this);
		
		integrator = new IntegratorMC(this, potentialMaster);
		mcMoveAtom = new MCMoveAtom(this, potentialMaster);
		mcMoveAtom.setStepSize(0.2);
		integrator.getMoveManager().addMCMove(mcMoveAtom);
		integrator.getMoveManager().setEquilibrating(false);
		ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
		species = new SpeciesSpheresMono(this);
        getSpeciesManager().addSpecies(species);
		box = new Box(this);
        addBox(box);
		box.setNMolecules(species, numAtoms);
		box.setDensity(0.65);
		potential = new P2Yukawa(this);
		double truncationRadius = 3.0*potential.getKappa();
		if(truncationRadius > 0.5*box.getBoundary().getDimensions().x(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*box.getBoundary().getDimensions().x(0));
		}
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
		potentialMaster.setCellRange(3);
		potentialMaster.setRange(potentialTruncated.getRange());
		potentialMaster.addPotential(potentialTruncated, new ISpecies[] {species, species});
			
		integrator.getMoveEventManager().addListener(potentialMaster.getNbrCellManager(box).makeMCMoveListener());
		
		new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(box);
		integrator.setBox(box);
		
		potentialMaster.getNbrCellManager(box).assignCellAll();
		
	}
	
	public static void main(String[] args){
		
		int numAtoms = 500;
		if (args.length > 0){
			numAtoms = Integer.valueOf(args[0]).intValue();
		}
		TestYukawaMC3D sim = new TestYukawaMC3D(numAtoms);
		
		MeterPressure pMeter = new MeterPressure(sim.getSpace());
		pMeter.setIntegrator(sim.integrator);
		AccumulatorAverageCollapsing pAccumulator = new AccumulatorAverageCollapsing();
		DataPump pPump = new DataPump(pMeter,pAccumulator);
        sim.integrator.addIntervalAction(pPump);
        sim.integrator.setActionInterval(pPump, 2*numAtoms);
		MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
		AccumulatorAverageCollapsing energyAccumulator = new AccumulatorAverageCollapsing();
		DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
		energyAccumulator.setBlockSize(50);
        sim.integrator.addIntervalAction(energyManager);

		final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
		Action repaintAction = simGraphic.getPaintAction(sim.box);

        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setPostAction(repaintAction);
        nSelector.setSpecies(sim.species);
        nSelector.setBox(sim.box);
        simGraphic.add(nSelector);
        simGraphic.getController().getReinitButton().setPostAction(repaintAction);
        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
		
	}
}
