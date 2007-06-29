package etomica.yukawa;

import etomica.action.Action;
import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.box.Box;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;

/**
 * Hard-core plus two Yukawa fluid Monte-Carlo simulation in 3D
 */

public class TestHC2YukawaMC3D extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "Test HC2 Yukawa MC3D";
    public IntegratorMC integrator;
	public MCMoveAtom mcMoveAtom;
	public SpeciesSpheresMono species;
	public Box box;
	public P2HC2Yukawa potential;
	public Controller controller;
	
	public TestHC2YukawaMC3D(){
		this(500);
	}
	
	public TestHC2YukawaMC3D(int numAtoms){
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
        box.getAgent(species).setNMolecules(numAtoms);
		box.setDensity(0.65);
		potential = new P2HC2Yukawa(this);
		double truncationRadius = 3.0*potential.getSigma();
		if(truncationRadius > 0.5*box.getBoundary().getDimensions().x(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*box.getBoundary().getDimensions().x(0));
		}
		
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
		potentialMaster.setCellRange(3);
		potentialMaster.setRange(potentialTruncated.getRange());
		potentialMaster.addPotential(potentialTruncated, new Species[] {species, species});
			
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
		TestHC2YukawaMC3D sim = new TestHC2YukawaMC3D(numAtoms);
		
		MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
		AccumulatorAverage energyAccumulator = new AccumulatorAverage(10);
		DataPump energyManager = new DataPump(energyMeter, energyAccumulator);
		energyAccumulator.setBlockSize(50);
        sim.integrator.addIntervalAction(energyManager);
		
		final SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
		Action repaintAction = simGraphic.getDisplayBoxPaintAction(sim.box);

        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setPostAction(repaintAction);
        simGraphic.getController().getReinitButton().setPostAction(repaintAction);

        nSelector.setSpeciesAgent(sim.box.getAgent(sim.species));
        simGraphic.add(nSelector);
        simGraphic.makeAndDisplayFrame(APP_NAME);
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);

/*
		double Z = ((DataDouble)((DataGroup)pAccumulator.getData()).getData(StatType.AVERAGE.index)).x*sim.box.volume()/(sim.box.moleculeCount().sim.integrator.getTemperature());
		double avgPE = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.AVERAGE.index)).x;
		avgPE /= numAtoms;
		System.out.println("Z="+Z);
		System.out.println("PE"+avgPE);
		double temp = sim.integrator.getTemperature();
		double Cv = ((DataDouble)((DataGroup)energyAccumulator.getData()).getData(StatType.STANDARD_DEVIATION.index)).x;
		Cv /= temp;
		Cv *= Cv/numAtoms;
		System.out.println("Cv/k="+Cv);
		
		if (Double.isNan(Z) || Math.abs(Z-0.15) > 0.15) {
			System.exit(1);
		}
		
		if (Double.isNaN(avgPE) || Math.abs(avgPE+4.56) > 0.03){
			System.exit(1);	
		}
		
		if (Double.isNaN(Cv) || Math.abs(Cv-0.61) > 0.45) {
			System.exit(1)
		}
*/
	}
}
