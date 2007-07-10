package etomica.densityofstate;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.box.Box;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistogramExpanding;
import etomica.yukawa.P2Yukawa;

/**
 * A Yukawa Monte-Carlo simulation in 3D
 */

public class DensityOfState extends Simulation{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public IntegratorMC integrator;
	public MCMoveAtom mcMoveAtom;
	public SpeciesSpheresMono species;
	public Box box;
	public P2Yukawa potential;
	public Controller controller;
	public ActivityIntegrate activityIntegrate; 
	
	public DensityOfState(){
		this(500);
	}
	
	
	
	public DensityOfState(int numAtoms){
		super(Space3D.getInstance(), false);
		
		
		potentialMaster = new PotentialMaster(space);
		integrator = new IntegratorMC(this, potentialMaster);
		mcMoveAtom = new MCMoveAtom(this, potentialMaster);
		mcMoveAtom.setAtomSource(new AtomSourceRandomLeaf());
		mcMoveAtom.setStepSize(0.2);
		integrator.getMoveManager().addMCMove(mcMoveAtom);
		integrator.getMoveManager().setEquilibrating(false);
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
		species = new SpeciesSpheresMono(this);
		box.setNMolecules(species, numAtoms);
		box = new Box(this);
		box.setDensity(0.65);
		potential = new P2Yukawa(this);
		double truncationRadius = 3.0*potential.getKappa();
		if(truncationRadius > 0.5*box.getBoundary().getDimensions().x(0)){
			throw new RuntimeException("Truncaiton radius too large.  Max allowed is "+0.5*box.getBoundary().getDimensions().x(0));
		}
		P2SoftSphericalTruncated potentialTruncated = new P2SoftSphericalTruncated(potential, truncationRadius);
		((PotentialMasterCell)potentialMaster).setCellRange(3);
		((PotentialMasterCell)potentialMaster).setRange(potentialTruncated.getRange());
		potentialMaster.addPotential(potentialTruncated, new Species[] {species, species});
			
		integrator.getMoveEventManager().addListener(((PotentialMasterCell)potentialMaster).getNbrCellManager(box).makeMCMoveListener());
		
		new ConfigurationLattice(new LatticeCubicFcc()).initializeCoordinates(box);
		integrator.setBox(box);
		
		((PotentialMasterCell)potentialMaster).getNbrCellManager(box).assignCellAll();
		
	}
	public PotentialMaster potentialMaster;
	
	public static void main(String[] args){
		
		int numAtoms = 500;
		long maxSteps = 10000000;
		double temperature = 6.0;
		
		if (args.length > 0){
			numAtoms = Integer.valueOf(args[0]).intValue();
		}
		if (args.length > 1){
			temperature = Double.parseDouble(args[1]);
		}
		if (args.length > 2) {
			maxSteps = Long.parseLong(args[2]);
		}
		 
		DensityOfState sim = new DensityOfState(numAtoms);
		sim.activityIntegrate.setMaxSteps(maxSteps);
		sim.integrator.setTemperature(temperature);
		sim.getController().actionPerformed();
				
		MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
		AccumulatorHistogram histogram = new AccumulatorHistogram(new HistogramExpanding(1), 100); //bin size =1
		
		histogram.setPushInterval(100);
		DataPump energyManager = new DataPump(energyMeter, histogram);
		

		sim.activityIntegrate.setMaxSteps(maxSteps);
		sim.integrator.addIntervalAction(energyManager);
		sim.getController().reset();
		sim.getController().actionPerformed();
		
		DataLogger dataLogger = new DataLogger();
		DataTableWriter dataTableWriter = new DataTableWriter();
		dataLogger.setFileName("histogram" + "@" + (int)(temperature*10) + ".dat");
		dataLogger.setDataSink(dataTableWriter);
		dataLogger.putDataInfo(histogram.getDataInfo());
		
		dataLogger.setWriteInterval(1);
		dataLogger.putData(histogram.getData());
		dataLogger.closeFile();

		
		
		
		
		//DisplayPlot HistogramDisplay = new DisplayPlot();
		//histogram.addDataSink(HistogramDisplay.getDataSet().makeDataSink());
		
		//SimulationGraphic simGraphic = new SimulationGraphic(sim);
        //simGraphic.makeAndDisplayFrame();
        //simGraphic.panel().add(HistogramDisplay.graphic());
        //ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        //colorScheme.setColor(sim.species.getMoleculeType(), java.awt.Color.red);
		
	}
}
