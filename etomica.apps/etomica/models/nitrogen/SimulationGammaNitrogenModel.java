package etomica.models.nitrogen;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicBcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.listener.IntegratorListenerAction;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMolecular;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Pixel;



/**
 * Simulation class for nitrogen molecules
 * gamma-N2 crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationGammaNitrogenModel extends Simulation{

	
	public SimulationGammaNitrogenModel(ISpace space, int numMolecule, double temperature, double pressure) {
		super(space);
		this.space = space;
		double a = 3.957;
		double c = 5.109;
		int nCell = (int)Math.round(Math.pow((numMolecule/2), 1.0/3.0));
	
		potentialMaster = new PotentialMaster();
				
		Basis basisBCC = new BasisCubicBcc();
		Basis basis = new BasisBigCell(space, basisBCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		//double L = Math.pow(4.0/1.0, 1.0/3.0);  // 4/1 is numAtom / density
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		int [] nCells = new int[]{1,1,1};
		
		IVector[] boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nCell*a, 0, 0});
		boxDim[1] = space.makeVector(new double[]{0, nCell*a, 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nCell*c});
		
		Boundary boundary = new BoundaryDeformablePeriodic(space,boxDim);
		primitive = new PrimitiveTetragonal(space, nCell*a, nCell*c);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsGamma();
		coordinateDef.setOrientationVectorGamma(space);
		coordinateDef.initializeCoordinates(nCells);
		
		box.setBoundary(boundary);
		double rC = 9.0;//box.getBoundary().getBoxSize().getX(0)*0.5;
		System.out.println("Truncation Radius: " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
		
		MCMoveVolume mcMoveVolume = new MCMoveVolume(this, potentialMaster, space);
		mcMoveVolume.setBox(box);
		pressure *= 1e9;
		mcMoveVolume.setPressure(Pascal.UNIT.toSim(pressure));
		
		integrator = new IntegratorMC(this, potentialMaster);
		integrator.getMoveManager().addMCMove(move);
		integrator.getMoveManager().addMCMove(rotate);
		integrator.getMoveManager().addMCMove(mcMoveVolume);
		integrator.setBox(box);
		
		integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
		
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
	}
	
	public static void main (String[] args){
		
		int numMolecule =54;
		double temperature = 35.0; // in Unit Kelvin
		double pressure = 0.32; //in Unit GPa
		long simSteps = 10000000;
		
		if(args.length > 1){
			simSteps = Long.parseLong(args[1]);
		}
		if(args.length > 2){
			temperature = Double.parseDouble(args[2]);
		}
		if(args.length > 3){
			pressure = Double.parseDouble(args[3]);
		}
		if(args.length > 4){
			numMolecule = Integer.parseInt(args[4]);
		}
		String filename = "gammaN2_nA"+numMolecule+"_T"+Math.round(temperature);
		
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running gamma-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa\n");

		
		SimulationGammaNitrogenModel sim = new SimulationGammaNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, pressure);

	    /*
	    // set up normal mode meter
//	    MeterNormalMode meterNormalMode = new MeterNormalMode();
//	    meterNormalMode.setCoordinateDefinition(sim.coordinateDef);
//	    WaveVectorFactory waveVectorFactory = new WaveVectorFactorySimple(sim.primitive, sim.space);
//	    meterNormalMode.setWaveVectorFactory(waveVectorFactory);
//	    meterNormalMode.setBox(sim.box);
//	    
//	    IntegratorListenerAction meterNormalModeListerner = new IntegratorListenerAction(meterNormalMode);
//	    meterNormalModeListerner.setInterval(numMolecule);
//	    sim.integrator.getEventManager().addListener(meterNormalModeListerner);
	    	    
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy (per molecule): "+ Kelvin.UNIT.fromSim(latticeEnergy)/numMolecule);
		//System.exit(1);
		
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyListener);
		
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();
//		meterNormalMode.reset();
		
//		WriteS sWriter = new WriteS(sim.space);
//		sWriter.setFilename(filename);
//		sWriter.setOverwrite(true);
//		sWriter.setMeter(meterNormalMode);
//		sWriter.setWaveVectorFactory(waveVectorFactory);
//		sWriter.setTemperature(temperature);
//		
//		IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
//		sWriterListener.setInterval((int)simSteps/10);
//		sim.integrator.getEventManager().addListener(sWriterListener);
		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//			PDBWriter pdbWriter = new PDBWriter(sim.box);
//			pdbWriter.setFileName(filename+".pdb");
//		    pdbWriter.actionPerformed();
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
	
//		double A = sWriter.getLastA();
//		System.out.println("A/N: " + A/numMolecule);
		System.out.println("Average energy (per molecule): "   + Kelvin.UNIT.fromSim(averageEnergy)/numMolecule  + " ;error: " + Kelvin.UNIT.fromSim(errorEnergy)/numMolecule);
	    
	    System.out.println("Box Dimension: " + sim.box.getBoundary().getBoxSize().toString());
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
		*/
		if(true){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
		    simGraphic.makeAndDisplayFrame("Alpha-Phase Nitrogen Crystal Structure");
			sim.activityIntegrate.setMaxSteps(simSteps);
			sim.getController().actionPerformed();
		}
			
	}

	protected Box box;
	protected ISpace space;
	protected PotentialMaster potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected PotentialMolecular potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	private static final long serialVersionUID = 1L;
}
