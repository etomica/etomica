package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
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
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
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
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.units.Pascal;
import etomica.units.Pixel;



/**
 * Simulation class for nitrogen molecules
 * beta-N2 crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationBetaNitrogenModel extends Simulation{

	
	public SimulationBetaNitrogenModel(ISpace space, int[] nC, double temperature, double pressure) {
		super(space);
		this.space = space;
		double a = 3.854;
		double c = 6.284; 
		int numMolecule = nC[0]*nC[1]*nC[2]*2;
	
		potentialMaster = new PotentialMaster();
		Basis basisHCP = new BasisHcp();
		Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC[0], nC[1], nC[2]});
		
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
		boxDim[0] = space.makeVector(new double[]{nC[0]*a, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC[1]*a*Math.cos(Degree.UNIT.toSim(60)), nC[1]*a*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC[2]*c});
		
		Boundary boundary = new BoundaryDeformablePeriodic(space,boxDim);
		primitive = new PrimitiveHexagonal(space, (nC[0])*a, nC[2]*c);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsGamma();
		coordinateDef.setOrientationVectorGamma(space);
		coordinateDef.initializeCoordinates(nCells);
		
		ConfigurationFile configFile = new ConfigurationFile("configFile+numMolecule");
		configFile.initializeCoordinates(box);
		
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
		//integrator.getMoveManager().addMCMove(mcMoveVolume);
		integrator.setBox(box);
		
		integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
		
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
	}
	
	public static void main (String[] args){
		
		int nC0 =7; 
		int nC1 =6; 
		int nC2 =6;
		double temperature = 50.0; // in Unit Kelvin
		double pressure = 0.0; //in Unit GPa
		long simSteps = 10000;
		
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
			nC0 = Integer.parseInt(args[4]);
		}
		if(args.length > 5){
			nC1 = Integer.parseInt(args[5]);
		}
		if(args.length > 6){
			nC2 = Integer.parseInt(args[6]);
		}
		
		int[] nC = new int []{nC0,nC1,nC2};
		int numMolecule =nC[0]*nC[1]*nC[2]*2;
		
		String filename = "betaN2_nA"+numMolecule+"_T"+Math.round(temperature);
		
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running beta-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa");
		System.out.println("Output file: " + filename + "\n");

		
		SimulationBetaNitrogenModel sim = new SimulationBetaNitrogenModel(Space3D.getInstance(3), nC, temperature, pressure);

	    
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
		
		WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
		writeConfig.setBox(sim.box);
		writeConfig.setDoApplyPBC(false);
		writeConfig.setFileName("configFile"+numMolecule+".pos");
		writeConfig.actionPerformed();
		
		writeConfig.setFileName("configBeta"+numMolecule+"T"+Math.round(temperature)+".pos");
		writeConfig.actionPerformed();
		
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(20));
		    simGraphic.makeAndDisplayFrame("Beta-Phase Nitrogen Crystal Structure");
			//sim.activityIntegrate.setMaxSteps(simSteps);
			//sim.getController().actionPerformed();
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
