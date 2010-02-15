package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterVolume;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.listener.IntegratorListenerAction;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
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
 * alpha-N2 crystal Structure
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationNitrogenModel extends Simulation{

	
	public SimulationNitrogenModel(ISpace space, int numMolecule, double temperature, double pressure, double scale) {
		super(space);
		this.space = space;
		double unitCellLength =scale*5.661;
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
			
		potentialMaster = new PotentialMaster();
				
		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		int [] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryDeformablePeriodic(space,nCell*unitCellLength);
		primitive = new PrimitiveCubic(space, nCell*unitCellLength);
	
		System.out.println("density: " + 4/(unitCellLength*unitCellLength*unitCellLength));
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
	
		box.setBoundary(boundary);
		double rCScale = 0.45;
		double rC =box.getBoundary().getBoxSize().getX(0)*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		potential = new P2Nitrogen(space, rCScale);
		potential.setBox(box);
		System.out.println("Box Dimension(before): " + box.getBoundary().getBoxSize().toString());
		final IVector initBox = space.makeVector(new double[]{box.getBoundary().getBoxSize().getX(0),
															  box.getBoundary().getBoxSize().getX(1),
															  box.getBoundary().getBoxSize().getX(2)});
		coordinateDef.setInitVolume(initBox);
		
		
		P0LatticeEnergyCorrec p0correc = new P0LatticeEnergyCorrec(space);
		p0correc.setSpecies(species);
		p0correc.setBox(box);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		potentialMaster.addPotential(p0correc, new ISpecies[]{species, species});
		
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
		
		int numMolecule =32;
		double temperature =35.0; // in Unit Kelvin
		double pressure = 0.0; // in Unit GPa
		long simSteps = 100000;
		double scale = 1.0;
		
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
		if(args.length > 5){
			scale = Double.parseDouble(args[5]);
		}
		String filename;
		if(temperature <1.0){
			filename = "alphaN2_nA"+numMolecule+"_T0"+(int)(temperature*10);
			
		} else {
			filename = "alphaN2_nA"+numMolecule+"_T"+Math.round(temperature);
		}
		
		filename = "NPT"+numMolecule+"T"+(int)temperature;
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running alpha-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa");
		System.out.println("With volume scaling of " + scale);
		System.out.println("Output file: " + filename + "\n");

		SimulationNitrogenModel sim = new SimulationNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, pressure, scale);
	 
//		MeterPressureByVolumeChange meterPressure = new MeterPressureByVolumeChange(sim.space);
//		//meterPressure.setIncludeLrc(true);
//		meterPressure.setIntegrator(sim.integrator);
//		AccumulatorAverageFixed pressureAverage = new AccumulatorAverageFixed();
//		pressureAverage.setBlockSize(10);
//		DataPump pumpPressure = new DataPump(meterPressure, pressureAverage);
		
//		IntegratorListenerAction meterPressureListener = new IntegratorListenerAction(pumpPressure);
//		meterPressureListener.setInterval(100);
//		sim.integrator.getEventManager().addListener(meterPressureListener);
		
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergySim = Kelvin.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar());
		System.out.println("1K = "+ Kelvin.UNIT.toSim(1)+ " sim unit");
		System.out.println("Lattice Energy Sim in K (per molecule): "+ latticeEnergySim/numMolecule);
		
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyListener);
		
	    // set up normal mode meter
	    MeterNormalMode meterNormalMode = new MeterNormalMode();
	    meterNormalMode.setCoordinateDefinition(sim.coordinateDef);
	    WaveVectorFactory waveVectorFactory = new WaveVectorFactorySimple(sim.primitive, sim.space);
	    meterNormalMode.setWaveVectorFactory(waveVectorFactory);
	    meterNormalMode.setBox(sim.box);
	    
	    IntegratorListenerAction meterNormalModeListerner = new IntegratorListenerAction(meterNormalMode);
	    meterNormalModeListerner.setInterval(numMolecule);
	    sim.integrator.getEventManager().addListener(meterNormalModeListerner);
		
		//MeterNormalizedCoord meterCoord = new MeterNormalizedCoord(sim.box, sim.coordinateDef, sim.species);
		MeterNormalizedCoord meterCoord = new MeterNormalizedCoord(sim.box, sim.coordinateDef, sim.species, true);
		
		IntegratorListenerAction meterCoordListener = new IntegratorListenerAction(meterCoord);
		meterCoordListener.setInterval(1);
		sim.integrator.getEventManager().addListener(meterCoordListener);
	
		MeterVolume meterVolume = new MeterVolume();
		meterVolume.setBox(sim.box);
		AccumulatorAverageFixed volumeAverage = new AccumulatorAverageFixed();
		volumeAverage.setBlockSize(10);
		DataPump pumpVolume = new DataPump(meterVolume, volumeAverage);
		IntegratorListenerAction meterVolumeListener = new IntegratorListenerAction(pumpVolume);
		meterVolumeListener.setInterval(20);
		sim.integrator.getEventManager().addListener(meterVolumeListener);		
		
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();
		meterNormalMode.reset();
		
		WriteS sWriter = new WriteS(sim.space);
		sWriter.setFilename(filename);
		sWriter.setOverwrite(true);
		sWriter.setMeter(meterNormalMode);
		sWriter.setWaveVectorFactory(waveVectorFactory);
		sWriter.setTemperature(Kelvin.UNIT.toSim(temperature));
		
		IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
		sWriterListener.setInterval((int)simSteps/10);
		sim.integrator.getEventManager().addListener(sWriterListener);
//		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//			PDBWriter pdbWriter = new PDBWriter(sim.box);
//			pdbWriter.setFileName(filename+".pdb");
//		    pdbWriter.actionPerformed();
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
	
		double A = Kelvin.UNIT.fromSim(sWriter.getLastA());
		System.out.println("Aharm/N in K: " + A/numMolecule);
		System.out.println("Helmholtz Free energy per molecule in K: " + ((A+latticeEnergySim)/numMolecule));
		
		System.out.println("Average energy in K (per molecule): " + (Kelvin.UNIT.fromSim(averageEnergy))/numMolecule  
				+ " ;error: " + Kelvin.UNIT.fromSim(errorEnergy)/numMolecule);
				
//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
//	
//		System.out.println("Average pressure in GPa: "   + Pascal.UNIT.fromSim(averagePressure)*1e-9  
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)*1e-9);
//	 	System.out.println("density: " + numMolecule/sim.box.getBoundary().volume());
 		sim.writeUdistribution(filename, meterCoord);
	 
		double averageVolume = ((DataGroup)volumeAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorVolume = ((DataGroup)volumeAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		System.out.println("Average Volume: " + averageVolume + " ;error: " + errorVolume);
		
		double scaling = Math.pow(averageVolume, 1.0/3.0)/(sim.primitive.getSize()[0]);
		System.out.println("scaling: " + scaling);
		
		System.out.println("Box Dimension(after): " + sim.box.getBoundary().getBoxSize().toString());
		System.out.println("Density: " + numMolecule/sim.box.getBoundary().volume());
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
		
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
		    simGraphic.makeAndDisplayFrame("Alpha-Phase Nitrogen Crystal Structure");
			sim.activityIntegrate.setMaxSteps(simSteps);
			//sim.getController().actionPerformed();
		}
	}
	
	void writeUdistribution(String filename, MeterNormalizedCoord meterCoord){
		DataGroup uData = (DataGroup)meterCoord.getData();
		for (int i=0; i<uData.getNData(); i++){
			String fName = filename+"U"+i+".coord";
			try {
				FileWriter fileWriter = new FileWriter(fName,false);
				
				DataDoubleArray uDistribution = (DataDoubleArray)uData.getData(i);
				
				for (int j=0; j<uDistribution.getLength()/uDistribution.getArrayDimension(); j++){
				
					fileWriter.write(uDistribution.getValue(new int[]{0,j})+" "+ uDistribution.getValue(new int[]{1,j}) + "\n");
				}
			
				fileWriter.close();
				
			} catch(IOException e){
				throw new RuntimeException("Failed to write coord data normalize coord U" + e);
			
			}
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
	protected SpeciesN2 species;
	private static final long serialVersionUID = 1L;
}
