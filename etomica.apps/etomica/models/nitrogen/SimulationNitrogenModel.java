package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
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

	
	public SimulationNitrogenModel(ISpace space, int numMolecule, double temperature, double pressure) {
		super(space);
		this.space = space;
		double unitCellLength =1.0*5.661;
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
			
		potentialMaster = new PotentialMaster();
				
		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		species = new SpeciesN2(space);
		species.setConformation(conformation);
		addSpecies(species);
		
		//double L = Math.pow(4.0/1.0, 1.0/3.0);  // 4/1 is numAtom / density
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
		double rC =9.0;//box.getBoundary().getBoxSize().getX(0)*0.5;
		System.out.println("Truncation Radius: " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		System.out.println("Box Dimension(before): " + box.getBoundary().getBoxSize().toString());
		final IVector initBox = space.makeVector(new double[]{box.getBoundary().getBoxSize().getX(0),
															  box.getBoundary().getBoxSize().getX(1),
															  box.getBoundary().getBoxSize().getX(2)});
		coordinateDef.setInitVolume(initBox);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
		
		MCMoveVolumeN2 mcMoveVolume = new MCMoveVolumeN2(this, potentialMaster, space);
		mcMoveVolume.setSpecies(species);
		mcMoveVolume.setBox(box);
		mcMoveVolume.setXYZChange();
		uLatticeCorrec = mcMoveVolume.getLatticeCorrec();
		
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
		double temperature = 25; // in Unit Kelvin
		double pressure = 0.0; // in Unit GPa
		long simSteps = 1000000;
		
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
		String filename;
		if(temperature <1.0){
			filename = "alphaN2_nA"+numMolecule+"_T0"+(int)(temperature*10);
			
		} else {
			filename = "alphaN2_nA"+numMolecule+"_T"+Math.round(temperature);
		}
		
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running alpha-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa\n");

		SimulationNitrogenModel sim = new SimulationNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, pressure);
	 
//		MeterPressure meterPressure = new MeterPressure(sim.space);
//		meterPressure.setIncludeLrc(true);
//		meterPressure.setIntegrator(sim.integrator);
//		AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
//		DataPump pumpPressure = new DataPump(meterPressure, pressureAverage);
//		
//		IntegratorListenerAction meterPressureListener = new IntegratorListenerAction(pumpPressure);
//		meterPressureListener.setInterval(100);
//		sim.integrator.getEventManager().addListener(meterPressureListener);
//		
//		System.out.println("Static Pressure: " + Pascal.UNIT.fromSim(meterPressure.getDataAsScalar())*1e-9);
//		
		
   
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy correction in K (per molecule): "+ sim.uLatticeCorrec/numMolecule);
		System.out.println("Lattice Energy in K (per molecule): "+ (Kelvin.UNIT.fromSim(latticeEnergy)+sim.uLatticeCorrec)/numMolecule);
		
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
		
		MeterNormalizedCoord meterCoord = new MeterNormalizedCoord(sim.box, sim.coordinateDef, sim.species);
		
		IntegratorListenerAction meterCoordListener = new IntegratorListenerAction(meterCoord);
		meterCoordListener.setInterval(1);
		sim.integrator.getEventManager().addListener(meterCoordListener);
		
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
		sWriter.setTemperature(temperature);
		
		IntegratorListenerAction sWriterListener = new IntegratorListenerAction(sWriter);
		sWriterListener.setInterval((int)simSteps/10);
		sim.integrator.getEventManager().addListener(sWriterListener);
		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//			PDBWriter pdbWriter = new PDBWriter(sim.box);
//			pdbWriter.setFileName(filename+".pdb");
//		    pdbWriter.actionPerformed();
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
	
		double A = sWriter.getLastA();
		System.out.println("Aharm/N in K: " + Kelvin.UNIT.fromSim(A)/numMolecule);
		System.out.println("Helmholtz Free energy per molecule in K: " + (Kelvin.UNIT.fromSim(A)+(Kelvin.UNIT.fromSim(latticeEnergy)+sim.uLatticeCorrec))/numMolecule);
		System.out.println("Average energy in K (per molecule): "   + (Kelvin.UNIT.fromSim(averageEnergy)+sim.uLatticeCorrec)/numMolecule  
				+ " ;error: " + Kelvin.UNIT.fromSim(errorEnergy)/numMolecule);
		
//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
//	
//		System.out.println("Average pressure: "   + Pascal.UNIT.fromSim(averagePressure)*1e-9  
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)*1e-9);
	 		
		//sim.writeUdistribution(filename, meterCoord);
	 
		System.out.println("Box Dimension(after): " + sim.box.getBoundary().getBoxSize().toString());
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
		
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
		    simGraphic.makeAndDisplayFrame("Alpha-Phase Nitrogen Crystal Structure");
			sim.activityIntegrate.setMaxSteps(simSteps);
			sim.getController().actionPerformed();
		}
		System.out.println("Box Dimension(after): " + sim.box.getBoundary().getBoxSize().toString());
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
	
	protected double uLatticeCorrec;
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
