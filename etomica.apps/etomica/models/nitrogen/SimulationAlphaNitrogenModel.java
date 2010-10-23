package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.atom.DiameterHashByType;
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
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMolecular;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.units.Pixel;



/**
 * Simulation class for nitrogen molecules
 * alpha-N2 crystal Structure
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationAlphaNitrogenModel extends Simulation{

	
	public SimulationAlphaNitrogenModel(ISpace space, int numMolecule, double temperature, double density) {
		super(space);
		this.space = space;
		int nCell = (int) Math.round(Math.pow((numMolecule/4), 1.0/3.0));
			
		double a = Math.pow(numMolecule/density, 1.0/3.0)/nCell;
		System.out.println("Unit Cell Length, a: " + a);
		
		potentialMaster = new PotentialMaster();
	//	potentialMaster = new PotentialMasterListMolecular(this, space);
				
		
		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nCell, nCell, nCell});
		
		species = new SpeciesN2(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		
		int [] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryRectangularPeriodic(space,nCell*a);
		primitive = new PrimitiveCubic(space, nCell*a);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
	
		box.setBoundary(boundary);
		double rCScale = 0.475;
		double rC =box.getBoundary().getBoxSize().getX(0)*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		PRotConstraint pRotConstraint = new PRotConstraint(space,coordinateDef,box);
		pRotConstraint.setConstraintAngle(70);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		//potentialMaster.addPotential(pRotConstraint,new ISpecies[]{species} );
		//potentialMaster.lrcMaster().isEnabled();
		
//	    int cellRange = 6;
//      potentialMaster.setRange(rC);
//      potentialMaster.setCellRange(cellRange); 
//      potentialMaster.getNeighborManager(box).reset();
//      
//      int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
//      if (potentialCells < cellRange*2+1) {
//          throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
//      }
//	
//      int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
//      System.out.println("numNeigh: " + numNeigh);
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
			
		integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature));
		integrator.getMoveManager().addMCMove(move);
		integrator.getMoveManager().addMCMove(rotate);
		integrator.setBox(box);
		
		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
		meterPE.setBox(box);
		System.out.println("lattice: " + meterPE.getDataAsScalar()/numMolecule);
		//System.exit(1);
		
//		double[] u = new double[coordinateDef.getCoordinateDim()];
//		
//		for (double i=-0.5; i<=0.51; i+=0.01){
//			u[0] = i;
//			coordinateDef.setToU(box.getMoleculeList(), u);
//			double energy = meterPE.getDataAsScalar();
//			
//			System.out.println(i+" "+ energy/numMolecule);
//		}
//		
//		System.exit(1);
		
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
	}
	
	public static void main (String[] args){
		
		int numMolecule = 4*4*4*4;
		double temperature = 1; // in Unit Kelvin
		long simSteps = 100000;
		double density = 0.025; //0.02204857502170207 (intial from literature with a = 5.661)
		
		if(args.length > 0){
			simSteps = Long.parseLong(args[0]);
		}
		if(args.length > 1){
			temperature = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			numMolecule = Integer.parseInt(args[2]);
		}
		String filename = "alphaN2_nA"+numMolecule+"_T"+temperature;
		System.out.println("Running alpha-N2 NVT simulation with " + simSteps + " steps" );
		System.out.println("number Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; density: "+ density +"\n");

		SimulationAlphaNitrogenModel sim = new SimulationAlphaNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, density);
	 
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergySim = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy Sim (per molecule): "+ latticeEnergySim);
					
//		MeterNormalizedCoord meterCoord = new MeterNormalizedCoord(sim.box, sim.coordinateDef, sim.species);
//        IntegratorListenerAction meterCoordListener = new IntegratorListenerAction(meterCoord);
//        meterCoordListener.setInterval(numMolecule);                                      
//        sim.integrator.getEventManager().addListener(meterCoordListener);       
//		
//		MeterOrientationDistribution meterOrient = new MeterOrientationDistribution(sim.box, sim.coordinateDef, sim.species);
//        IntegratorListenerAction meterOrientListener = new IntegratorListenerAction(meterOrient);
//        meterOrientListener.setInterval(numMolecule);                                      
//        sim.integrator.getEventManager().addListener(meterOrientListener);
		
//		MeterPressureMolecular meterPressure = new MeterPressureMolecular(sim.space);
//		meterPressure.setIntegrator(sim.integrator);
//						
//		double staticPressure = meterPressure.getDataAsScalar();
//		System.out.println("Static Pressure (GPa): " + Pascal.UNIT.fromSim(staticPressure)/1e9);

		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    	    
		    
			DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
			diameter.setDiameter(sim.species.getPType(), 0.0);
			
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
		    simGraphic.makeAndDisplayFrame("Alpha-Phase Nitrogen Crystal Structure");
		    
			
			IAction output = new IAction(){

				public void actionPerformed() {
					System.out.println("energy: " + (meterPotentialEnergy.getDataAsScalar()+82512.1706428608));
					
				}
				
			};
			
			IntegratorListenerAction outListener = new IntegratorListenerAction(output);
			outListener.setInterval(1);
			sim.integrator.getEventManager().addListener(outListener);
			
//			double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//			double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
			
			//sim.activityIntegrate.setMaxSteps(1000000);
			//sim.getController().actionPerformed();
		    return;
		    
		}
		//System.exit(1);
			
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval((int)simSteps/100);
		sim.integrator.getEventManager().addListener(energyListener);
			
//		AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
//		DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//		IntegratorListenerAction pressureListener = new IntegratorListenerAction(pressurePump);
//		pressureListener.setInterval((int)simSteps/100);
//		sim.integrator.getEventManager().addListener(pressureListener);
			
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
			
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();
	
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();

//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
		System.out.println("Average energy (per molecule): " + (averageEnergy)/numMolecule  
				+ " ;error: " + errorEnergy/numMolecule);
//		System.out.println("Average pressure (GPa): " + Pascal.UNIT.fromSim(averagePressure)/1e9 
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)/1e9);
//		sim.writeUdistribution(filename, meterCoord);
//		meterOrient.writeUdistribution(filename+"Ort");
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken(s): " + (endTime - startTime)/1000);
		
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
