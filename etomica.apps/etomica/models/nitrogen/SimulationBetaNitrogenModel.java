package etomica.models.nitrogen;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
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
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.nbr.list.molecule.NeighborListManagerSlantyMolecular;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
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

	
	public SimulationBetaNitrogenModel(ISpace space, int numMolecule, double temperature, double pressure, double newScale, double density) {
		super(space);
		this.space = space;
		
		BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
	    BoxAgentManager boxAgentManager = new BoxAgentManager(boxAgentSource);
	     
		
    	double ratio = 1.631;
		double a = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double c = a*ratio;
		System.out.println("\na: " + a + " ;cDim: " + c);
	
		int nC = (int)Math.pow(numMolecule/1.999999999, 1.0/3.0);
	
		Basis basisHCP = new BasisHcp();
		Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});
		
		species = new SpeciesN2(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		int [] nCells = new int[]{1,1,1};
		
		IVector[] boxDim = new IVector[3];
		boxDim[0] = space.makeVector(new double[]{nC*a, 0, 0});
		boxDim[1] = space.makeVector(new double[]{-nC*a*Math.cos(Degree.UNIT.toSim(60)), nC*a*Math.sin(Degree.UNIT.toSim(60)), 0});
		boxDim[2] = space.makeVector(new double[]{0, 0, nC*c});
		
		Boundary boundary = new BoundaryDeformablePeriodic(space,boxDim);
		primitive = new PrimitiveHexagonal(space, (nC)*a, nC*c);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);

		box.setBoundary(boundary);
		double rC = a*nC*0.475;
		System.out.println("Truncation Radius: " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		PRotConstraint pRot= new PRotConstraint(space, coordinateDef, box);
		pRot.setBox(box);
		pRot.setConstraintAngle(0.1);
		
		BoxInflate boxInflate = new BoxInflate(box, space);
		boxInflate.setScale(newScale);
		boxInflate.actionPerformed();

//		ConfigurationFile configFile = new ConfigurationFile("configFile"+numMolecule);
//		configFile.initializeCoordinates(box);
		
		//potentialMaster = new PotentialMaster();
		potentialMaster = new PotentialMasterListMolecular(this, rC, boxAgentSource, boxAgentManager, new NeighborListManagerSlantyMolecular.NeighborListSlantyAgentSourceMolecular(rC, space), space);
	    
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		//potentialMaster.addPotential(pRot, new ISpecies[]{species});
		
	    int cellRange = 6;
      potentialMaster.setRange(rC);
      potentialMaster.setCellRange(cellRange); 
      potentialMaster.getNeighborManager(box).reset();
      
      potential.setRange(Double.POSITIVE_INFINITY);
      int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
      if (potentialCells < cellRange*2+1) {
          throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
      }
	
      int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
      System.out.println("numNeigh: " + numNeigh);
		
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
		//integrator.getMoveManager().addMCMove(rotate);
		//integrator.getMoveManager().addMCMove(mcMoveVolume);
		integrator.setBox(box);
		
//		NormalModesFromFile nm = new NormalModesFromFile("beta"+numMolecule+"_2ndDer_d"+density, 3);
//		MeterHarmonicEnergy meterHarm = new MeterHarmonicEnergy(coordinateDef, nm);
//		
//		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
//		meterPE.setBox(box);
//		double latticeEnergy = meterPE.getDataAsScalar();
//		double[] u = new double[coordinateDef.getCoordinateDim()];
//		
//		for (double i=-0.5; i<=0.51; i+=0.01){
//			u[0] = i;
//			coordinateDef.setToU(box.getMoleculeList(), u);
//			double pe = meterPE.getDataAsScalar();
//			double he = meterHarm.getDataAsScalar();
//			System.out.println(i+" "+ (pe-latticeEnergy)/numMolecule + " " + he/numMolecule);
//		}
//		
//		System.exit(1);
//		
		
		
		integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
		
		activityIntegrate = new ActivityIntegrate(integrator);
		getController().addAction(activityIntegrate);
	}
	
	public static void main (String[] args){
		
		double temperature =1; // in Unit Kelvin
		double pressure = 0.0; //in Unit GPa
		long simSteps = 100000;
		double newScale = 1.0;
		double density = 0.025;
		int numMolecule = 6*6*6*2;
		if(args.length > 1){
			simSteps = Long.parseLong(args[1]);
		}
		if(args.length > 2){
			temperature = Double.parseDouble(args[2]);
		}
		if(args.length > 3){
			numMolecule = Integer.parseInt(args[3]);
		}
		String filename = "betaN2_nA"+numMolecule+"_T"+temperature;
		
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running beta-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa");
		System.out.println("With volume scaling of " + newScale);
		System.out.println("Output file: " + filename + "\n");

		
		SimulationBetaNitrogenModel sim = new SimulationBetaNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, pressure, newScale, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		final double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy (per molecule): "+ latticeEnergy/numMolecule);
		//System.exit(1);
		System.out.println("Lattice Energy: "+ latticeEnergy);
		//System.exit(1);
		
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyListener);
		
//		MeterPressureMolecular meterPressure = new MeterPressureMolecular(sim.space);
//		meterPressure.setIntegrator(sim.integrator);
//						
//		AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
//		DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
//		IntegratorListenerAction pressureListener = new IntegratorListenerAction(pressurePump);
//		pressureListener.setInterval((int)simSteps/100);
//		sim.integrator.getEventManager().addListener(pressureListener);
			
//		double staticPressure = meterPressure.getDataAsScalar();
//		System.out.println("Static Pressure (GPa): " + Pascal.UNIT.fromSim(staticPressure)/1e9);
		
		
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    simGraphic.makeAndDisplayFrame("Beta-Phase Nitrogen Crystal Structure");
		    
		    DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
			diameter.setDiameter(sim.species.getPType(), 0.0);
			
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
			IAction output = new IAction(){

				public void actionPerformed() {
					System.out.println("energy: " + (meterPotentialEnergy.getDataAsScalar()-latticeEnergy));
					
				}
				
			};
			
			IntegratorListenerAction outListener = new IntegratorListenerAction(output);
			outListener.setInterval(numMolecule);
			sim.integrator.getEventManager().addListener(outListener);
			return;
		}
		
//		MeterOrientationDistribution meterOrient = new MeterOrientationDistribution(sim.box, sim.coordinateDef, sim.species);
//        IntegratorListenerAction meterOrientListener = new IntegratorListenerAction(meterOrient);
//        meterOrientListener.setInterval(numMolecule);                                      
//        sim.integrator.getEventManager().addListener(meterOrientListener);       
		
        sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();

		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//		sim.writeUdistribution(filename, meterOrient);
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
		System.out.println("Average energy (per molecule): "   + averageEnergy/numMolecule  
				+ " ;error: " + errorEnergy/numMolecule);
//		System.out.println("Average pressure (GPa): " + Pascal.UNIT.fromSim(averagePressure)/1e9 
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)/1e9);

	    long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
			

			
	}
	
	void writeUdistribution(String filename, MeterOrientationDistribution meterOrient){
		DataGroup uData = (DataGroup)meterOrient.getData();
		
		for (int i=0; i<uData.getNData(); i++){
			String fName = filename+"U"+i+".orient";
			try {
				FileWriter fileWriter = new FileWriter(fName,false);
				
				DataDoubleArray uDistribution = (DataDoubleArray)uData.getData(i);
				
				for (int j=0; j<uDistribution.getLength()/uDistribution.getArrayDimension(); j++){
					fileWriter.write(uDistribution.getValue(new int[]{0,j})+" "+ uDistribution.getValue(new int[]{1,j}) + "\n");
				}
			
				fileWriter.close();
				
			} catch(IOException e){
				throw new RuntimeException("Failed to write coord data orientation U" + e);
			
			}
		}
		
	}
	protected Box box;
	protected ISpace space;
	protected PotentialMasterListMolecular potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	private static final long serialVersionUID = 1L;
}
