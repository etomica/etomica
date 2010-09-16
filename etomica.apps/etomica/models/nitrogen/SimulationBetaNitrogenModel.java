package etomica.models.nitrogen;

import etomica.action.BoxInflate;
import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.DiameterHashByType;
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

	
	public SimulationBetaNitrogenModel(ISpace space, int[] nC, double temperature, double pressure, double newScale, double density) {
		super(space);
		this.space = space;
		
    	double ratio = 1.631;
		double a = Math.pow(4.0/(Math.sqrt(3.0)*ratio*density), 1.0/3.0);
		double c = a*ratio;
		System.out.println("\n\na: " + a + " ;cDim: " + c);
		int numMolecule = nC[0]*nC[1]*nC[2]*2;
	
		potentialMaster = new PotentialMaster();
		Basis basisHCP = new BasisHcp();
		Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC[0], nC[1], nC[2]});
		
		species = new SpeciesN2(space);
		addSpecies(species);
		
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
		coordinateDef.setIsBeta();
		coordinateDef.setOrientationVectorBeta(space);
		coordinateDef.initializeCoordinates(nCells);

		box.setBoundary(boundary);
		double rC = a*nC[0]*0.475;
		System.out.println("Truncation Radius: " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);
		
		PRotConstraint pRot= new PRotConstraint(space, coordinateDef, box);
		pRot.setBox(box);
		pRot.setConstraintAngle(70);
		
		BoxInflate boxInflate = new BoxInflate(box, space);
		boxInflate.setScale(newScale);
		boxInflate.actionPerformed();

//		ConfigurationFile configFile = new ConfigurationFile("configFile"+numMolecule);
//		configFile.initializeCoordinates(box);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		potentialMaster.addPotential(pRot, new ISpecies[]{species});
		
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
		
		int nC0 =6; 
		int nC1 =6; 
		int nC2 =6;
		double temperature =40; // in Unit Kelvin
		double pressure = 0.0; //in Unit GPa
		long simSteps = 100000;
		double newScale = 1.0;
		double density = 0.025;
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
		if(args.length > 7){
			newScale = Double.parseDouble(args[7]);
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
		System.out.println("With volume scaling of " + newScale);
		System.out.println("Output file: " + filename + "\n");

		
		SimulationBetaNitrogenModel sim = new SimulationBetaNitrogenModel(Space3D.getInstance(3), nC, temperature, pressure, newScale, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		System.out.println("Lattice Energy (per molecule): "+ meterPotentialEnergy.getDataAsScalar());
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
		
		
		if(true){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    simGraphic.makeAndDisplayFrame("Beta-Phase Nitrogen Crystal Structure");
		    
		    DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
			diameter.setDiameter(sim.species.getPType(), 0.0);
			
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
			IAction output = new IAction(){

				public void actionPerformed() {
					System.out.println("energy: " + (meterPotentialEnergy.getDataAsScalar()+329916.40502540825));
					
				}
				
			};
			
//			IntegratorListenerAction outListener = new IntegratorListenerAction(output);
//			outListener.setInterval(numMolecule);
//			sim.integrator.getEventManager().addListener(outListener);
			return;
		}
		
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();

		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
//		double averagePressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//		double errorPressure = ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
		System.out.println("Average energy (per molecule): "   + averageEnergy/numMolecule  
				+ " ;error: " + errorEnergy/numMolecule);
//		System.out.println("Average pressure (GPa): " + Pascal.UNIT.fromSim(averagePressure)/1e9 
//				+ " ;error: " + Pascal.UNIT.fromSim(errorPressure)/1e9);

		double  a = sim.box.getBoundary().getEdgeVector(0).getX(0)/nC0;
		double  c = sim.box.getBoundary().getEdgeVector(2).getX(2)/nC2;
		System.out.println("\na: " +a + " ;c: "+c +" ;c/a: " + (c/a));
		double scaling = sim.box.getBoundary().getEdgeVector(0).getX(0)/(nC0*3.854);
	    System.out.println("scaling: " + scaling);
	    long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
			

			
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
