package etomica.models.nitrogen;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.ISpecies;
import etomica.box.Box;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
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
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationNitrogenModel extends Simulation{

	
	public SimulationNitrogenModel(ISpace space, int nCell, double temperature) {
		super(space);
		this.space = space;
		double unitCellLength = 5.661;
		int numAtoms = 4*nCell*nCell*nCell;
		
		potentialMaster = new PotentialMaster();
				
		
		Basis basisFCC = new BasisCubicFcc();
		Primitive primitiveImag = new PrimitiveCubic(space, unitCellLength);
		Basis basis = new BasisBigCell(space, primitiveImag, basisFCC, new int[]{nCell, nCell, nCell});
		
		ConformationNitrogen conformation = new ConformationNitrogen(space);
		SpeciesN2 species = new SpeciesN2(space);
		species.setConformation(conformation);
		getSpeciesManager().addSpecies(species);
		
		//double L = Math.pow(4.0/1.0, 1.0/3.0);  // 4/1 is numAtom / density
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numAtoms);		
		
		int [] nCells = new int[]{1,1,1};
		Boundary boundary = new BoundaryRectangularPeriodic(space,nCell*unitCellLength);
		Primitive primitive = new PrimitiveCubic(space, nCell*unitCellLength);
		
		CoordinateDefinitionNitrogen coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setOrientationVector(space);
		coordinateDef.initializeCoordinates(nCells);
		
		box.setBoundary(boundary);
		potential = new P2Nitrogen(space);
		potential.setBox(box);
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		
		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster,getRandom(),space);
		move.setBox(box);
		move.setPotential(potential);
		((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);
		
		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);
		((MCMoveStepTracker)rotate.getTracker()).setNoisyAdjustment(true);
		
		MCMoveVolume mcMoveVolume = new MCMoveVolume(this, potentialMaster, space);
		mcMoveVolume.setBox(box);
		mcMoveVolume.setPressure(0.0);
		
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
		
		int nCell = 2;
		double temperature = 35.0; // in Unit Kelvin
		long simSteps = 10000000;
		
		SimulationNitrogenModel sim = new SimulationNitrogenModel(Space3D.getInstance(3), nCell, temperature);
	    int numMolecule = sim.box.getMoleculeList().getMoleculeCount();
	/*
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy (per molecule): "+ Kelvin.UNIT.fromSim(latticeEnergy)/numMolecule);
		
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyListener);
						
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
		double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
		double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);
		
		System.out.println("Average energy (per molecule): "   + Kelvin.UNIT.fromSim(averageEnergy)/numMolecule  + " ;error: " + Kelvin.UNIT.fromSim(errorEnergy)/numMolecule);
		*/
		
		
		SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
	    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
	    simGraphic.makeAndDisplayFrame("Test");
		//sim.activityIntegrate.setMaxSteps(simSteps);
		//sim.getController().actionPerformed();
		
	}

	
	protected Box box;
	protected ISpace space;
	protected PotentialMaster potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected PotentialMolecular potential;
	private static final long serialVersionUID = 1L;
}
