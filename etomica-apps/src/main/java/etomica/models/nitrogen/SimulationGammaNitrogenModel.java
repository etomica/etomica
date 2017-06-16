/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
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
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
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

	
	public SimulationGammaNitrogenModel(Space space, int numMolecule, double temperature, double pressure) {
		super(space);
		this.space = space;
		double inita = 3.957;
		double initc = 5.109;
		double ratio = initc/inita; 
		int nCell = (int)Math.round(Math.pow((numMolecule/2), 1.0/3.0));
	
		double density = numMolecule/(nCell*nCell*nCell*inita*inita*initc)*1.1;
		
		//MinimizeGammaNitrogenLatticeParameter minimizer = new MinimizeGammaNitrogenLatticeParameter
		//	(space, numMolecule, density, ratio);
		
		double a = 3.8778; //minimizer.getA();
		double c = 5.3198; //minimizer.getC();
		
		potentialMaster = new PotentialMaster();
				
		Basis basisBCC = new BasisCubicBcc();
		Basis basis = new BasisBigCell(space, basisBCC, new int[]{nCell, nCell, nCell});
		
		species = new SpeciesN2ShellModel(space);
		addSpecies(species);
		
		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numMolecule);		
		int [] nCells = new int[]{1,1,1};
				
		Boundary boundary = new BoundaryRectangularPeriodic(space, new double[]{nCell*a, nCell*a, nCell*c});
		primitive = new PrimitiveTetragonal(space, nCell*a, nCell*c);
		
		double volume = nCell*a*nCell*a*nCell*c;
		System.out.println("density: " + numMolecule/volume);
		//System.exit(1);
		
		coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
		coordinateDef.setIsGamma();
		coordinateDef.setOrientationVectorGamma(space);
		coordinateDef.initializeCoordinates(nCells);
		
		double [] u = new double[coordinateDef.getCoordinateDim()];
		for (int i=3; i<coordinateDef.getCoordinateDim(); i+=5){
			u[i] = 0.00001;
			u[i+1] = 0.0;
		}
		coordinateDef.setToU(coordinateDef.getBox().getMoleculeList(), u);
		
		double[] devi = coordinateDef.calcU(coordinateDef.getBox().getMoleculeList());
		for (int i = 0; i<devi.length; i++){
			System.out.println("devi[" +i +"]: " + devi[i] );
			
		}
		
		
		System.exit(1);
		
		
		box.setBoundary(boundary);
		double rCScale = 0.45;
		double rC = box.getBoundary().getBoxSize().getX(0)*rCScale;
		System.out.println("Truncation Radius (" + rCScale +" Box Length): " + rC);
		potential = new P2NitrogenShellModel(space, rC);
		potential.setBox(box);
		
		PRotConstraint pRotConstraint = new PRotConstraint(space, coordinateDef, box);
		pRotConstraint.setConstraintAngle(90);
		
		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
		//potentialMaster.addPotential(pRotConstraint, new ISpecies[]{species});
		
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
		
		int nCell = 2;
		int numMolecule =nCell*nCell*nCell*2;
		double temperature = 10.0; // in Unit Kelvin
		double pressure = 0.21; //in Unit GPa
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
		String filename = "gammaN2_nA"+numMolecule+"_T"+Math.round(temperature);
		
		if(args.length > 0){
			filename = args[0];
		} 
		System.out.println("Running gamma-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; pressure: "+ pressure+"GPa\n");

		
		SimulationGammaNitrogenModel sim = new SimulationGammaNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, pressure);

		if(true){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
		    
			DiameterHashByType diameter = new DiameterHashByType(sim);
			diameter.setDiameter(sim.species.getNitrogenType(), 3.1);
			diameter.setDiameter(sim.species.getPType(), 0.0);
			
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
			System.out.println("Diameter is: " + diameter.getDiameter(sim.species.getPType()));
		    
				
		
		    
		    ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
		    //colorScheme.setColor(sim.species.getNitrogenType(),java.awt.Color.red);
		    simGraphic.makeAndDisplayFrame("Gamma-Phase Nitrogen Crystal Structure");
		    
			sim.activityIntegrate.setMaxSteps(simSteps);
			//sim.getController().actionPerformed();
			return;
		}
	    

	    	    
		MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterPotentialEnergy.setBox(sim.box);
		double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy (per molecule) in K: "+ Kelvin.UNIT.fromSim(latticeEnergy)/numMolecule);
		
		MeterPressure meterPressure = new MeterPressure(sim.space);
		meterPressure.setIntegrator(sim.integrator);
		
		double staticPressure = meterPressure.getDataAsScalar();
		System.out.println("Static Pressure: " + staticPressure);
		
		double volume = sim.getBox(0).getBoundary().volume();
		System.out.println("Enthaply, H: " + (latticeEnergy + staticPressure*volume)/numMolecule);
		
		//System.exit(1);
		
		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyListener);
		
		sim.activityIntegrate.setMaxSteps(simSteps/5);
		//sim.getController().actionPerformed();
		System.out.println("****System Equilibrated (20% of SimSteps)****");
		
		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);
		sim.integrator.getMoveManager().setEquilibrating(false);
		sim.getController().reset();

		sim.activityIntegrate.setMaxSteps(simSteps);
		//sim.getController().actionPerformed();

		
		double averageEnergy = energyAverage.getData().getValue(energyAverage.AVERAGE.index);
		double errorEnergy = energyAverage.getData().getValue(energyAverage.ERROR.index);
	
		System.out.println("Average energy (per molecule): "   + Kelvin.UNIT.fromSim(averageEnergy)/numMolecule  + " ;error: " + Kelvin.UNIT.fromSim(errorEnergy)/numMolecule);
	    
	    System.out.println("Box Dimension: " + sim.box.getBoundary().getBoxSize().toString());
		long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));
		

		

			
	}

	protected Box box;
	protected Space space;
	protected PotentialMaster potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected PotentialMolecular potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2ShellModel species;
	private static final long serialVersionUID = 1L;
}
