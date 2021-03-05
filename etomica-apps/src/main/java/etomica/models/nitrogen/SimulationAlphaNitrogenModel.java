/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureMolecular;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import java.io.FileWriter;
import java.io.IOException;



/**
 * Simulation class for nitrogen molecules
 * alpha-N2 crystal Structure
 * 
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationAlphaNitrogenModel extends Simulation{

	
	public SimulationAlphaNitrogenModel(Space space, int[] nC, double temperature, double density) {
		super(space);
		this.space = space;

		double a = Math.pow(4.0 / density, 1.0 / 3.0);
		System.out.println("Unit Cell Length, a: " + a);

		Basis basisFCC = new BasisCubicFcc();
		Basis basis = new BasisBigCell(space, basisFCC, new int[]{nC[0], nC[1], nC[2]});

		species = SpeciesN2.create(false);
		addSpecies(species);

        potentialMaster = new PotentialMasterListMolecular(this, space);

		int numMolecule = 4 * nC[0] * nC[1] * nC[2];
		double[] boxSize = new double[]{nC[0] * a, nC[1] * a, nC[2] * a};
		Boundary boundary = new BoundaryRectangularPeriodic(space, boxSize);
		box = this.makeBox(boundary);
		box.setNMolecules(species, numMolecule);

		int[] nCells = new int[]{1, 1, 1};

		primitive = new PrimitiveTetragonal(space, nC[0] * a, nC[2] * a);

		coordinateDef = new CoordinateDefinitionNitrogen(getSpeciesManager(), box, primitive, basis, space);
		coordinateDef.setIsAlpha();
		coordinateDef.setOrientationVectorAlpha(space);
		coordinateDef.initializeCoordinates(nCells);
		double rCScale = 0.475;
		double rC = box.getBoundary().getBoxSize().getX(0) * rCScale;
		System.out.println("Truncation Radius (" + rCScale + " Box Length): " + rC);
		potential = new P2Nitrogen(space, rC);
		potential.setBox(box);

		PRotConstraint pRotConstraint = new PRotConstraint(space, coordinateDef, box);
		pRotConstraint.setConstraintAngle(65);

		potentialMaster.addPotential(potential, new ISpecies[]{species, species});
//		potentialMaster.addPotential(pRotConstraint,new ISpecies[]{species} );
		//potentialMaster.lrcMaster().isEnabled();

		int cellRange = 6;
		potentialMaster.setRange(rC);
		potentialMaster.setCellRange(cellRange);
		potentialMaster.getNeighborManager(box).reset();

		int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
		if (potentialCells < cellRange * 2 + 1) {
			throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
		}
		potential.setRange(Double.POSITIVE_INFINITY);

		int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().get(0))[0].size();
		System.out.println("numNeigh: " + numNeigh);

		MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
		move.setBox(box);
		move.setPotential(potential);
		move.setDoExcludeNonNeighbors(true);

		MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
		rotate.setBox(box);

//		((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);	
//		((MCMoveStepTracker)rotate.getTracker()).setNoisyAdjustment(true);

		integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature), box);
		integrator.getMoveManager().addMCMove(move);
		integrator.getMoveManager().addMCMove(rotate);

//		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
//		meterPE.setBox(box);
//		System.out.println("lattice energy (sim unit): " + meterPE.getDataAsScalar()/numMolecule);
//		//System.exit(1);

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

		this.getController().addActivity(new ActivityIntegrate(integrator));
	}
	
	public static void main (String[] args){
		int nCx = 6;
		int nCy = 6;
		int nCz = 6;
		
		double temperature = 0.002; // in Unit Kelvin
		long simSteps = 100000;
		double density = 0.0230; //0.02204857502170207 (intial from literature with a = 5.661)
		
		if(args.length > 0){
			simSteps = Long.parseLong(args[0]);
		}
		if(args.length > 1){
			temperature = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			nCx = Integer.parseInt(args[2]);
		}
		if(args.length > 3){
			nCy = Integer.parseInt(args[3]);
		}
		if(args.length > 4){
			nCz = Integer.parseInt(args[4]);
		}
		if(args.length > 5){
			density = Double.parseDouble(args[5]);
		}
		
		
		int[] nC = new int[]{nCx, nCy, nCz};
		int numMolecule = 4*nC[0]*nC[1]*nC[2];
		
		String filename = "alphaN2_nA"+numMolecule+"_T"+temperature;
		System.out.println("Running alpha-N2 NVT simulation with " + simSteps + " steps" );
		System.out.println("number Molecules: " + numMolecule+ " ; temperature: " + temperature
				+"K ; density: "+ density +"\n");

		SimulationAlphaNitrogenModel sim = new SimulationAlphaNitrogenModel(Space3D.getInstance(3), nC, temperature, density);
	 
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
		final double latticeEnergySim = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy per molecule (sim unit): "+ latticeEnergySim/numMolecule);
		System.out.println("Lattice Energy: "+ latticeEnergySim);			
		
//		MeterNormalizedCoord meterCoord = new MeterNormalizedCoord(sim.box, sim.coordinateDef, sim.species);
//        IntegratorListenerAction meterCoordListener = new IntegratorListenerAction(meterCoord);
//        meterCoordListener.setInterval(numMolecule);                                      
//        sim.integrator.getEventManager().addListener(meterCoordListener);       
//		
//		MeterOrientationDistribution meterOrient = new MeterOrientationDistribution(sim.box, sim.coordinateDef, sim.species);
//        IntegratorListenerAction meterOrientListener = new IntegratorListenerAction(meterOrient);
//        meterOrientListener.setInterval(numMolecule);                                      
//        sim.integrator.getEventManager().addListener(meterOrientListener);
		
		MeterPressureMolecular meterPressure = new MeterPressureMolecular(sim.space);
		meterPressure.setIntegrator(sim.integrator);
						
		double staticPressure = meterPressure.getDataAsScalar();
		System.out.println("Static Pressure (sim unit): " + staticPressure);
		
		double volume = sim.box.getBoundary().volume();
		System.out.println("volume: " + volume);
		
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim);
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    	    
			DiameterHashByType diameter = new DiameterHashByType();
			diameter.setDiameter(sim.species.getTypeByName("N"), 3.1);
			diameter.setDiameter(sim.species.getTypeByName("P"), 0.0);
			simGraphic.getDisplayBox(sim.box).setDiameterHash(diameter);
			
		    simGraphic.makeAndDisplayFrame("Alpha-Phase Nitrogen Crystal Structure");
		    
			
			IAction output = new IAction(){
				public void actionPerformed() {
					System.out.println("energy: " + (meterPotentialEnergy.getDataAsScalar()-latticeEnergySim));
					
				}
				
			};
			
			IntegratorListenerAction outListener = new IntegratorListenerAction(output);
			outListener.setInterval(100);
			sim.integrator.getEventManager().addListener(outListener);
			
//			double averageEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index);
//			double errorEnergy = ((DataGroup)energyAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index);

            sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, 1000000));
		    return;
		    
		}

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, simSteps / 5));
System.out.println("****System Equilibrated (20% of SimSteps)****");

		long startTime = System.currentTimeMillis();
		System.out.println("\nStart Time: " + startTime);



		AccumulatorAverage energyAverage = new AccumulatorAverageFixed();
		DataPump energyPump = new DataPump(meterPotentialEnergy, energyAverage);
		IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
		energyListener.setInterval(numMolecule);
		sim.integrator.getEventManager().addListener(energyListener);

		AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
		DataPump pressurePump = new DataPump(meterPressure, pressureAverage);
		IntegratorListenerAction pressureListener = new IntegratorListenerAction(pressurePump);
		pressureListener.setInterval((int)simSteps/200);
		sim.integrator.getEventManager().addListener(pressureListener);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, simSteps));

		double averageEnergy = energyAverage.getData().getValue(energyAverage.AVERAGE.index);
		double errorEnergy = energyAverage.getData().getValue(energyAverage.ERROR.index);

		double averagePressure = pressureAverage.getData().getValue(energyAverage.AVERAGE.index);
		double errorPressure = pressureAverage.getData().getValue(energyAverage.ERROR.index);
		
		System.out.println("Average energy (per molecule): " + (averageEnergy)/numMolecule  
				+ " ;error: " + errorEnergy/numMolecule);
		System.out.println("Average pressure (sim unit): " + averagePressure
				+ " ;error: " + errorPressure);
		
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
	protected Space space;
	protected PotentialMasterListMolecular potentialMaster;
	protected IntegratorMC integrator;
	
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesGeneral species;
	private static final long serialVersionUID = 1L;
}
