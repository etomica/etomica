/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressureMolecular;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisHcp;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveHexagonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.nbr.list.molecule.BoxAgentSourceCellManagerListMolecular;
import etomica.nbr.list.molecule.NeighborListManagerSlantyMolecular;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Degree;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import java.io.FileWriter;
import java.io.IOException;



/**
 * Simulation class for nitrogen molecules
 * beta-N2 crystal Structure
 * 
 * @author Tai Boon Tan
 *
 */
public class SimulationBetaNitrogenModel extends Simulation{

	
	public SimulationBetaNitrogenModel(Space space, int numMolecule, double temperature, double density) {
        super(space);
        this.space = space;

        BoxAgentSourceCellManagerListMolecular boxAgentSource = new BoxAgentSourceCellManagerListMolecular(this, null, space);
        BoxAgentManager<NeighborCellManagerMolecular> boxAgentManager = new BoxAgentManager<NeighborCellManagerMolecular>(boxAgentSource, this);


        double ratio = 1.631;
        double a = Math.pow(4.0 / (Math.sqrt(3.0) * ratio * density), 1.0 / 3.0);
        double c = a * ratio;
        System.out.println("\na: " + a + " ;cDim: " + c);

        int nC = (int) Math.pow(numMolecule / 1.999999999, 1.0 / 3.0);

        Basis basisHCP = new BasisHcp();
        Basis basis = new BasisBigCell(space, basisHCP, new int[]{nC, nC, nC});

        species = new SpeciesN2(space);
        addSpecies(species);

        box = this.makeBox();
        box.setNMolecules(species, numMolecule);
        int[] nCells = new int[]{1, 1, 1};

        Vector[] boxDim = new Vector[3];
        boxDim[0] = space.makeVector(new double[]{nC * a, 0, 0});
        boxDim[1] = space.makeVector(new double[]{-nC * a * Math.cos(Degree.UNIT.toSim(60)), nC * a * Math.sin(Degree.UNIT.toSim(60)), 0});
        boxDim[2] = space.makeVector(new double[]{0, 0, nC * c});

        Boundary boundary = new BoundaryDeformablePeriodic(space, boxDim);
        primitive = new PrimitiveHexagonal(space, (nC) * a, nC * c);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsBeta();
        coordinateDef.setOrientationVectorBeta(space);
        coordinateDef.initializeCoordinates(nCells);

        box.setBoundary(boundary);
        double rC = a * nC * 0.475;
        System.out.println("Truncation Radius: " + rC);
        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        double[] u = new double[20];
        if (true) {
//			BetaPhaseLatticeParameter parameters = new BetaPhaseLatticeParameter();
//			double[][] param = parameters.getParameter(density);

            BetaPhaseLatticeParameterNA parameters = new BetaPhaseLatticeParameterNA();
            double[][] param = parameters.getParameter(numMolecule);

//			BetaPhaseLatticeParameterLS parameters = new BetaPhaseLatticeParameterLS();
//			double[][] param = parameters.getParameter(density);

            int kParam = 0;
            for (int i = 0; i < param.length; i++) {
                for (int j = 0; j < param[0].length; j++) {
                    u[kParam] = param[i][j];
                    kParam++;
                }
            }

            int numDOF = coordinateDef.getCoordinateDim();
            double[] newU = new double[numDOF];

            if (true) {
                for (int j = 0; j < numDOF; j += 10) {
                    if (j > 0 && j % (nC * 10) == 0) {
                        j += nC * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k];
                    }
                }

                for (int j = nC * 10; j < numDOF; j += 10) {
                    if (j > nC * 10 && j % (nC * 10) == 0) {
                        j += nC * 10;
                        if (j >= numDOF) {
                            break;
                        }
                    }
                    for (int k = 0; k < 10; k++) {
                        newU[j + k] = u[k + 10];
                    }
                }
            }

            coordinateDef.setToU(box.getMoleculeList(), newU);
            coordinateDef.initNominalU(box.getMoleculeList());

        }

//		PRotConstraint pRot= new PRotConstraint(space, coordinateDef, box);
//		pRot.setBox(box);
//		pRot.setConstraintAngle(0.1);

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
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        int numNeigh = potentialMaster.getNeighborManager(box).getUpList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
        System.out.println("numNeigh: " + numNeigh);

        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        move.setBox(box);
        move.setPotential(potential);
        move.setDoExcludeNonNeighbors(true);

        MCMoveRotateMolecule3D rotate = new MCMoveRotateMolecule3D(potentialMaster, getRandom(), space);
        rotate.setBox(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature), box);
        integrator.getMoveManager().addMCMove(move);
        integrator.getMoveManager().addMCMove(rotate);

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


//		MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
//		meterPE.setBox(box);
//		double latticeEnergy = meterPE.getDataAsScalar();
//		double[][] eigenvector = ArrayReader1D.getFromFile("/tmp/NoFindPair432eVecALL");
//		
//		int numComp = eigenvector[0].length;
//
////		for (int i=0; i<numComp; i++){
////			System.out.println(i+" "+ eigenvector[1][i]);
////		}
////		System.exit(1);
//		
//		double interval = 0.01;
//		double[] dev = new double[numComp];
//		for (int i=-50; i<=50; i++){
//			for(int y=0; y<numComp; y++){
//				dev[y] = i*interval*eigenvector[0][y];
//			}
//			
//			coordinateDef.setToU(box.getMoleculeList(), dev);
//			double pe = meterPE.getDataAsScalar();
//			System.out.println((i*interval)+" "+ (pe-latticeEnergy)/numMolecule);
//		}
//		
//		System.exit(1);


        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
    }
	
	public static void main (String[] args){
		
		double temperature =0.000001; // in Unit Kelvin
		long simSteps = 100000;
		double density = 0.023;
		int nC = 6;
		int numMolecule = nC*nC*nC*2;
		if(args.length > 0){
			simSteps = Long.parseLong(args[0]);
		}
		if(args.length > 1){
			temperature = Double.parseDouble(args[1]);
		}
		if(args.length > 2){
			numMolecule = Integer.parseInt(args[2]);
		}
		if(args.length > 3){
			density = Double.parseDouble(args[3]);
		}
		
		System.out.println("Running beta-N2 crystal structure simulation with " + simSteps + " steps" );
		System.out.println("num Molecules: " + numMolecule+ " ; temperature: " + temperature+"K ");
		System.out.println("density: " + density);
		
		
		SimulationBetaNitrogenModel sim = new SimulationBetaNitrogenModel(Space3D.getInstance(3), numMolecule, temperature, density);
	    
		final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
		final double latticeEnergy = meterPotentialEnergy.getDataAsScalar();
		System.out.println("Lattice Energy per molecule (sim unit): "+ latticeEnergy/numMolecule);
		System.out.println("Lattice Energy: "+ latticeEnergy);
//		System.exit(1);
		MeterPressureMolecular meterPressure = new MeterPressureMolecular(sim.space);
		meterPressure.setIntegrator(sim.integrator);
			
		double staticPressure = meterPressure.getDataAsScalar();
		System.out.println("Static Pressure (sim unit): " + staticPressure);
		
		double volume = sim.box.getBoundary().volume();
		System.out.println("volume: " + volume);
		
		if(true){
			SimulationGraphic simGraphic = new SimulationGraphic(sim);
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
		    simGraphic.makeAndDisplayFrame("Beta-Phase Nitrogen Crystal Structure");
		    
		    DiameterHashByType diameter = new DiameterHashByType();
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

		sim.getController().reset();

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
		
		sim.activityIntegrate.setMaxSteps(simSteps);
		sim.getController().actionPerformed();
		
//		sim.writeUdistribution(filename, meterOrient);
		
		double averageEnergy = energyAverage.getData().getValue(energyAverage.AVERAGE.index);
		double errorEnergy = energyAverage.getData().getValue(energyAverage.ERROR.index);
		
		double averagePressure = pressureAverage.getData().getValue(energyAverage.AVERAGE.index);
		double errorPressure = pressureAverage.getData().getValue(energyAverage.ERROR.index);
		
		System.out.println("Average energy (per molecule): "   + averageEnergy/numMolecule  
				+ " ;error: " + errorEnergy/numMolecule);
		System.out.println("Average pressure (sim unit): " + averagePressure 
				+ " ;error: " + errorPressure);

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
	protected Space space;
	protected PotentialMasterListMolecular potentialMaster;
	protected IntegratorMC integrator;
	protected ActivityIntegrate activityIntegrate;
	protected P2Nitrogen potential;
	protected CoordinateDefinitionNitrogen coordinateDef;
	protected Primitive primitive;
	protected SpeciesN2 species;
	private static final long serialVersionUID = 1L;
}
