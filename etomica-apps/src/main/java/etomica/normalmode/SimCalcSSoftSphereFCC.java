/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeRandom;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Pixel;

/**
 * MC simulation of FCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * @author Tai Boon Tan
 */
public class SimCalcSSoftSphereFCC extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public Primitive primitive, primitiveUnitCell;
    public Basis basis;
    public int[] nCells;
    public P1Constraint p1Constraint;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMaster potentialMaster;

    public SimCalcSSoftSphereFCC(Space _space, int numAtoms, double density, double temperature, int exponent) {
		super(_space);

		potentialMaster = new PotentialMasterList(this, space);
		//potentialMaster = new PotentialMasterMonatomic(this);

		SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
		addSpecies(species);

		box = new Box(space);
		addBox(box);
		box.setNMolecules(species, numAtoms);

		integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);

		if (space.D() == 1) {
			primitive = new PrimitiveCubic(space, 1.0 / density);
			boundary = new BoundaryRectangularPeriodic(space, numAtoms/ density);
			nCells = new int[] { numAtoms };
			basis = new BasisMonatomic(space);
		} else {
			double L = Math.pow(4.0 / density, 1.0 / 3.0);
			int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
			primitive = new PrimitiveCubic(space, n * L);
			primitiveUnitCell = new PrimitiveCubic(space, L);

			nCells = new int[] {n, n, n};
			boundary = new BoundaryRectangularPeriodic(space, n * L);
			Basis basisFCC = new BasisCubicFcc();
			basis = new BasisBigCell(space, basisFCC, nCells);
		}

		Potential2SoftSpherical potential = new P2SoftSphere(space);
		double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;

		/*
		 * When we consider the interaction of the neighborlist, we are safe to use <P2SoftSphereTruncated>;
		 * if not, we will have to use <P2SoftSphereTruncatedShifted>
		 * This is because the possible occurrence of atom "jumping" in and out of the truncation radius
		 * 	might potentially causes our energy average to be off
		 */
		if(potentialMaster instanceof PotentialMasterList){
			potential = new P2SoftSphericalTruncated(space, potential, truncationRadius);

        } else {
			potential = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);

        }
		/*
		 * this is to make sure we don't include the long-tailed correction to
		 * our solid calculation; and it is only important for fluid simulation
		 */
        potentialMaster.lrcMaster().setEnabled(false);

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});
        box.setBoundary(boundary);


        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
		coordinateDefinition.initializeCoordinates(new int[]{1,1,1});

        MCMoveAtomCoupled move = new MCMoveAtomCoupled(potentialMaster, new MeterPotentialEnergy(potentialMaster), getRandom(), space);
        move.setPotential(potential);
		move.setDoExcludeNonNeighbors(true);

        integrator.getMoveManager().addMCMove(move);
		((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);

        activityIntegrate = new ActivityIntegrate(integrator);

		/*
		 * 1-body Potential to Constraint the atom from moving too far away from
		 * its lattice-site
		 */
		p1Constraint = new P1Constraint(space, primitiveUnitCell.getSize()[0], box, coordinateDefinition);
        potentialMaster.addPotential(p1Constraint, new AtomType[]{sphereType});

		if (potentialMaster instanceof PotentialMasterList) {
            double neighborRange = truncationRadius;
            int cellRange = 7;
            ((PotentialMasterList)potentialMaster).setRange(neighborRange);
            ((PotentialMasterList)potentialMaster).setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
            // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
            ((PotentialMasterList)potentialMaster).getNeighborManager(box).reset();
            int potentialCells = ((PotentialMasterList)potentialMaster).getNbrCellManager(box).getLattice().getSize()[0];
            if (potentialCells < cellRange*2+1) {
                throw new RuntimeException("oops ("+potentialCells+" < "+(cellRange*2+1)+")");
            }
            if (potentialCells > cellRange*2+1) {
                System.out.println("could probably use a larger truncation radius ("+potentialCells+" > "+(cellRange*2+1)+")");
            }
            ((P2SoftSphericalTruncated)potential).setTruncationRadius(0.6*boundary.getBoxSize().getX(0));
		}

        integrator.setBox(box);
		getController().addAction(activityIntegrate);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// defaults
		int D = 3;
		int nA = 32;
		double density = 1.256;
		double temperature =1.4;
		int exponent = 12;
		if (D == 1) {
			nA = 3;
			density = 1.0;
		}
		long simSteps = 1000000;

		// parse arguments
		if (args.length > 1) {
			density = Double.parseDouble(args[1]);
		}
		if (args.length > 2) {
			simSteps = Long.parseLong(args[2]);
		}
		if (args.length > 3) {
			nA = Integer.parseInt(args[3]);
		}
		if (args.length > 4) {
			temperature = Double.parseDouble(args[4]);
		}
		if (args.length > 5) {
			exponent = Integer.parseInt(args[5]);
		}
		String filename = "Neigh_n" + exponent + "_T14";//+ (int) Math.round(temperature);
		if (args.length > 0) {
			filename = args[0];
		}

		System.out.println("Running " + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
				+ " soft sphere simulation");
		System.out.println(nA + " atoms with exponent " + exponent
				+ " and density " + density);
		System.out.println("isotherm temperature at " + temperature);
		System.out.println(simSteps + " steps");
		System.out.println("output data to " + filename);


        // construct simulation
		SimCalcSSoftSphereFCC sim = new SimCalcSSoftSphereFCC(Space.getInstance(D), nA, density, temperature, exponent);

		/*
         * Graphical Simulation
		 */
		if(false){
			SimulationGraphic simGraphic = new SimulationGraphic(sim);
			simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));
			ColorSchemeRandom colorRandom = new ColorSchemeRandom(sim.box, sim.random);
			simGraphic.getDisplayBox(sim.box).setColorScheme(colorRandom);
			simGraphic.makeAndDisplayFrame();
			return;
		}

		/*
		 * Quantify the atomic displacement from its lattice site
		 */
		if(false){
			MeterAtomicDisplacement meterDisplacement = new MeterAtomicDisplacement(sim.space, sim.coordinateDefinition);
			DataTableWriter histogramSink = new DataTableWriter();

            IntegratorListenerAction meterDisplacementListener = new IntegratorListenerAction(meterDisplacement);
			meterDisplacementListener.setInterval(100);
			sim.integrator.getEventManager().addListener(meterDisplacementListener);

            sim.activityIntegrate.setMaxSteps(simSteps);
			sim.activityIntegrate.actionPerformed();

            DataLogger dataLogger = new DataLogger();
	        dataLogger.setFileName(filename+"hist_dist");
	        dataLogger.setDataSink(histogramSink);
	        histogramSink.setIncludeHeader(false);
	        dataLogger.setAppending(false);
	        dataLogger.setWriteInterval(1);
	        dataLogger.putDataInfo(meterDisplacement.getDataInfo());
	        dataLogger.putData(meterDisplacement.getData());
	        dataLogger.closeFile();

            return;
		}

        // set up normal-mode meter
		MeterNormalMode meterNormalMode = new MeterNormalMode();
		meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
		WaveVectorFactory waveVectorFactory;
		if (D == 1) {
			waveVectorFactory = new WaveVectorFactory1D();
		} else if (D == 2) {
			waveVectorFactory = null;
		} else {
			waveVectorFactory = new WaveVectorFactorySimple(sim.primitive, sim.space);
		}
		meterNormalMode.setWaveVectorFactory(waveVectorFactory);
		meterNormalMode.setBox(sim.box);

		IntegratorListenerAction meterNormalModeListener = new IntegratorListenerAction(meterNormalMode);
		meterNormalModeListener.setInterval(nA);
		sim.integrator.getEventManager().addListener(meterNormalModeListener);

        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
		meterEnergy.setBox(sim.box);
		double latticeEnergy = meterEnergy.getDataAsScalar();
		System.out.println("Lattice Energy per particle: " + (latticeEnergy /nA));
		System.out.println(" ");

		AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
		DataPump energyPump = new DataPump(meterEnergy, energyAverage);

		IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
		energyPumpListener.setInterval(100);
		sim.integrator.getEventManager().addListener(energyPumpListener);

		sim.activityIntegrate.setMaxSteps(simSteps / 10); // simSteps/10
		sim.getController().actionPerformed();
		System.out.println("equilibrated");

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

        double A = sWriter.getLastA();
		System.out.println("A/N: " + A/nA);
        System.out.println("Average Energy: " + energyAverage.getData().getValue(AccumulatorAverage.AVERAGE.index)
                + " ,Error: " + energyAverage.getData().getValue(AccumulatorAverage.ERROR.index));
        System.out.println(" ");

        long endTime = System.currentTimeMillis();
		System.out.println("End Time: " + endTime);
		System.out.println("Time taken: " + (endTime - startTime));

	}

}
