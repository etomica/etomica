/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.*;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimCalcSLJ extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public Primitive primitive, primitiveUnitCell;
    public Basis basis;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMaster potentialMaster;
    protected P1Constraint p1Constraint;
    public SimCalcSLJ(Space _space, int numAtoms, double density, double temperature) {
        super(_space);

        potentialMaster = new PotentialMasterMonatomic(this);
        potentialMaster.lrcMaster().setEnabled(false);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature);
        MCMoveAtomCoupled move = new MCMoveAtomCoupled(potentialMaster, new MeterPotentialEnergy(potentialMaster), getRandom(), space);
        move.setStepSize(0.1);
        move.setStepSizeMax(0.5);
        move.setDoExcludeNonNeighbors(true);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            primitive = new PrimitiveCubic(space, n*L);
            primitiveUnitCell = new PrimitiveCubic(space, L);

            nCells = new int[]{n,n,n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }

        Potential2SoftSpherical potential = new P2LennardJones(space, 1.0, 1.0);
        double truncationRadius = boundary.getBoxSize().getX(0) * 0.45;
        if (potentialMaster instanceof PotentialMasterList) {
            // with neighborlisting, we can use an unshifted potential
            potential = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        }
        else {
            // without neighborlisting, we should use a shifted potential
            // it should provide a better approximation of the actual variation in
            // potential energy as atoms move in and out of the cutoff
            potential = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        }
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});
        move.setPotential(potential);

        box.setBoundary(boundary);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1,1,1});


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
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 0.962;
        double temperature = 0.1378;
        if (D == 1) {
            nA = 3;
            density = 0.5;
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

        String filename = "LJCB_FCC_nA"+nA+ "_d"+density
		+"T"+ (int) Math.round(temperature);
		if (args.length > 0) {
			filename = args[0];
		}

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " Lennard-Jones simulation");
        System.out.println(nA + " atoms at density " + density+" and temperature "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSLJ sim = new SimCalcSLJ(Space.getInstance(D), nA, density, temperature);

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

        sim.activityIntegrate.setMaxSteps(simSteps/10);
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
