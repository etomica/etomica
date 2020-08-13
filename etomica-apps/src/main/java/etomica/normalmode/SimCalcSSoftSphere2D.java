/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.PDBWriter;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPressure;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveOrthorhombicHexagonal;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of FCC soft-sphere model in 2D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * @author Tai Boon Tan
 */
public class SimCalcSSoftSphere2D extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public Primitive primitive;
    public Basis basis;
    public int[] nCells;
    public CoordinateDefinition coordinateDefinition;
    public PotentialMasterMonatomic potentialMaster;
    public SimCalcSSoftSphere2D(Space _space, int numAtoms, int[] nCells, double temperature, int exponent) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        potentialMaster = new PotentialMasterMonatomic(this);

        primitive = new PrimitiveOrthorhombicHexagonal(space, 1);
        Vector[] dimension = space.makeVectorArray(2);
        for (int i = 0; i < space.D(); i++) {
            dimension[i].Ea1Tv1(nCells[i], primitive.vectors()[i]);
        }
        boundary = new BoundaryDeformablePeriodic(space, dimension);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MCMoveAtomCoupled move = new MCMoveAtomCoupled(potentialMaster, new MeterPotentialEnergy(potentialMaster), getRandom(), space);
        move.setStepSize(0.2);
        move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker) move.getTracker()).setNoisyAdjustment(true);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);
        // activityIntegrate.setMaxSteps(nSteps);

        basis = new BasisOrthorhombicHexagonal();

        Potential2SoftSpherical potential = new P2SoftSphere(space);

        double truncationRadius = boundary.getBoxSize().getX(0) * 0.495;
        P2SoftSphericalTruncatedShifted pTruncated = new P2SoftSphericalTruncatedShifted(space, potential, truncationRadius);
        //potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008

        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});
        move.setPotential(pTruncated);

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(nCells);


        /*
         * 1-body Potential to Constraint the atom from moving too far
         * 	away from its lattice-site
         */
        P1Constraint p1Constraint = new P1Constraint(space, primitive.getSize()[0], box, coordinateDefinition);
        potentialMaster.addPotential(p1Constraint, new AtomType[]{sphereType});
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
    	int D = 2;
        int dim2D = 2;
        int[] nCells = new int[] {2*dim2D, dim2D} ;
        int nA = 2*nCells[0]*nCells[1];
        double temperature = 0.1;
        int exponent = 12 ;

        long simSteps =100000;

        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            nA = Integer.parseInt(args[2]);
        }
        if (args.length > 3) {
            temperature = Double.parseDouble(args[3]);
        }
        if (args.length > 4) {
        	exponent = Integer.parseInt(args[4]);
        }
        String filename = "Test_2D_CB_FCC_n"+exponent+"_T0"+ (int)Math.round(temperature*10);
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " soft sphere simulation");
        System.out.println(nA + " atoms with exponent " + exponent);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSSoftSphere2D sim = new SimCalcSSoftSphere2D(Space.getInstance(D), nA, nCells, temperature, exponent);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.primitive;

        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory= new WaveVectorFactory2D(primitive, sim.space);

        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        IntegratorListenerAction meterNormalModeListener = new IntegratorListenerAction(meterNormalMode);
        meterNormalModeListener.setInterval(nA);
        sim.integrator.getEventManager().addListener(meterNormalModeListener);


        MeterPressure meterPressure = new MeterPressure(sim.space);
        meterPressure.setIntegrator(sim.integrator);
        System.out.println("\nPressure Lattice: "+ meterPressure.getDataAsScalar());
        //System.exit(1);


        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        System.out.println("Lattice Energy per particle: "+ meterEnergy.getDataAsScalar()/nA);
        System.out.println(" ");
        //System.exit(1);

        AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);

        IntegratorListenerAction pressurePumpListener = new IntegratorListenerAction(pressurePump);
        pressurePumpListener.setInterval(100);
	    sim.integrator.getEventManager().addListener(pressurePumpListener);

        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);

        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        energyPumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(energyPumpListener);

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), simSteps/10);  //simSteps/10
System.out.println("equilibrated");
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
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), simSteps);
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("calcS_nA"+nA+"_n"+exponent+"_T"+temperature+".pdb");
        pdbWriter.actionPerformed();

        System.out.println("Average Energy: " + energyAverage.getData().getValue(energyAverage.AVERAGE.index));
        System.out.println("Error Energy: " + energyAverage.getData().getValue(energyAverage.ERROR.index));
        System.out.println(" ");

        System.out.println("Average Pressure: " + pressureAverage.getData().getValue(pressureAverage.AVERAGE.index));
        System.out.println("Error Pressure: " + pressureAverage.getData().getValue(pressureAverage.ERROR.index));
    }
}
