/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.PDBWriter;
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
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncatedShifted;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesSpheresMono;

/**
 * MC simulation of FCC soft-sphere model in 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 * 
 * with SuperBox
 * 
 * @author Tai Boon Tan
 */
public class SimCalcSSoftSphereFCCSuperBox extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Boundary boundary;
    public Primitive primitive;
    public Basis basis;
    public int[] nCells;
    public SpeciesSpheresMono speciesA, speciesB;
    public CoordinateDefinitionLeafSuperBox coordinateDefinition;
    protected P2SoftSphericalTruncatedShifted pTruncatedAA;
    protected P2SoftSphericalTruncatedShifted pTruncatedAB;
    protected PotentialMasterMonatomic potentialMaster;
    public SimCalcSSoftSphereFCCSuperBox(Space _space, int numAtoms, double density, double temperature, int exponent) {
        super(_space);


        potentialMaster = new PotentialMasterMonatomic(this);

        speciesA = new SpeciesSpheresMono(this, space);
        speciesB = new SpeciesSpheresMono(this, space);
        addSpecies(speciesA);
        addSpecies(speciesB);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(speciesA, numAtoms/8);
        box.setNMolecules(speciesB, numAtoms*7/8);
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0/density);
            boundary = new BoundaryRectangularPeriodic(space, numAtoms/density);
            nCells = new int[]{numAtoms};
            basis = new BasisMonatomic(space);
        } else {
            double L = Math.pow(4.0/density, 1.0/3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int)Math.round(Math.pow(numAtoms/4, 1.0/3.0));
            nCells = new int[]{n,n,n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicFcc();
        }

        Potential2SoftSpherical potentialAA = new P2SoftSphere(space, 1.0, 1.0, exponent);
        Potential2SoftSpherical potentialAB = new P2SoftSphere(space, 1.0, 0.5, exponent);
        //Potential2SoftSpherical potentialABafter = new P2SoftSphere(space, 1.0, 1.0, exponent);

        double truncationRadius = boundary.getBoxSize().getX(0)* 0.495;
        pTruncatedAA = new P2SoftSphericalTruncatedShifted(space, potentialAA, truncationRadius);
        pTruncatedAB = new P2SoftSphericalTruncatedShifted(space, potentialAB, truncationRadius);
        //pTruncatedABafter = new P2SoftSphericalTruncatedShifted(space, potentialABafter, truncationRadius);
        potentialMaster.lrcMaster().setEnabled(false); //turn off the long-range correction ::updated 7/4/2008

        AtomType sphereTypeA = speciesA.getLeafType();
        AtomType sphereTypeB = speciesB.getLeafType();
        potentialMaster.addPotential(pTruncatedAA, new AtomType[]{sphereTypeA, sphereTypeA});
        potentialMaster.addPotential(pTruncatedAB, new AtomType[]{sphereTypeA, sphereTypeB});

        box.setBoundary(boundary);
        coordinateDefinition = new CoordinateDefinitionLeafSuperBox(box, primitive, basis, space);
        coordinateDefinition.setSpecies(speciesA, speciesB);
        coordinateDefinition.setIs256();
        coordinateDefinition.initializeCoordinates(nCells);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);
        MCMoveAtomSuperBox move = new MCMoveAtomSuperBox(potentialMaster, getRandom(), space, coordinateDefinition);
        move.setStepSize(0.2);
        move.setStepSizeMax(0.5);
        integrator.getMoveManager().addMCMove(move);
        ((MCMoveStepTracker)move.getTracker()).setNoisyAdjustment(true);

        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        move.setPotential(pTruncatedAA);

        // activityIntegrate.setMaxSteps(nSteps);

        /*
         * 1-body Potential to Constraint the atom from moving too far
         * 	away from its lattice-site
         */
       P1Constraint p1Constraint = new P1Constraint(space, primitive.getSize()[0], box, coordinateDefinition);
        potentialMaster.addPotential(p1Constraint, new AtomType[]{sphereTypeA});

       integrator.setBox(box);
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 1.256;
        double temperature = 0.2;
        int exponent = 12;
        if (D == 1) {
            nA = 3;
            density = 1.0;
        }
        long simSteps =1000000;

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
        String filename = "Super_CB_FCC_n" + exponent + "_T"+temperature;
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " soft sphere superbox simulation");
        System.out.println(nA + " atoms with exponent " + exponent+" and density "+density);
        System.out.println("isotherm temperature at "+temperature);
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcSSoftSphereFCCSuperBox sim = new SimCalcSSoftSphereFCCSuperBox(Space.getInstance(D), nA*8, density, temperature, exponent);

        // set up initial configuration and save nominal positions
        Primitive primitive = sim.primitive;

        // set up normal-mode meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactory waveVectorFactory;
        if (D == 1) {
            waveVectorFactory = new WaveVectorFactory1D();
        } else if (D == 2) {
            waveVectorFactory = null;
        } else {
            waveVectorFactory = new WaveVectorFactorySuperBox(primitive, sim.space);
        }
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        IntegratorListenerAction meterNormalModeListener = new IntegratorListenerAction(meterNormalMode);
        meterNormalModeListener.setInterval(nA);
        sim.integrator.getEventManager().addListener(meterNormalModeListener);

        //MeterPressureSuperBox meterPressure = new MeterPressureSuperBox(sim.space);
        //meterPressure.setIntegrator(sim.integrator);
        //System.out.println("\nPressure Lattice: "+ meterPressure.getDataAsScalar());
        //System.exit(1);


        MeterPotentialEnergy meterEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterEnergy.setBox(sim.box);
        System.out.println("Lattice Energy per particle: "+ meterEnergy.getDataAsScalar()/nA);
        System.out.println(" ");
//        ((P2SoftSphere)pTruncatedABbefore.getWrappedPotential()).setEpsilon(1);
//        sim.potentialMaster.removePotential(pTruncatedABbefore);
//        sim.potentialMaster.addPotential(pTruncatedABafter, new IAtomType[] {sim.speciesA.getLeafType(), sim.speciesB.getLeafType()});
        /*
        AccumulatorAverage pressureAverage = new AccumulatorAverageCollapsing();
	    DataPump pressurePump = new DataPump(meterPressure, pressureAverage);

	    sim.integrator.addIntervalAction(pressurePump);
	    sim.integrator.setActionInterval(pressurePump, 100);
	    */
        AccumulatorAverage energyAverage = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(meterEnergy, energyAverage);

        IntegratorListenerAction energyPumpListener = new IntegratorListenerAction(energyPump);
        energyPumpListener.setInterval(100);
        sim.integrator.getEventManager().addListener(energyPumpListener);


        /*
        SimulationGraphic simGraphic = new SimulationGraphic(sim, Space.getInstance(3),sim.getController());
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(50));

        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.speciesA.getLeafType(), java.awt.Color.BLACK);
        colorScheme.setColor(sim.speciesB.getLeafType(), java.awt.Color.WHITE);

        simGraphic.makeAndDisplayFrame("Sim CalcS Super Box");
        */

        sim.activityIntegrate.setMaxSteps(simSteps/10);  //simSteps/10
        sim.getController().actionPerformed();
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
        sWriterListener.setInterval((int)simSteps/20);
        sim.integrator.getEventManager().addListener(sWriterListener);

        sim.activityIntegrate.setMaxSteps(simSteps);
        sim.getController().actionPerformed();
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("calcS_n"+exponent+"_T"+temperature+".pdb");
        pdbWriter.actionPerformed();

        System.out.println("Average Energy: " + energyAverage.getData().getValue(AccumulatorAverage.AVERAGE.index) / nA);
        System.out.println("Error Energy: " + energyAverage.getData().getValue(AccumulatorAverage.ERROR.index) / nA);
        System.out.println(" ");
        /*
        System.out.println("Average Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.AVERAGE.index));
        System.out.println("Error Pressure: "+ ((DataGroup)pressureAverage.getData()).getValue(AccumulatorAverage.StatType.ERROR.index));
    	*/
    }
}
