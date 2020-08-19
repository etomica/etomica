/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.DataSourceScalar;
import etomica.data.IDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.crystal.*;
import etomica.overlap.DataOverlap;
import etomica.overlap.IntegratorOverlap;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Energy;
import etomica.util.ParameterBase;
import etomica.util.ReadParameters;
import etomica.virial.overlap.AccumulatorVirialOverlapSingleAverage;
import etomica.virial.overlap.DataSourceVirialOverlap;

/**
 * Main method that constructs an overlap sampling module to calculate the free
 * energy difference between a harmonic system and a lennard jones solid.
 *
 * @author Andrew Schultz
 */
public class SimOverlapLJModule {

    private static final long serialVersionUID = 1L;
    public IntegratorOverlap integratorOverlap;
    public DataSourceVirialOverlap dsvo;
    public IntegratorBox[] integrators;

    public Box boxTarget, boxHarmonic;
    public Boundary boundaryTarget, boundaryHarmonic;
    public int[] nCells;
    public Basis basis;
    public NormalModes normalModes;
    public Primitive primitive;
    public double refPref;
    public AccumulatorVirialOverlapSingleAverage[] accumulators;
    public DataPump[] accumulatorPumps;
    public IDataSource[] meters;

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapLJModule.SimOverlapParam
     */
    public static void main(String[] args) {

        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        String inputFilename = null;
        if (args.length > 0) {
            inputFilename = args[0];
        }
        if (inputFilename != null) {
            ReadParameters readParameters = new ReadParameters(inputFilename, params);
            readParameters.readParameters();
        }
        double density = params.density;
        long numSteps = params.numSteps;
        int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        int D = params.D;
        String filename = params.filename;
        if (filename.length() == 0) {
            filename = "normal_modes_LJ_3D_" + numMolecules;
        }
        String refFileName = args.length > 0 ? filename + "_ref" : null;

        System.out.println("Running " + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal")) + " hard sphere overlap simulation");
        System.out.println(numMolecules + " atoms at density " + density + " and temperature " + temperature);
        System.out.println((numSteps / 100) + " total steps of 1000");
        System.out.println("output data to " + filename);

        //instantiate simulation
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(sim, space);
        sim.addSpecies(species);

        PotentialMaster potentialMasterTarget = new PotentialMasterMonatomic(sim);
        IntegratorBox[] integrators = new IntegratorBox[2];

        NormalModes normalModes = new NormalModesFromFile(filename, space.D());
        normalModes.setTemperature(temperature);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();

        // HARMONIC
        Boundary boundaryHarmonic = new BoundaryRectangularPeriodic(space);
        Primitive primitive;
        int[] nCells;
        Basis basis;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            nCells = new int[]{numMolecules};
            basis = new BasisMonatomic(space);
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, numMolecules / density);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            primitive = new PrimitiveCubic(space, L);
            int n = (int) Math.round(Math.pow(numMolecules / 4, 1.0 / 3.0));
            nCells = new int[]{n, n, n};
            boundaryHarmonic = new BoundaryRectangularPeriodic(space, n * L);
            basis = new BasisCubicFcc();
        }
        Box boxHarmonic = new Box(boundaryHarmonic, space);
        sim.addBox(boxHarmonic);
        boxHarmonic.setNMolecules(species, numMolecules);

        IntegratorMC integratorHarmonic = new IntegratorMC(potentialMasterTarget, sim.getRandom(), 1.0, boxHarmonic);

        MCMoveHarmonic move = new MCMoveHarmonic(sim.getRandom());
        integratorHarmonic.getMoveManager().addMCMove(move);
        integrators[0] = integratorHarmonic;


        CoordinateDefinitionLeaf coordinateDefinitionHarmonic = new CoordinateDefinitionLeaf(boxHarmonic, primitive, basis, space);
        coordinateDefinitionHarmonic.initializeCoordinates(nCells);

        move.setOmegaSquared(normalModes.getOmegaSquared());
        move.setEigenVectors(normalModes.getEigenvectors());
        move.setWaveVectors(waveVectorFactory.getWaveVectors());
        move.setWaveVectorCoefficients(waveVectorFactory.getCoefficients());
        move.setCoordinateDefinition(coordinateDefinitionHarmonic);
        move.setTemperature(temperature);

        move.setBox(boxHarmonic);


        // TARGET


        Boundary boundaryTarget;
        if (space.D() == 1) {
            boundaryTarget = new BoundaryRectangularPeriodic(space, numMolecules / density);
        } else {
            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numMolecules / 4, 1.0 / 3.0));
            boundaryTarget = new BoundaryRectangularPeriodic(space, n * L);
        }
        Box boxTarget = new Box(boundaryTarget, space);
        sim.addBox(boxTarget);
        boxTarget.setNMolecules(species, numMolecules);

        IntegratorMC integratorTarget = new IntegratorMC(potentialMasterTarget, sim.getRandom(), temperature, boxTarget);
        MCMoveAtomCoupled atomMove = new MCMoveAtomCoupled(potentialMasterTarget, new MeterPotentialEnergy(potentialMasterTarget), sim.getRandom(), space);
        atomMove.setStepSize(0.1);
        atomMove.setStepSizeMax(0.5);
        integratorTarget.getMoveManager().addMCMove(atomMove);
        ((MCMoveStepTracker) atomMove.getTracker()).setNoisyAdjustment(true);

        integrators[1] = integratorTarget;

        waveVectorFactory.makeWaveVectors(boxTarget);

        CoordinateDefinitionLeaf coordinateDefinitionTarget = new CoordinateDefinitionLeaf(boxTarget, primitive, basis, space);
        coordinateDefinitionTarget.initializeCoordinates(nCells);

        Potential2SoftSpherical potential = new P2LennardJones(space, 1.0, 1.0);
        double truncationRadius = boundaryTarget.getBoxSize().getX(0) * 0.45;
        P2SoftSphericalTruncated pTruncated = new P2SoftSphericalTruncated(space, potential, truncationRadius);
        AtomType sphereType = species.getLeafType();
        potentialMasterTarget.addPotential(pTruncated, new AtomType[]{sphereType, sphereType});
        atomMove.setPotential(pTruncated);

        potentialMasterTarget.lrcMaster().setEnabled(false);
        integratorTarget.reset();
        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(integratorTarget);
        double latticeEnergy = meterPE.getDataAsScalar();


        MeterHarmonicEnergy meterReference = new MeterHarmonicEnergy(coordinateDefinitionHarmonic, normalModes);
        MeterHarmonicEnergy meterHarmonicTarget = new MeterHarmonicEnergy(coordinateDefinitionTarget, normalModes);

        APIPotentialReference potentialReference = new APIPotentialReference(new MeterHarmonicEnergy[]{meterReference, meterHarmonicTarget});

        APIPotentialTarget potentialTarget = new APIPotentialTarget(potentialMasterTarget);
        potentialTarget.setLatticeEnergy(latticeEnergy);

        SimOverlapModule module = new SimOverlapModule(boxHarmonic, boxTarget, integrators[0], integrators[1],
                potentialReference, potentialTarget, temperature);
        module.setTargetDataInterval(numMolecules);
        module.setReferenceDataInterval(1);
        module.integratorOverlap.setDoAdjustOnTime(false);
        module.integratorOverlap.setAggressiveAdjustStepFraction(true);

        //start simulation
        module.getIntegratorOverlap().setNumSubSteps(1000);
        numSteps /= 1000;

        module.initRefPref(refFileName, numSteps / 20);
        double refPref = module.getAlphaCenter();
        if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }
        System.out.flush();

        module.equilibrate(numSteps / 10);
        refPref = module.getAlphaCenter();
        if (Double.isNaN(refPref) || refPref == 0 || Double.isInfinite(refPref)) {
            throw new RuntimeException("Simulation failed to find a valid ref pref");
        }

        System.out.println("equilibration finished");
        System.out.flush();

        module.controller.runActivityBlocking(new ActivityIntegrate(module.integratorOverlap, numSteps));

        System.out.println("final reference optimal step fraction " + module.getIntegratorOverlap().getIdealRefStepFraction() + " (actual: " + module.getIntegratorOverlap().getRefStepFraction() + ")");

        double[][] omega2 = normalModes.getOmegaSquared();
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double AHarmonic = 0;
        for (int i = 0; i < omega2.length; i++) {
            for (int j = 0; j < omega2[0].length; j++) {
                if (!Double.isInfinite(omega2[i][j])) {
                    AHarmonic += coeffs[i] * Math.log(omega2[i][j] * coeffs[i] / (temperature * Math.PI));
                }
            }
        }

        int totalCells = 1;
        for (int i = 0; i < D; i++) {
            totalCells *= nCells[i];
        }
        int basisSize = basis.getScaledCoordinates().length;
        double fac = 1;
        if (totalCells % 2 == 0) {
            fac = Math.pow(2, D);
        }
        AHarmonic -= Math.log(Math.pow(2.0, basisSize * D * (totalCells - fac) / 2.0) / Math.pow(totalCells, 0.5 * D));
        System.out.println("Harmonic-reference free energy: " + AHarmonic * temperature);

        DataOverlap dataOverlap = module.getDataOverlap();
        double ratio = dataOverlap.getOverlapAverageAndError()[0];
        double lnError = dataOverlap.getLogAverageAndError(ratio)[1];
        System.out.println("ratio average: " + ratio);
        System.out.println("free energy difference: " + (-temperature * Math.log(ratio)) + ", error: " + temperature * lnError);
        System.out.println("target free energy: " + temperature * (AHarmonic - Math.log(ratio)));
        double[] refData = dataOverlap.getAverageAndError(true, ratio);
        System.out.println("harmonic ratio average: " + refData[0] + " error: " + refData[1]);

        double[] targetData = dataOverlap.getAverageAndError(false, ratio);
        System.out.println("target ratio average: " + targetData[0] + " error: " + targetData[1]);
    }

    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 32;
        public double density = 1.0;
        public int D = 3;
        public long numSteps = 1000000;
        public String filename = "d10_T02";
        public double temperature = 0.186;
    }

    public static class MeterPotentialEnergyDifference extends DataSourceScalar {

        private static final long serialVersionUID = 1L;
        protected final DataSourceScalar meter;
        protected final double offset;
        public MeterPotentialEnergyDifference(DataSourceScalar meter, double offset) {
            super("energy", Energy.DIMENSION);
            this.meter = meter;
            this.offset = offset;
        }

        public double getDataAsScalar() {
            return meter.getDataAsScalar() - offset;
        }
    }

    public static class APIPotentialTarget implements IAPIPotential {

        protected final PotentialMaster wrappedPotentialMaster;
        private final IteratorDirective id = new IteratorDirective();
        private final PotentialCalculationEnergySum pc = new PotentialCalculationEnergySum();
        protected double latticeEnergy;

        public APIPotentialTarget(PotentialMaster wrappedPotentialMaster) {
            this.wrappedPotentialMaster = wrappedPotentialMaster;
        }

        public double calculateEnergy(Box box) {
            pc.zeroSum();
            wrappedPotentialMaster.calculate(box, id, pc);
            return pc.getSum() - latticeEnergy;
        }

        public double getLatticeEnergy() {
            return latticeEnergy;
        }

        public void setLatticeEnergy(double newLatticeEnergy) {
            latticeEnergy = newLatticeEnergy;
        }
    }

    public static class APIPotentialReference implements IAPIPotential {

        protected final MeterHarmonicEnergy[] meters;

        public APIPotentialReference(MeterHarmonicEnergy[] meters) {
            this.meters = meters;
        }

        public double calculateEnergy(Box box) {
            return meters[box.getIndex()].getDataAsScalar();
        }
    }

}
