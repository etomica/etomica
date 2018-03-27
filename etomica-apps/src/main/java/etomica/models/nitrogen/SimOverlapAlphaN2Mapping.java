/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.listener.IntegratorListenerAction;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.nbr.list.molecule.PotentialMasterListMolecular;
import etomica.normalmode.BasisBigCell;
import etomica.normalmode.MCMoveMoleculeCoupled;
import etomica.normalmode.MoleculeSiteSourceNitrogen;
import etomica.simulation.Simulation;
import etomica.space.*;
import etomica.space3d.Orientation3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;

import java.awt.*;
import java.io.File;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * Temperature-perturbation simulation for alpha-phase Nitrogen
 *
 * @author Weisong Lin
 */
public class SimOverlapAlphaN2Mapping extends Simulation {

    protected final MoleculeAgentManager latticeCoordinates;
    public IntegratorMC integrator;
    public ActivityIntegrate activityIntegrate;
    public Box box;
    public Primitive primitive;
    protected PotentialMasterListMolecular potentialMaster;
    protected double latticeEnergy;
    protected SpeciesN2 species;
    protected CoordinateDefinitionNitrogen coordinateDef;
    protected P2Nitrogen potential;
    protected PRotConstraint pRotConstraint;
    protected MCMoveRotateMolecule3DN2AveCosThetaConstraint rotate;

    public SimOverlapAlphaN2Mapping(Space space, int[] nC, double density, double temperature,
                                    long numSteps, double rcScale, double constraintAngle, boolean noRotScale, boolean doRotation, boolean doTranslation) {
        super(space);

        int numMolecules = nC[0] * nC[1] * nC[2] * 4;

        double a = Math.pow(4.0 / density, 1.0 / 3.0);
//        System.out.println("Unit Cell Length, a: " + a);

        Basis basisFCC = new BasisCubicFcc();
        Basis basis = new BasisBigCell(space, basisFCC, new int[]{nC[0], nC[1], nC[2]});

        species = new SpeciesN2(space);
        addSpecies(species);

        box = new Box(space);
        addBox(box);
        box.setNMolecules(species, numMolecules);

        int[] nCells = new int[]{1, 1, 1};
        double[] boxSize = new double[]{nC[0] * a, nC[1] * a, nC[2] * a};
        Boundary boundary = new BoundaryRectangularPeriodic(space, boxSize);
        primitive = new PrimitiveTetragonal(space, nC[0] * a, nC[2] * a);

        coordinateDef = new CoordinateDefinitionNitrogen(this, box, primitive, basis, space);
        coordinateDef.setIsAlpha();
        coordinateDef.setOrientationVectorAlpha(space);
        coordinateDef.initializeCoordinates(nCells);
        box.setBoundary(boundary);

        latticeCoordinates = new MoleculeAgentManager(this, box, new MoleculeSiteSourceNitrogen(space, new MoleculePositionCOM(space), new NitrogenOrientationDefinition(space)));
        double rC = box.getBoundary().getBoxSize().getX(0) * rcScale;
//        System.out.println("Truncation Radius (" + rcScale + " Box Length): " + rC);
        potential = new P2Nitrogen(space, rC);
        potential.setBox(box);

        pRotConstraint = new PRotConstraint(space, coordinateDef, box);
        pRotConstraint.setConstraintAngle(constraintAngle);

        //potentialMaster = new PotentialMaster();
        potentialMaster = new PotentialMasterListMolecular(this, space);
        potentialMaster.addPotential(potential, new ISpecies[]{species, species});
        if (!noRotScale) {
//            System.out.println("set constraint angle to = " + constraintAngle);
            potentialMaster.addPotential(pRotConstraint, new ISpecies[]{species});
        }

        int cellRange = 6;
        potentialMaster.setRange(rC);
        potentialMaster.setCellRange(cellRange);
        potentialMaster.getNeighborManager(box).reset();

        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

//        int numNeigh = potentialMaster.getNeighborManager(box).getDownList(box.getMoleculeList().getMolecule(0))[0].getMoleculeCount();
//        System.out.println("numNeigh: " + numNeigh);

        MCMoveMoleculeCoupled move = new MCMoveMoleculeCoupled(potentialMaster, getRandom(), space);
        move.setBox(box);
        move.setPotential(potential);
        move.setDoExcludeNonNeighbors(true);

        rotate = new MCMoveRotateMolecule3DN2AveCosThetaConstraint(potentialMaster, getRandom(), space, coordinateDef, 0.8);
        rotate.setBox(box);

        integrator = new IntegratorMC(potentialMaster, getRandom(), Kelvin.UNIT.toSim(temperature));
        if (doTranslation) integrator.getMoveManager().addMCMove(move);
        if (doRotation) integrator.getMoveManager().addMCMove(rotate);
        integrator.setBox(box);

        potential.setRange(Double.POSITIVE_INFINITY);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster);
        meterPE.setBox(box);
        latticeEnergy = meterPE.getDataAsScalar();
//        System.out.println("lattice energy per molecule (sim unit): " + latticeEnergy / numMolecules);

        potential.setRange(rC);
        potential.setRange(Double.POSITIVE_INFINITY);


        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

    }

    /**
     * @param args filename containing simulation parameters
     * @see SimOverlapAlphaN2Mapping.SimOverlapParam
     */
    public static void main(String[] args) {
        final long startTime = System.currentTimeMillis();
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        ParseArgs.doParseArgs(params, args);
        double density = params.density;
        int numSteps = params.numSteps;
        final int numMolecules = params.numMolecules;
        double temperature = params.temperature;
        int[] nC = params.nC;
        double rcScale = params.rcScale;
        double constraintAngle = params.constraintAngle;
        boolean noRotScale = params.noRotScale;
        boolean doTranslation = params.doTranslation;
        boolean doRotation = params.doRotation;
        boolean runGraphic = params.runGraphic;
        boolean doMapping = params.doMapping;
        boolean doConventional = params.doConventional;
        String configFileName = "configT" + temperature;
        String filename = "alphaN2d" + density + "_T" + temperature + "Cons0.8";


        //instantiate simulation
        final SimOverlapAlphaN2Mapping sim = new SimOverlapAlphaN2Mapping(Space.getInstance(3), nC, density, temperature,
                numSteps, rcScale, constraintAngle, noRotScale, doRotation, doTranslation);

        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE2.setBox(sim.box);
        final double latticeEnergy = meterPE2.getDataAsScalar();
        System.out.println("latticeEnergy = " + latticeEnergy);

        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator(sim.integrator);

        MeterDADBNitrogen meterDADB = null;
        if (doMapping) {
            meterDADB = new MeterDADBNitrogen(sim, meterPE, sim.potentialMaster, Kelvin.UNIT.toSim(temperature), sim.latticeCoordinates);
        }
        long initStep = (1 + (numMolecules / 500)) * 100 * numMolecules;
        sim.initialize(initStep);
        System.out.flush();


        if (runGraphic) {
            sim.activityIntegrate.setMaxSteps(Long.MAX_VALUE);
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, sim.space, sim.getController());
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors == null) {
                        allColors = new Color[768];
                        for (int i = 0; i < 256; i++) {
                            allColors[i] = new Color(255 - i, i, 0);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 256] = new Color(0, 255 - i, i);
                        }
                        for (int i = 0; i < 256; i++) {
                            allColors[i + 512] = new Color(i, 0, 255 - i);
                        }
                    }
                    return allColors[(2 * a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());

            simGraphic.makeAndDisplayFrame("TP Alpha Nitrogen");
            return;
        }


        if (doMapping) {
            meterDADB.doTranslation = doTranslation;
            meterDADB.doRotation = doRotation;
        }


        AccumulatorAverageFixed accumulatorAverageFixedDADB = null;
        DataPumpListener dataPumpListenerDADB = null;
        if (doMapping) {
            int DADBBlockSize = numSteps / (4 * numMolecules * 100);
            if (DADBBlockSize == 0) DADBBlockSize = 1;
            accumulatorAverageFixedDADB = new AccumulatorAverageFixed(DADBBlockSize);
            dataPumpListenerDADB = new DataPumpListener(meterDADB, accumulatorAverageFixedDADB, numMolecules * 4);
        }

        AccumulatorAverageFixed accumulatorAverageFixedPE = null;
        DataPumpListener dataPumpListenerPE = null;
        if (doConventional) {
            int PEBlockSize = numSteps / (10 * 100);
            if (PEBlockSize == 0) PEBlockSize = 1;
            accumulatorAverageFixedPE = new AccumulatorAverageFixed(PEBlockSize);
            dataPumpListenerPE = new DataPumpListener(meterPE, accumulatorAverageFixedPE, 1);
        }
        sim.activityIntegrate.setMaxSteps(numSteps / 10);
        sim.getController().actionPerformed();

        sim.activityIntegrate.setMaxSteps(numSteps);

        if (doMapping) sim.integrator.getEventManager().addListener(dataPumpListenerDADB);
        if (doConventional) sim.integrator.getEventManager().addListener(dataPumpListenerPE);
        sim.getController().reset();
        sim.getController().actionPerformed();


        long endTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println(date.format(cal.getTime()));
        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
        System.out.println("Time taken (in mins): " + totalTime);
        System.out.println("numSteps = " + numSteps);
        System.out.println("temperature = " + temperature);

        if (doMapping) {
            double MappingAverage = accumulatorAverageFixedDADB.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double MappingErr = accumulatorAverageFixedDADB.getData(AccumulatorAverage.ERROR).getValue(0);
            double MappingCor = accumulatorAverageFixedDADB.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
            System.out.println("MappingAverage= " + MappingAverage + "\tMappingErr= "
                    + MappingErr + "\tMappingCor= " + MappingCor + "\ttotalTime= " + totalTime);
        }

        int N = sim.box.getMoleculeList().getMoleculeCount();
        double fac = (doRotation ? 1.0 : 0) * N + (doTranslation ? 1.5 : 0) * (N - 1);
        double harmonicEnergy = fac * Kelvin.UNIT.toSim((temperature));
//        System.out.println("doRotation:" + doRotation + " doTranslation:" + doTranslation);


        if (doConventional) {
            double PEAverage = accumulatorAverageFixedPE.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double PEErr = accumulatorAverageFixedPE.getData(AccumulatorAverage.ERROR).getValue(0);
            double PECor = accumulatorAverageFixedPE.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);
            System.out.println("PEAverage= " + (PEAverage - latticeEnergy - harmonicEnergy)
                    + "\tPEErr= " + PEErr + "\tPECor= " + PECor + "\ttotalTime= " + totalTime);
        }
//        System.out.println("PE full = " + PEAverage);
//        System.out.println("Harmonic = " + harmonicEnergy);
        final double endLatticeEnergy = meterPE2.getDataAsScalar();
        System.out.println("endLE = " + endLatticeEnergy);
    }

    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomolous contributions
        System.out.println("\nEquilibration Steps: " + initSteps + " + numSteps/10");
        activityIntegrate.setMaxSteps(initSteps);
        getController().actionPerformed();
        getController().reset();
    }



    public static class NitrogenOrientationDefinition implements etomica.normalmode.MoleculeSiteSourceNitrogen.MoleculeOrientationDefinition {
        protected final Orientation3D or;
        protected final Vector v1;

        public NitrogenOrientationDefinition(Space space) {
            or = new Orientation3D(space);
            v1 = space.makeVector();

        }

        public IOrientation getOrientation(IMolecule molecule) {
            IAtomList leafList = molecule.getChildList();
            Vector n1 = leafList.getAtom(0).getPosition();
            Vector n2 = leafList.getAtom(1).getPosition();
            v1.Ev1Mv2(n2, n1);
            v1.normalize();
            or.setDirection(v1);
            return or;
        }
    }

    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 256;
        public int[] nC = new int[]{4, 4, 4};
        public double density = 0.023; //0.02204857502170207 (intial from literature with a = 5.661)
        public int numSteps = 100000;
        public double temperature = 20; // in unit Kelvin
        public double rcScale = 0.475;
        public double constraintAngle = 70;
        public boolean noRotScale = false;
        public boolean runGraphic = false;
        public boolean doRotation = true;
        public boolean doTranslation = true;
        public boolean doMapping = true;
        public boolean doConventional = false;
    }
}
