/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMD;
import etomica.lattice.crystal.*;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.SpeciesGeneral;

/**
 * MD simulation of hard spheres in 1D or 3D with tabulation of the
 * collective-coordinate S-matrix. No graphic display of simulation.
 */
public class SimCalcS extends Simulation {

    private static final long serialVersionUID = 1L;
    public IntegratorMD integrator;

    public Box box;
    public Boundary bdry;
    public Primitive primitive;
    public CoordinateDefinition coordinateDefinition;
    public SimCalcS(Space _space, int numAtoms, double density) {
        super(_space);

        SpeciesGeneral species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        PotentialMaster potentialMaster = new PotentialMasterMonatomic(this);


        Potential potential = new P2HardSphere(space, 1.0, false);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType,
                sphereType});

        int n;
        Basis basis;
        if (space.D() == 1) {
            primitive = new PrimitiveCubic(space, 1.0 / density);
            n = numAtoms;
            bdry = new BoundaryRectangularPeriodic(space, numAtoms / density);
            ((IntegratorHard) integrator).setNullPotential(new P1HardPeriodic(space), sphereType);
            basis = new BasisMonatomic(space);
        } else {
            primitive = new PrimitiveCubic(space, 1);
            double v = primitive.unitCell().getVolume();
            primitive.scaleSize(Math.pow(v * density / 4, -1.0 / 3.0));
            n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            int[] nCells = new int[]{n, n, n};
            bdry = new BoundaryDeformableLattice(primitive, new int[]{n, n, n});
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
            primitive.scaleSize(n);

        }
        box = this.makeBox(bdry);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorHard(random, potentialMaster, box);
        integrator.setTimeStep(0.04);
        integrator.setTemperature(1.0);

        integrator.setIsothermal(false);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        // activityIntegrate.setMaxSteps(nSteps);



        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        if (space.D() == 1) {
            coordinateDefinition.initializeCoordinates(new int[]{n});
        } else {

            coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});
        }
    }

    /**
     * @param args
     */
    public static void main(String[] args) {

        // defaults
        int D = 3;
        int nA = 32;
        double density = 1.3;
        if (D == 1) {
            nA = 3;
            density = 0.5;
        }

        double simTime = 10000;

        // parse arguments
        if (args.length > 1) {
            density = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            simTime = Double.parseDouble(args[2]);
        }
        if (args.length > 3) {
            nA = Integer.parseInt(args[3]);
        }
        String filename = "normal_modes" + D + "D_"+nA+"_"+((int)(density*100))+"_cubic";
        if (args.length > 0) {
            filename = args[0];
        }

        System.out.println("Running "
                + (D == 1 ? "1D" : (D == 3 ? "FCC" : "2D hexagonal"))
                + " hard sphere simulation");
        System.out.println(nA + " atoms at density " + density);
        System.out.println(simTime + " time units");
        System.out.println("output data to " + filename);

        // construct simulation
        SimCalcS sim = new SimCalcS(Space.getInstance(D), nA, density);

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
            waveVectorFactory = new WaveVectorFactorySimple(primitive, sim.space);
        }
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);

        IntegratorListenerAction meterListener = new IntegratorListenerAction(meterNormalMode);
        meterListener.setInterval(2);
        sim.integrator.getEventManager().addListener(meterListener);

        // MeterMomentumCOM meterCOM = new MeterMomentumCOM(sim.space);
        // MeterPositionCOM meterCOM = new MeterPositionCOM(sim.space);
        // DataSinkConsole console = new DataSinkConsole();
        // DataPump comPump = new DataPump(meterCOM,console);
        // IntervalActionAdapter comAdapter = new
        // IntervalActionAdapter(comPump);
        // sim.integrator.addListener(comAdapter);
        // meterCOM.setBox(sim.box);

        // start simulation
        int nSteps = (int) (simTime / sim.integrator.getTimeStep());
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps / 10));
System.out.println("equilibration finished");
        meterNormalMode.reset();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, nSteps));

        WriteS sWriter = new WriteS(sim.space);
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setOverwrite(true);
        sWriter.actionPerformed();

    }
}
