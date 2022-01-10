/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorPT;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.util.ParameterBase;

import java.awt.*;


public class SimHSMSD extends Simulation {
    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public IntegratorHard integrator;
    public PotentialMaster potentialMaster;
    public IntegratorPT.MCMoveSwap swapMove;
    public P2HardSphere potential;
    public SpeciesGeneral species;
    public CoordinateDefinitionLeaf coordinateDefinition;


    public SimHSMSD(Space _space, int numAtoms, double density, double tStep) {
        super(Space3D.getInstance());
        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this), true);
        addSpecies(species);

        double L = Math.pow(4.0 / density, 1.0 / 3.0);
        int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
        Boundary boundary = new BoundaryRectangularPeriodic(space, n * L);
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        PrimitiveCubic primitive = new PrimitiveCubic(space, n * L);
        int[] nCells = new int[]{n, n, n};
        Basis basisFCC = new BasisCubicFcc();
        BasisBigCell basis = new BasisBigCell(space, basisFCC, nCells);

        double neighborRangeFac = 1.6;
        double sigma = 1.0;

        potentialMaster = new PotentialMasterList(this, sigma * neighborRangeFac, space);
        integrator = new IntegratorHard(potentialMaster, random, tStep, 1.0, box);
        integrator.setIsothermal(false);

        getController().addActivity(new ActivityIntegrate(integrator));

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();

        potential = new P2HardSphere(space, 1.0, false);
        AtomType leafType = species.getLeafType();

        potentialMaster.addPotential(potential, new AtomType[]{leafType, leafType});

        coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

        NeighborListManager nbrManager = ((PotentialMasterList) potentialMaster).getNeighborManager(box);
        integrator.getEventManager().addListener(nbrManager);
    }

    /**
     * @param args filename containing simulation parameters
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();

        double density = params.density;
        long numSteps = params.numSteps;
        final int numAtoms = params.numAtoms;
        double tStep = params.tStep;

        System.out.println(numAtoms+" atoms at density "+density);
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimHSMSD sim = new SimHSMSD(Space.getInstance(3), numAtoms, density, tStep);
        if (false) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            simGraphic.setPaintInterval(sim.box, 1000);
            ColorScheme colorScheme = new ColorScheme() {
                protected Color[] allColors;

                public Color getAtomColor(IAtom a) {
                    if (allColors==null) {
                        allColors = new Color[768];
                        for (int i=0; i<256; i++) {
                            allColors[i] = new Color(255-i,i,0);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+256] = new Color(0,255-i,i);
                        }
                        for (int i=0; i<256; i++) {
                            allColors[i+512] = new Color(i,0,255-i);
                        }
                    }
                    return allColors[(2*a.getLeafIndex()) % 768];
                }
            };
            simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            DisplayTextBox timer = new DisplayTextBox();
            DataSourceCountSteps counter = new DataSourceCountSteps(sim.integrator);
            DataPumpListener counterPump = new DataPumpListener(counter, timer, 100);
            sim.integrator.getEventManager().addListener(counterPump);
            simGraphic.getPanel().controlPanel.add(timer.graphic());
            simGraphic.makeAndDisplayFrame("HS");
            return;
        }

        //start simulation
        double L = Math.pow(numAtoms, 1.0/3.0);

        //Mapped
        MeterSolidPropsLJ meterSolid = new  MeterSolidPropsLJ(sim.space, new MeterPotentialEnergyFromIntegrator(sim.integrator),  sim.potentialMaster, sim.coordinateDefinition, 1.0,new double[11]);

//Initialization
        System.out.flush();
        if (args.length == 0) {
            // quick initialization
            sim.initialize(numSteps/10);
        }
        else {
            long nSteps = numSteps/20 + 50*numAtoms + numAtoms*numAtoms*3;
            if (nSteps > numSteps/2) nSteps = numSteps/2;
            sim.initialize(nSteps);
        }

        int numBlocks = 100;
        int interval = 2*numAtoms;
        int intervalLS = 5*interval;
        long blockSize = numSteps/(numBlocks*interval);
        if (blockSize == 0) blockSize = 1;
        long blockSizeLS = numSteps/(numBlocks*intervalLS);
        if (blockSizeLS == 0) blockSizeLS = 1;
        int o=2;
        while (blockSize<numSteps/5 && (numSteps != numBlocks*intervalLS*blockSizeLS || numSteps != numBlocks*interval*blockSize)) {
            interval = 2*numAtoms+(o%2==0 ? (o/2) : -(o/2));
            if (interval < 1 || interval > numSteps/5) {
                throw new RuntimeException("oops interval "+interval);
            }
            // only need to enforce intervalLS if nCutoffsLS>0.  whatever.
            intervalLS = 5*interval;
            blockSize = numSteps/(numBlocks*interval);
            if (blockSize == 0) blockSize = 1;
            blockSizeLS = numSteps/(numBlocks*intervalLS);
            if (blockSizeLS == 0) blockSizeLS = 1;
            o++;
        }
        if (numSteps != numBlocks*intervalLS*blockSizeLS || numSteps != numBlocks*interval*blockSize) {
            throw new RuntimeException("unable to find appropriate intervals");
        }
        System.out.println("block size "+blockSize+" interval "+interval);


        //Mapped
        AccumulatorAverageFixed accumulatorUP = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorUPPump = new DataPumpListener(meterSolid, accumulatorUP, interval);
        sim.integrator.getEventManager().addListener(accumulatorUPPump);

        final long startTime = System.currentTimeMillis();

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, numSteps));

        System.out.println();

        IData dataUPAvg = accumulatorUP.getData(accumulatorUP.AVERAGE);
        IData dataUPErr = accumulatorUP.getData(accumulatorUP.ERROR);
        IData dataUPCorr = accumulatorUP.getData(accumulatorUP.BLOCK_CORRELATION);

        double u2   = dataUPAvg.getValue(0);
        double u2Err = dataUPErr.getValue(0);
        double u2Corr = dataUPCorr.getValue(0);

        double u2_alpha   = dataUPAvg.getValue(1);
        double u2_alphaErr = dataUPErr.getValue(1);
        double u2_alphaCorr = dataUPCorr.getValue(1);

        System.out.println();
        System.out.println("************************************************************************");
        System.out.println("1 1 1");
        System.out.println("u2        : " + u2 + " err: " + u2Err + " corr: " + u2Corr);
        System.out.println("u2_alpha  : " + u2_alpha + " err: " + u2_alphaErr + " corr: " + u2_alphaCorr);


        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);
    }



    public void initialize(long initSteps) {
        // equilibrate off the lattice to avoid anomalous contributions
        this.getController().runActivityBlocking(new ActivityIntegrate(this.integrator, initSteps));
    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numAtoms = 500;
        public double density = 1.2;
        public long numSteps = 1000000;
        public double tStep = 0.005;
    }
}
