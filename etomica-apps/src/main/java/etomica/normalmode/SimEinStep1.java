/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;


import etomica.action.activity.ActivityIntegrate2;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingDiscard;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorScheme;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.lattice.crystal.*;
import etomica.nbr.list.NeighborListManagerSlanty;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformablePeriodic;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Simulation that samples a composite energy function (soft sphere and
 * Einstein crystal) and perturbs into overlap regions shared by systems with
 * more or less soft sphere contributions.
 * 
 * @author Andrew Schultz
 */
public class SimEinStep1 extends Simulation {

    public final PotentialMasterList potentialMaster;
    public IntegratorMC integrator;

    public Box box;
    public Boundary boundary;
    public int[] nCells;
    public Basis basis;
    public Primitive primitive;
    public MCMoveEinsteinCrystal atomMove;
    public PotentialMasterMonatomic potentialMasterHarmonic;
    public SimEinStep1(Space _space, final int numAtoms, double density, final double temperature, double spring, int exponent, double rc, boolean slanty) {
        super(_space);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        addSpecies(species);

        if (slanty) {
            potentialMaster = new PotentialMasterList(this, rc, new NeighborListManagerSlanty.NeighborListSlantyAgentSource(rc), space);
        } else {
            potentialMaster = new PotentialMasterList(this, space);
        }

        // TARGET
        if (slanty) {
            int c = (int) Math.round(Math.pow(numAtoms, 1.0 / 3.0));
            nCells = new int[]{c, c, c};

            double L = Math.pow(Math.sqrt(2) / density, 1.0 / 3.0);
            double angle = Math.PI / 3;

//            primitive = new PrimitiveFcc(space, L*c);
            primitive = new PrimitiveTriclinic(space, L * c, L * c, L * c, angle, angle, angle);

            boundary = new BoundaryDeformablePeriodic(space, primitive.vectors());
            ((BoundaryDeformablePeriodic) boundary).setTruncationRadius(rc);
            Basis basisSimple = new Basis(new Vector3D[]{new Vector3D(0.0, 0.0, 0.0)});
            basis = new BasisBigCell(space, basisSimple, nCells);
        } else {

            double L = Math.pow(4.0 / density, 1.0 / 3.0);
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            primitive = new PrimitiveCubic(space, n * L);

            nCells = new int[]{n, n, n};
            boundary = new BoundaryRectangularPeriodic(space, n * L);
            Basis basisFCC = new BasisCubicFcc();
            basis = new BasisBigCell(space, basisFCC, nCells);
        }
        box = this.makeBox(boundary);
        box.setNMolecules(species, numAtoms);

        integrator = new IntegratorMC(potentialMaster, getRandom(), temperature, box);

        CoordinateDefinitionLeaf coordinateDefinition = new CoordinateDefinitionLeaf(box, primitive, basis, space);
        coordinateDefinition.initializeCoordinates(new int[]{1, 1, 1});

//        ConfigurationFile configEC = new ConfigurationFile("ec");
//        configEC.initializeCoordinates(box);

        Potential2SoftSpherical potential = exponent > 0 ? new P2SoftSphere(space, 1.0, 1.0, exponent) : new P2LennardJones(space);
        potential = new P2SoftSphericalTruncated(space, potential, rc);
        AtomType sphereType = species.getLeafType();
        potentialMaster.addPotential(potential, new AtomType[]{sphereType, sphereType});


        potentialMaster.lrcMaster().setEnabled(false);

        int cellRange = 7;
        potentialMaster.setRange(rc);
        potentialMaster.setCellRange(cellRange); // insanely high, this lets us have neighborRange close to dimensions/2
        // find neighbors now.  Don't hook up NeighborListManager (neighbors won't change)
        potentialMaster.getNeighborManager(box).reset();
        int potentialCells = potentialMaster.getNbrCellManager(box).getLattice().getSize()[0];
        if (potentialCells < cellRange * 2 + 1) {
            throw new RuntimeException("oops (" + potentialCells + " < " + (cellRange * 2 + 1) + ")");
        }

        if (false) {
            P1HarmonicSite p1Harmonic = new P1HarmonicSite(space);
            p1Harmonic.setSpringConstant(spring);
            p1Harmonic.setAtomAgentManager(box, coordinateDefinition.siteManager);
            potentialMasterHarmonic = new PotentialMasterMonatomic(this);
            potentialMasterHarmonic.addPotential(p1Harmonic, new AtomType[]{sphereType, sphereType});
        }

        atomMove = new MCMoveEinsteinCrystal(space, random);
        atomMove.setCoordinateDefinition(coordinateDefinition);
        atomMove.setAlphaEin(spring);
        atomMove.setTemperature(temperature);
        integrator.getMoveManager().addMCMove(atomMove);

        this.getController2().addActivity(new ActivityIntegrate2(integrator));

        // extend potential range, so that atoms that move outside the truncation range will still interact
        // atoms that move in will not interact since they won't be neighbors
        //XXX we don't want to do this because our potential is shifted!
        ((P2SoftSphericalTruncated) potential).setTruncationRadius(0.6 * boundary.getBoxSize().getX(0));
    }

    /**
     * @param args filename containing simulation parameters
     * @see SimEinStep1.SimOverlapParam
     */
    public static void main(String[] args) {
        //set up simulation parameters
        SimOverlapParam params = new SimOverlapParam();
        if (args.length == 0) {
            params.numMolecules = 500;
            params.slanty = false;
            params.rc = 3.0;
            params.exponentN = 0;
            params.density = 1;
            params.temperature = 1.0;
            params.numSteps = 1000000;
            params.spring = 50000;
//            params.spring = 5.694733676007302e+03;
        }
        else {
            ParseArgs.doParseArgs(params, args);
        }
        double density = params.density;
        int exponentN = params.exponentN;
        final int numMolecules = params.numMolecules;
        long numSteps = params.numSteps/numMolecules;
        final double temperature = params.temperature;
        double rc = params.rc;
        double spring = params.spring;
        boolean slanty = params.slanty;
        double x0 = params.x0;
        double f = params.f;

        double c = Math.exp(x0);
        double xf = Math.log(spring + c);
        double x=x0+(xf-x0)*f;
        spring=(Math.exp(x)-c);


        System.out.println("Running Einstein crystal simulation (step 1) with spring="+spring);
        System.out.println(numMolecules+" atoms at density "+density+" and temperature "+temperature);
        if (exponentN > 0) {
            System.out.println("Soft spheres, exponent N: "+ exponentN);
        }
        else {
            System.out.println("Lennard-Jones");
        }
        System.out.println(numSteps+" steps");

        //instantiate simulation
        final SimEinStep1 sim = new SimEinStep1(Space3D.getInstance(), numMolecules, density, temperature, spring, exponentN, rc, slanty);

        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        final double latticeEnergy = meterPE.getDataAsScalar();
        System.out.println("uLat "+latticeEnergy/numMolecules);
        System.out.println("buLat "+latticeEnergy/numMolecules/temperature);

        DataSourceScalar meter = new DataSourceScalar("foo", Null.DIMENSION) {

            public double getDataAsScalar() {
                double pe = meterPE.getDataAsScalar() - latticeEnergy;
//                System.out.println(pe/numAtoms+" "+pe/temperature+" "+Math.exp(-pe/temperature));
                return Math.exp(-pe/temperature);
            }

        };

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

            AccumulatorHistory peHist = new AccumulatorHistory(new HistoryCollapsingDiscard());
            DataPumpListener accumulatorPump = new DataPumpListener(meter, peHist);
            sim.integrator.getEventManager().addListener(accumulatorPump);
            DisplayPlot pePlot = new DisplayPlot();
            peHist.setDataSink(pePlot.getDataSet().makeDataSink());
            pePlot.setLabel("PE");
            simGraphic.add(pePlot);

            simGraphic.makeAndDisplayFrame();
            return;
        }

        //start simulation

        System.out.flush();

        int numBlocks = 100;
        long blockSize = numSteps/numBlocks;
        if (blockSize == 0) blockSize = 1;
        System.out.println("block size "+blockSize);
        AccumulatorAverageFixed accumulator = new AccumulatorAverageFixed(blockSize);
        DataPumpListener accumulatorPump = new DataPumpListener(meter, accumulator);
        sim.integrator.getEventManager().addListener(accumulatorPump);

        final long startTime = System.currentTimeMillis();

        sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), numSteps);

        //MeterTargetTP.openFW("x"+numMolecules+".dat");
        //MeterTargetTP.closeFW();

        DataGroup data = (DataGroup)accumulator.getData();
        IData dataErr = data.getData(accumulator.ERROR.index);
        IData dataAvg = data.getData(accumulator.AVERAGE.index);
        IData dataCorrelation = data.getData(accumulator.BLOCK_CORRELATION.index);
        double avg = dataAvg.getValue(0);
        double err = dataErr.getValue(0);
        double cor = dataCorrelation.getValue(0);
        System.out.println("Qratio  "+avg+" "+err+" "+cor);
        System.out.println("dbetaA  "+(-Math.log(avg))+" "+err/avg);
        System.out.println("dbetaA/N  "+(-Math.log(avg)/numMolecules)+" "+err/avg/numMolecules);
        System.out.println("betaA/N "+(latticeEnergy/numMolecules/temperature-Math.log(avg)/numMolecules));

        long endTime = System.currentTimeMillis();
        System.out.println("time: " + (endTime - startTime)/1000.0);
    }
    
    /**
     * Inner class for parameters understood by the HSMD3D constructor
     */
    public static class SimOverlapParam extends ParameterBase {
        public int numMolecules = 256;
        public double density = 1.28;
        public int exponentN = 0;
        public long numSteps = 1000000;
        public double temperature = 2;
        public double rc = 2.7;
        public double spring = 2500*240*(3.405*3.405)/120;
        public boolean slanty = false;
        public double f = 1.0;
        public double x0 = 3.5;
    }
}
