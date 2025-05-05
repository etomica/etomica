/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.data.AccumulatorHistogram;
import etomica.data.histogram.HistogramExpanding;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.models.water.P2WaterSPC;
import etomica.models.water.SpeciesWater3P;
import etomica.potential.P2Electrostatic;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalSum;
import etomica.potential.PotentialMoleculePair;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.units.Kelvin;
import etomica.units.dimensions.Length;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.ClusterBonds;
import etomica.virial.cluster.ClusterCoupledFlippedPoints;
import etomica.virial.cluster.ClusterSum;
import etomica.virial.cluster.Standard;
import etomica.virial.mcmove.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.wheatley.ClusterWheatleyHS;

import java.awt.*;

/**
 * Computes CO2-H2O mixture virial coefficients using ab-initio potentials
 * for both components.  Classical and semiclassical coefficients can be
 * computed.
 * 
 * For water, only pairwise-additive contributions are computed.
 * 
 * @author Andrew Schultz
 */
public class VirialFlipDemo {


    public static void main(String[] args) {

        VirialCO2SCParam params = new VirialCO2SCParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs.doParseArgs(params, args);
        }
        else {
            // customize parameters here
            params.temperature = 800;
        }

        final double temperatureK = params.temperature;
        double sigmaHSRef = 3;
        int nPoints = 3;

        final double HSB = Standard.BHS(nPoints, sigmaHSRef);

        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);

        ISpecies species = SpeciesWater3P.create();

        SpeciesManager.Builder sb = SpeciesManager.builder();

        sb.addSpecies(species);
        SpeciesManager sm = sb.build();

        AtomType typeH = species.getAtomType(0);
        AtomType typeO = species.getAtomType(1);
        PotentialMoleculePair p2 = new PotentialMoleculePair(space, sm);
        p2.setAtomPotential(typeH, typeH, new P2Electrostatic(P2WaterSPC.chargeH, P2WaterSPC.chargeH));
        P2SoftSphericalSum p2OO = new P2SoftSphericalSum(new P2Electrostatic(P2WaterSPC.chargeO, P2WaterSPC.chargeO),
                                                         new P2LennardJones(P2WaterSPC.sigmaOO, P2WaterSPC.epsilonOO));
        p2.setAtomPotential(typeO, typeO, p2OO);
        p2.setAtomPotential(typeO, typeH, new P2Electrostatic(P2WaterSPC.chargeO, P2WaterSPC.chargeH));

        MayerGeneral fTarget = new MayerGeneral(p2);
        ClusterBonds cb = new ClusterBonds(5, new int[][][]{{{0, 1}, {0, 2}}});
        ClusterSum targetClusterSingle = new ClusterSum(new ClusterBonds[]{cb}, new double[]{1}, new MayerFunction[]{fTarget});
        targetClusterSingle.setCaching(false);
        int[][] flipPoints = new int[][]{{0,1}};
        ClusterCoupledFlippedPoints targetCluster = new ClusterCoupledFlippedPoints(targetClusterSingle, space, flipPoints, 0);

        targetCluster.setTemperature(temperature);

        ClusterWheatleyHS refCluster = new ClusterWheatleyHS(nPoints, fRef);

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, 3, temperature, refCluster, targetCluster);
        sim.init();

        ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).startMolecule = 1;

        sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{40, 40, 40}));
        sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{40, 40, 40}));
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]);
        DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
        displayBox0.setShowBoundary(false);
        displayBox1.setShowBoundary(false);
        ((DisplayBoxCanvasG3DSys) displayBox0.canvas).setBackgroundColor(Color.WHITE);
        ((DisplayBoxCanvasG3DSys) displayBox1.canvas).setBackgroundColor(Color.WHITE);

        ColorScheme colorScheme = new ColorScheme() {
            @Override
            public Color getAtomColor(IAtom a) {
                if (a.getIndex() < 2) return Color.WHITE;
                Color[] colorO = new Color[]{Color.RED, Color.BLUE, Color.GREEN};
                return colorO[a.getParentGroup().getIndex()];
            }
        };
        displayBox0.setColorScheme(colorScheme);
        displayBox1.setColorScheme(colorScheme);

        sim.integratorOS.setNumSubSteps(1000);
        sim.setAccumulatorBlockSize(1000);

        // if running interactively, set filename to null so that it doens't read
        // (or write) to a refpref file
        sim.initRefPref(null, 10, false);
        sim.equilibrate(null, 20, false);
        sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }

        AccumulatorHistogram accHist01 = new AccumulatorHistogram(new HistogramExpanding(1), 100);
        accHist01.putDataInfo(new DataDouble.DataInfoDouble("histogram 01", Length.DIMENSION));
        AccumulatorHistogram accHist02 = new AccumulatorHistogram(new HistogramExpanding(1), 100);
        accHist02.putDataInfo(new DataDouble.DataInfoDouble("histogram 02", Length.DIMENSION));
        BoxCluster box1 = sim.box[1];
        IntegratorListener histListenerTarget = new IntegratorListener() {

            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                CoordinatePairSet cPairs = box1.getCPairSet();

                double r01 = Math.sqrt(cPairs.getr2(0, 1));
                double r02 = Math.sqrt(cPairs.getr2(0, 2));
                DataDouble d = new DataDouble();
                d.x = r01;
                accHist01.putData(d);
                d.x = r02;
                accHist02.putData(d);
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        sim.integrators[1].getEventManager().addListener(histListenerTarget);
        DisplayPlotXChart plotHist = new DisplayPlotXChart();
        accHist01.addDataSink(plotHist.makeSink("histogram 01"));
        accHist02.addDataSink(plotHist.makeSink("histogram 02"));
        plotHist.setLabel("histograms");
        plotHist.setYLog(true);
        simGraphic.add(plotHist);


        simGraphic.makeAndDisplayFrame();
    }
    /**
     * Inner class for parameters
     */
    public static class VirialCO2SCParam extends ParameterBase {
        // don't change these
        public double temperature = 100;
    }
}
