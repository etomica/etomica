/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import java.awt.Color;

import etomica.action.activity.ActivityIntegrate;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.SimulationGraphic;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.potential.PotentialMolecular;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSumPolarizable;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.cluster.Standard;


public class VirialGCPMGraphic {

    public static void main(String[] args) {

        int nPoints = 5;
        double temperature = Kelvin.UNIT.toSim(300);
        int numSubSteps = 1000;
        double deltaCut = 100;

        double sigmaHSRef = 3.2;
        final double[] HSB = new double[7];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);

        System.out.println("sigmaHSRef: " + sigmaHSRef);
        System.out.println("B2HS: " + HSB[2]);
        System.out.println("B3HS: " + HSB[3] + " = " + (HSB[3] / (HSB[2] * HSB[2])) + " B2HS^2");
        System.out.println("B4HS: " + HSB[4] + " = " + (HSB[4] / (HSB[2] * HSB[2] * HSB[2])) + " B2HS^3");
        System.out.println("B5HS: " + HSB[5] + " = 0.110252 B2HS^4");
        System.out.println("B6HS: " + HSB[6] + " = 0.03881 B2HS^5");
        System.out.println("Water Direct sampling B" + nPoints + " at T=" + Kelvin.UNIT.fromSim(temperature));

        Space space = Space3D.getInstance();

        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        final PotentialMolecular pTarget = new PNWaterGCPM(space);

        MayerGeneral fTarget = new MayerGeneral(pTarget);
        ClusterAbstract targetCluster = Standard.virialClusterPolarizable(nPoints, fTarget, nPoints > 3, false);
        ((ClusterSumPolarizable) targetCluster).setDeltaCut(deltaCut);

        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints > 3, null, false);


        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        SpeciesWater4P species = new SpeciesWater4P(space);
        species.setConformation(new ConformationWaterGCPM(space));

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, temperature, refCluster, targetCluster, false);
        sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
        sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
        simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);

        ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
        ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(0), Color.WHITE);
        ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[0]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
        ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box[1]).getColorScheme()).setColor(species.getAtomType(1), Color.RED);
        simGraphic.makeAndDisplayFrame();

        sim.integratorOS.setNumSubSteps(numSubSteps);
        sim.setAccumulatorBlockSize(1000);

        // if running interactively, set filename to null so that it doens't read
        // (or write) to a refpref file
        sim.initRefPref(null, 10);
        sim.equilibrate(null, 20);
        sim.getController().addActivity(new ActivityIntegrate(sim.integratorOS));
        if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
            throw new RuntimeException("Oops");
        }

//        sim.ai.setMaxSteps(Long.MAX_VALUE);
//        sim.getController().actionPerformed();
//
//        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
//        double ratio = sim.dsvo.getDataAsScalar();
//        double error = sim.dsvo.getError();
//        System.out.println("ratio average: "+ratio+", error: "+error);
//        System.out.println("abs average: "+ratio*HSB[nPoints]+", error: "+error*HSB[nPoints]);
//
//        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
//        System.out.println("hard sphere ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[1]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index)).getData()[1]);
//        System.out.println("hard sphere   average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
//                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
//        System.out.println("hard sphere overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
//                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);
//
//        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.dsvo.minDiffLocation());
//        System.out.println("water ratio average: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO.index)).getData()[1]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorRatioAverageCovariance.StatType.RATIO_ERROR.index)).getData()[1]);
//        System.out.println("water average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[0]
//                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[0]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[0]);
//        System.out.println("water overlap average: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.AVERAGE.index)).getData()[1]
//                          +" stdev: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.STANDARD_DEVIATION.index)).getData()[1]
//                          +" error: "+((DataDoubleArray)allYourBase.getData(AccumulatorAverage.StatType.ERROR.index)).getData()[1]);

    }

}



