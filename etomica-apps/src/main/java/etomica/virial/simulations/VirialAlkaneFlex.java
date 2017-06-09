/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.graphics.*;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2LennardJones;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.*;
import etomica.units.Dimension;
import etomica.util.Arrays;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.util.ReadParameters;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import javax.swing.*;
import java.awt.*;
import java.io.File;

/**
 * Mayer sampling simulation for alkanes using the TraPPE force field.
 *   M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase
 *   Equilibria. 1. United-Atom Description of n-Alkanes," J. Phys. Chem. B
 *   102, 2569-2577 (1998)
 */
public class VirialAlkaneFlex {


    public static void main(String[] args) {
        VirialSiepmannSpheresParam params = new VirialSiepmannSpheresParam();
        boolean isCommandline = false;
        if (args.length > 0) {
            isCommandline = true;
            if (new File(args[0]).exists()) {
                ReadParameters paramReader = new ReadParameters(args[0], params);
                paramReader.readParameters();
                args = (String[])Arrays.removeObject(args, args[0]);
            }
            if (args.length > 0) {
                ParseArgs parseArgs = new ParseArgs(params);
                parseArgs.parseArgs(args);
            }
        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;
        int nRealPoints = params.nRealPoints;
        double sigmaCH2 = 3.95;
        double sigmaCH3 = 3.75;
        double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;
        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        PotentialGroup pTargetGroup = new PotentialGroup(2);
        System.out.println("Siepman "+nSpheres+"-mer chains B"+nPoints+" flexible correction at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);

        ClusterAbstract targetCluster = null;
        
        ClusterAbstract refCluster = null;
        
        double refIntegral = 0;

        if (nPoints == 3) {
            // 0 - (1,3) - 2
            ClusterBonds bonds1a = new ClusterBonds(4, new int[][][]{{{0,1},{1,2}}});
            ClusterBonds bonds1b = new ClusterBonds(4, new int[][][]{{{0,3},{2,3}}});
            ClusterBonds bonds1as = new ClusterBonds(4, new int[][][]{{{0,1},{2,3}}});
            ClusterBonds bonds1bs = new ClusterBonds(4, new int[][][]{{{0,3},{1,2}}});
            // coefficient is 0.5 (because of the permutations)
            targetCluster = new ClusterSum(new ClusterBonds[]{bonds1a, bonds1b, bonds1as, bonds1bs}, new double[]{1.0/2.0, 1.0/2.0, -1.0/2.0, -1.0/2.0}, new MayerFunction[]{fTarget});
            refCluster = new ClusterSum(new ClusterBonds[]{bonds1a}, new double[]{1}, new MayerFunction[]{fRef});
            refIntegral = 4*HSB[2]*HSB[2];
        }
        else if (nPoints == 4) {
            // -3/2     T4 = I3I2 - I3*I2
            double w1 = -3.0/2.0/2.0;
            ClusterBonds bonds1a = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{1,3},{2,3}}});
            ClusterBonds bonds1b = new ClusterBonds(5, new int[][][]{{{0,4},{4,2},{4,3},{2,3}}});
            ClusterBonds bonds1as = new ClusterBonds(5, new int[][][]{{{0,1},{4,2},{4,3},{2,3}}});
            ClusterBonds bonds1bs = new ClusterBonds(5, new int[][][]{{{0,4},{1,2},{1,3},{2,3}}});
            double ref1 = 3*HSB[3]*2*HSB[2];

            // -3/2     T3 = C4 - I3*I2
            //   C4 = 0-1-2-3   4-membered chain
            double w2 = -3.0/2.0/2.0;
            ClusterBonds bonds2a = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{2,3}}});
            ClusterBonds bonds2b = new ClusterBonds(5, new int[][][]{{{0,4},{4,2},{2,3}}});
            ClusterBonds bonds2as = new ClusterBonds(5, new int[][][]{{{0,4},{1,2},{2,3}}});
            ClusterBonds bonds2bs = new ClusterBonds(5, new int[][][]{{{0,1},{4,2},{2,3}}});
            double ref2 = -8*HSB[2]*HSB[2]*HSB[2];

            // -1/2     T2 = S3 - I3*I2
            //   S3 = 1-{0,2,3}  3-arm star
            double w3 = -1.0/2.0/2.0;
            ClusterBonds bonds3a = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{1,3}}});
            ClusterBonds bonds3b = new ClusterBonds(5, new int[][][]{{{0,4},{4,2},{4,3}}});
            ClusterBonds bonds3as = new ClusterBonds(5, new int[][][]{{{0,1},{4,2},{4,3}}});
            ClusterBonds bonds3bs = new ClusterBonds(5, new int[][][]{{{0,4},{1,2},{1,3}}});
            ClusterBonds bonds3cs = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{4,3}}});
            ClusterBonds bonds3ds = new ClusterBonds(5, new int[][][]{{{0,4},{4,2},{1,3}}});
            ClusterBonds bonds3es = new ClusterBonds(5, new int[][][]{{{0,4},{4,2},{1,3}}});
            ClusterBonds bonds3fs = new ClusterBonds(5, new int[][][]{{{0,1},{1,2},{4,3}}});
            double ref3 = -8*HSB[2]*HSB[2]*HSB[2];

            if (nRealPoints == 5) {
                System.out.println("4-body contribution to B5");
                // we want the I2*foo contribution to B5.  This is the foo calculation.
                // foo = (3 T4 + 6 T3 + 1/2 T2)
                w1 = 3.0/2.0;
                w2 = 6.0/2.0;
                w3 = 1.0/2.0/2.0;
            }
            
            targetCluster = new ClusterSum(new ClusterBonds[]{
                    bonds1a, bonds1b, bonds1as, bonds1bs, bonds2a, bonds2b, bonds2as, bonds2bs,
                    bonds3a, bonds3b, bonds3as, bonds3bs, bonds3cs, bonds3ds, bonds3es, bonds3fs},
                    new double[]{w1,w1,-w1,-w1, w2,w2,-w2,-w2, w3,w3,-w3/3,-w3/3,-w3/3,-w3/3,-w3/3,-w3/3},
                    new MayerFunction[]{fTarget});
            refCluster = new ClusterSum(new ClusterBonds[]{bonds1a, bonds2a, bonds3a},
                    new double[]{Math.abs(w1),-Math.abs(w2),-Math.abs(w3)}, new MayerFunction[]{fRef});
            refIntegral = Math.abs(w1)*ref1 - Math.abs(w2)*ref2 - Math.abs(w3)*ref3;
        }
        else if (nPoints == 5) {
            ClusterBonds[] allBonds = new ClusterBonds[0];
            double[] allWeights = new double[0];
            
            
            // -2/3      I2F4      (subtracting I2 * F4)
            double w = -2.0/3.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{1,4},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{5,4},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{1,4},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{5,4},{2,4},{3,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});
            
            // -2        I2(R4D)   4-ring+diagonal, articulation point at less-bonded corner of ring
            //                     (subtracting I2 * R4D)
            w = -2.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{2,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{2,4},{3,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});
            
            // -2        I2(R4D)   4-ring+diagonal, articulation point at more-bonded corner of ring
            //                     (subtracting I2 * R4D)
            w = -2.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{1,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{5,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{1,4},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{5,4},{3,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});

            // -1/2        I3I3   (subtracting I3 * I3)
            w = -1.0/2.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{0,4},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{0,4},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{0,4},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{0,4},{1,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});

            // -2          I3C3   (subtracting I3I2 * I2)
            w = -2/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{4,2},{4,3},{2,3},{0,1},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{4,2},{4,3},{2,3},{0,5},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{4,2},{4,3},{2,3},{0,5},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{4,2},{4,3},{2,3},{0,1},{5,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});

            // -2          I3I2I2  (I2 attached at different points of I3)
            //             (subtracting I3I2 * I2)
            w = -2/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{1,2},{1,3},{2,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{5,2},{5,3},{2,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{1,2},{1,3},{2,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{5,2},{5,3},{2,3},{1,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});
            
            // -2      I2(R4)   4-ring  (subtracting I2 * R4)
            w = -2.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{1,2},{0,3},{1,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{5,2},{0,3},{5,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{1,2},{0,3},{1,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,2},{5,2},{0,3},{5,3},{1,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});

            // -1          I3I2I2  (both I2 attached at the same point of I3)
            //             (subtracting I3I2 * I2)
            w = -1.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{5,4}}}));
            
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{2,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{2,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{2,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{2,3},{5,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w/3, -w/3, -w/3, -w/3, -w/3, -w/3});

            // 3 I3I2 * I2       => 3 I2 * T4
            // -3 I3 * I2 * I2

            // 3 I3 * C3         => 3 I3 * T1
            // -3 I3 * I2 * I2

            // -2        C5      (subtracting C4 * I2)
            w = -2/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{2,3},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{2,3},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{2,3},{3,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{2,3},{3,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});
            
            // -2        C4I2 (I2 attached at a middle point of C4)
            //                (subtracting S3 * I2)
            w = -2.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{2,3},{2,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{2,3},{2,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{2,3},{2,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{2,3},{2,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w, -w});
            
            // -1/6      S4       (subtracting S3 * I2)
            w = -1.0/6.0/2.0;
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{5,4}}}));
            
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{1,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{1,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{5,3},{1,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{1,2},{1,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,1},{5,2},{5,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{1,2},{5,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{1,3},{5,4}}}));
            allBonds = (ClusterBonds[])Arrays.addObject(allBonds, new ClusterBonds(6, new int[][][]{{{0,5},{5,2},{5,3},{1,4}}}));
            allWeights = append(allWeights, new double[]{w, w, -w/4, -w/4, -w/4, -w/4, -w/4, -w/4, -w/4, -w/4});
            
            // 6    C4 * I2
            // 1/2  S3 * I2
            // 9/2  C3 * C3

            // -18  C3 * I2 * I2
            // 7    I2 * I2 * I2 * I2

            // 6    I2 * (C4 - C3 * I2) = 6 I2 * T3
            // 1/2  I2 * (S3 - C3 * I2) = 1/2 I2 * T2
            // 9/2  C3 * (C3 - I2 * I2) = 9/2 (T1 + I2 * I2) * T1
            //                          = 9/2 (T1 * T1 + I2 * I2 * T1)
            // -7   I2 * I2 * (C3 - I2 * I2) = -7 I2 * I2 * T1
            
            // analytic:
            // 3 I3*T1 + 9/2 T1^2 + 9/2 T1*I2^2 - 7 T1*I2^2
            //   = 3 I3*T1 + 9/2 T1^2 - 5/2 T1*I2^2
            // re-simulate
            // 3 I2*T4 + 6 I2*T3 + 1/2 I2*T2
            //  = I2 * (3 T4 + 6 T3 + 1/2 T2)

            targetCluster = new ClusterSum(allBonds, allWeights, new MayerFunction[]{fTarget});
            
            refCluster = Standard.virialCluster(5, fRef, true, eRef, false);
            refIntegral = HSB[5];
        }
        else {
            throw new RuntimeException("not yet");
        }
        
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // overerr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), ClusterWeightAbs.makeWeightCluster(targetCluster)};
        
        SpeciesAlkane species = new SpeciesAlkane(space, nSpheres);
        
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{species},
                new int[]{nPoints+1},temperature, new ClusterAbstract[]{refCluster, targetCluster}, sampleClusters, true);

        int[] constraintMap = new int[nPoints+1];
        for (int i=0; i<nPoints; i++) {
            constraintMap[i] = i;
        }
        constraintMap[nPoints] = 1;
        ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
        ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
        ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
        ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);
        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }

        AtomType typeCH3 = species.getAtomType(0);
        AtomType typeCH2 = species.getAtomType(1);
        pTargetGroup.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));
        // CH2 on molecule1 to CH3 on molecule2
        pTargetGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));
        pTargetGroup.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH2}));
        pTargetGroup.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));
        
        sim.integratorOS.setNumSubSteps(1000);

        // create the intramolecular potential here, add to it and add it to
        // the potential master if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        if (nSpheres > 2) {
            P3BondAngle p3 = new P3BondAngle(space);
            p3.setAngle(Math.PI*114.0/180.0);
            p3.setEpsilon(Kelvin.UNIT.toSim(62500));
            int[][] triplets = new int[nSpheres-2][3];
            for (int i=0; i<nSpheres-2; i++) {
                triplets[i][0] = i;
                triplets[i][1] = i+1;
                triplets[i][2] = i+2;
            }
            pIntra.addPotential(p3, new Atomset3IteratorIndexList(triplets));
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(0)});
        }
        MCMoveClusterTorsionMulti[] torsionMoves = null;
        if (nSpheres > 3) {
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres-3][4];
            for (int i=0; i<nSpheres-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            pIntra.addPotential(p4, new Atomset4IteratorIndexList(quads));
            torsionMoves = new MCMoveClusterTorsionMulti[2];
            torsionMoves[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves[0].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMoves[0]);
            torsionMoves[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves[1].setTemperature(temperature);
            sim.integrators[1].getMoveManager().addMCMove(torsionMoves[1]);
        }
        if (nSpheres > 4) {
            pIntra.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres-1}}));
        }
        if (nSpheres > 5) {
            int[][] pairs = new int[2*(nSpheres-5)][2];
            for (int i=0; i<nSpheres-5; i++) {
                pairs[2*i][0] = 0;
                pairs[2*i][1] = nSpheres-2-i;
                pairs[2*i+1][0] = nSpheres-1;
                pairs[2*i+1][1] = i+1;
            }
            pIntra.addPotential(p2CH2CH3,new ApiIndexList(pairs));
        }
        if (nSpheres > 6) {
            int[][] pairs = new int[(nSpheres-6)*(nSpheres-5)/2][2];
            int k = 0;
            for (int i=1; i<nSpheres-5; i++) {
                for (int j=i+4; j<nSpheres-1; j++) {
                    pairs[k][0] = i;
                    pairs[k][1] = j;
                    k++;
                }
            }
            pIntra.addPotential(p2CH2,new ApiIndexList(pairs));
        }

        for (int j=0; j<10000 && sim.box[1].getSampleCluster().value(sim.box[1]) < 1e-10; j++) {
            sim.integrators[1].doStep();
        }
        if (sim.box[1].getSampleCluster().value(sim.box[1]) == 0) {
            throw new RuntimeException("could not find a configuration for target system");
        }
        sim.accumulators[1].reset();
        
        if (false) {
            double size = (nSpheres+5)*1.5;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0/size));
            displayBox1.setPixelUnit(new Pixel(300.0/size));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
            
            
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
            diameterManager.setDiameter(typeCH2, 0.5*sigmaCH2);
            diameterManager.setDiameter(typeCH3, 0.5*sigmaCH3);
            displayBox1.setDiameterHash(diameterManager);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
            displayBox0.setColorScheme(colorScheme);
            colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
            displayBox1.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 10);
                    sim.equilibrate(null, 20);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }
            
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("B"+nPoints+" (L/mol)^"+(nPoints-1));
            final JPanel panelParentGroup = new JPanel(new java.awt.BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), java.awt.BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), java.awt.BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());
            
            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();
                
                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
            };
            IEtomicaDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
            Unit unit = new CompoundUnit(new Unit[]{new UnitRatio(Liter.UNIT, Mole.UNIT)}, new double[]{nPoints-1});
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            averageBox.setUnit(unit);
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            errorBox.setUnit(unit);
            sim.integratorOS.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));
            
            return;
        }
        
        // if running interactively, don't use the file
        String refFileName = isCommandline ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);
        
        sim.setAccumulatorBlockSize((int)steps);
        
        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "
                +sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "
                +sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));
        if (nSpheres > 3) {
            System.out.println("Torsion move acceptance "+torsionMoves[0].getTracker().acceptanceRatio()+" "+
                    torsionMoves[1].getTracker().acceptanceRatio());
        }

        if (false) {
            final double refIntegralF = refIntegral;
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    System.out.println("abs average: "+ratio*refIntegralF+", error: "+error*refIntegralF);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        
        sim.printResults(refIntegral);
	}
    
    public static ClusterBonds[] append(ClusterBonds[] inArray, ClusterBonds[] newBonds) {
        ClusterBonds[] outArray = new ClusterBonds[inArray.length + newBonds.length];
        System.arraycopy(inArray, 0, outArray, 0, inArray.length);
        System.arraycopy(newBonds, 0, outArray, inArray.length, newBonds.length);
        return outArray;
    }

    public static double[] append(double[] inArray, double[] newWeights) {
        double[] outArray = new double[inArray.length + newWeights.length];
        System.arraycopy(inArray, 0, outArray, 0, inArray.length);
        System.arraycopy(newWeights, 0, outArray, inArray.length, newWeights.length);
        return outArray;
    }

    /**
     * Inner class for parameters
     */
    public static class VirialSiepmannSpheresParam extends ParameterBase {
        public int nPoints = 3;
        public int nRealPoints = 3;
        public int nSpheres = 4;
        public double temperature = 300.0;   // Kelvin
        public long numSteps = 10000;
        public double refFreq = -1;
        public boolean flex = true;
    }
}
