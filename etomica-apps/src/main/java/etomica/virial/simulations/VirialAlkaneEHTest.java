/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;


import etomica.AlkaneEH.SpeciesAlkaneEH;
import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.*;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.IntegratorListenerAction;
import etomica.molecule.MoleculePositionGeometricCenterAlkaneEH;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.*;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.DimensionRatio;
import etomica.units.dimensions.Quantity;
import etomica.units.dimensions.Volume;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;

import javax.swing.*;
import java.awt.*;
import java.util.Map;
import java.util.Set;


/**
 * Mayer sampling simulation for alkanes using the TraPPE-Explicit Hydrogen force field.
 * M.G. Martin and J.I. Siepmann, "Transferable Potentials for Phase Equilibria. 3. Explicit-Hydrogen Description of Normal Alkanes, 1999
 *   
 * @author shu
 * 01-31-2013
 */

public class VirialAlkaneEHTest {

    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagrams flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            byte nc = gs.nodeCount();
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                str += " "+gs.nodeCount()+"bc";
            }
            else {
                str += " "+gs.getStore().toNumberString();
                if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
            first = false;
        }
        return str;
    }
	
	public static void main(String[] args) {

        VirialParam params = new VirialParam();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
        } else {

        }
        final int nPoints = params.nPoints;
        int nSpheres = params.nSpheres;// number of carbons
        double temperature = params.temperature;
        long steps = params.numSteps;
        double refFreq = params.refFreq;

        double sigmaH = 3.31; 
        double sigmaCH2 = 3.65;//C in CH2
        double sigmaCH3 = 3.30;//C in CH3
        double sigmaCH2H = 0.5 * (sigmaH + sigmaCH2);
        double sigmaCH3H = 0.5 * (sigmaH + sigmaCH3);
        double sigmaCH2CH3 = 0.5 * (sigmaCH2 + sigmaCH3);

        double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;////?????????

        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);

        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
//        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        PotentialGroup pTargetGroup = new PotentialGroup(2);
        
        if (nSpheres < 2){
        	throw new RuntimeException("invalid nSpheres parameter, you mean methane?");

        }
        System.out.println("Test, TraPPE-EH, "+nSpheres+"-mer chains B"+nPoints+" flexible correction at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        
        double epsilonH = Kelvin.UNIT.toSim(15.3);
        double epsilonCH2 = Kelvin.UNIT.toSim(5.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(4.0);
        double epsilonCH2H = Math.sqrt(epsilonH*epsilonCH2);
        double epsilonCH3H = Math.sqrt(epsilonH*epsilonCH3);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
                             
        P2LennardJones pHH = new P2LennardJones(space, sigmaH, epsilonH);// H-H
        P2LennardJones pCH2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);//CH2-CH2
        P2LennardJones pCH3CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);//CH3-CH3
        P2LennardJones pCH2H = new P2LennardJones(space, sigmaCH2H, epsilonCH2H);//H-CH2
        P2LennardJones pCH3H = new P2LennardJones(space, sigmaCH3H, epsilonCH3H);//H-CH3
        P2LennardJones pCH2CH3 = new P2LennardJones(space, sigmaCH2CH3, epsilonCH2CH3);//CH2-CH3

        MayerGeneral fTarget = new MayerGeneral(pTargetGroup);

        boolean alkaneFlex = nSpheres > 1 && nPoints > 2;
        VirialDiagrams alkaneDiagrams = new VirialDiagrams(nPoints, false, alkaneFlex);
        alkaneDiagrams.setDoReeHoover(false);
        ClusterSum targetCluster = alkaneDiagrams.makeVirialCluster(fTarget);
        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
        
        double refIntegral = HSB[nPoints];

        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];

        int[] targetDiagramNumbers = new int[0];
        
        targetDiagrams = alkaneDiagrams.makeSingleVirialClusters(targetCluster, null, fTarget);
        targetDiagramNumbers = new int[targetDiagrams.length];
        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = alkaneDiagrams.getMSMCGraphs(true, false);
        Map<Graph,Graph> cancelMap = alkaneDiagrams.getCancelMap();
        int iGraph = 0;
        boolean[] diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
        	System.out.print(iGraph+" ("+g.coefficient()+") "+g.getStore().toNumberString()); // toNumberString: its corresponding number
        	targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());

            Graph cancelGraph = cancelMap.get(g);
            if (cancelGraph != null) {
                diagramFlexCorrection[iGraph] = true;
                String gnStr = cancelGraph.getStore().toNumberString();
                Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(cancelGraph);

                System.out.print(" - "+getSplitGraphString(gSplit, alkaneDiagrams, false));

            }
        	System.out.println();
        	iGraph++;
        }
        System.out.println();
        Set<Graph> disconnectedGraphs = alkaneDiagrams.getExtraDisconnectedVirialGraphs();
        if (disconnectedGraphs.size() > 0) {
        	System.out.println("extra clusters:");
            
        	for (Graph g : disconnectedGraphs) {
        		Set<Graph> gSplit = alkaneDiagrams.getSplitDisconnectedVirialGraphs(g);
        		System.out.println(g.coefficient()+" "+getSplitGraphString(gSplit, alkaneDiagrams, true));
        	}
        	System.out.println();
        }
        
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }

        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // eovererr expects this string, BnHS
        System.out.println("B"+nPoints+"HS: "+refIntegral);
        
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), ClusterWeightAbs.makeWeightCluster(targetCluster)};

        SpeciesAlkaneEH species = new SpeciesAlkaneEH(space, nSpheres);
        
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{species},
                new int[]{alkaneFlex ?(nPoints+1):nPoints },temperature, new ClusterAbstract[]{refCluster, targetCluster},targetDiagrams, sampleClusters, false);
        
        if (alkaneFlex) {
            int[] constraintMap = new int[nPoints+1];
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }

            constraintMap[nPoints] = 0;
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);

        }

//        ((MCMoveStepTracker)sim.mcMoveTranslate[0].getTracker()).setNoisyAdjustment(true);
//        ((MCMoveStepTracker)sim.mcMoveTranslate[1].getTracker()).setNoisyAdjustment(true);

        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }
        
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;


        AtomType typeCH3 = species.getC_3Type();// C in CH3
        AtomType typeCH2 = species.getC_2Type();// C in CH2
        AtomType typeH = species.getHType();  // H

        pTargetGroup.addPotential(pHH, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeH}));//H-H
        pTargetGroup.addPotential(pCH2H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH2}));//H-CH2
        pTargetGroup.addPotential(pCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH3}));//H-CH3
        pTargetGroup.addPotential(pCH2H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH2}));//CH2-H
        pTargetGroup.addPotential(pCH2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));//CH2-CH2
        pTargetGroup.addPotential(pCH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));//CH2-CH3
        pTargetGroup.addPotential(pCH3H, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeH, typeCH3}));//CH3-H
        pTargetGroup.addPotential(pCH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));//CH3-CH2
        pTargetGroup.addPotential(pCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));//CH3-CH3

        // create the intramolecular potential here, add to it and add it to the "potential master" if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        
        //*********************************** torsional potential (XCCH) *********************************//
       
        double torsionHCCH_c3 = 717; //u(HCCH) = C3 * [1 - cos(3phi)], c3/kB = 717K
    	double torsionCCCH_c3 = 854;//u(CCCH) = C3 * [1 - cos(3phi)], c3/kB = 854K

    	P4BondTorsionAlkaneXCCH p4HCCH = new P4BondTorsionAlkaneXCCH(space, 0.0, 0.0, 0.0, Kelvin.UNIT.toSim(torsionHCCH_c3));
    	P4BondTorsionAlkaneXCCH p4CCCH = new P4BondTorsionAlkaneXCCH(space, 0.0, 0.0, 0.0, Kelvin.UNIT.toSim(torsionCCCH_c3));
    	
        //*********************************** Set Geometric center *********************************//
    	// geometric center is based on all carbons, hydrogens not included
    	MoleculePositionGeometricCenterAlkaneEH center = new MoleculePositionGeometricCenterAlkaneEH(space,species);
    	((MCMoveRotateMolecule3D)sim.mcMoveRotate[0]).setPositionDefinition(center);
    	((MCMoveRotateMolecule3D)sim.mcMoveRotate[1]).setPositionDefinition(center);
    	((CoordinatePairMoleculeSet)sim.box[0].getCPairSet()).setPositionDefinition(center);
    	((CoordinatePairMoleculeSet)sim.box[1].getCPairSet()).setPositionDefinition(center);
    	
        //*********************************** Add moves *********************************//
        MCMoveClusterRotateCH3[] rotateCH3Move = new MCMoveClusterRotateCH3[2];
        MCMoveClusterWiggleAlkaneEH[] wiggleMove = new MCMoveClusterWiggleAlkaneEH[2];
        MCMoveClusterTorsionAlkaneEH[] torsionMove = null;

        if (nSpheres == 2){// C2H6, only H-C-C-H torsion potential, 9 sets in total 
            System.out.println("C2H6, add H--C(in CH3)--C(in CH3)--H torsion intramolecular potential !");  
            
        	int[][] quadsHCCH = new int[9][0];
        	quadsHCCH[0]= new int[]{2,0,1,3};
        	quadsHCCH[1]= new int[]{2,0,1,5};
        	quadsHCCH[2]= new int[]{2,0,1,7};
        	quadsHCCH[3]= new int[]{4,0,1,3};
        	quadsHCCH[4]= new int[]{4,0,1,5};
        	quadsHCCH[5]= new int[]{4,0,1,7};
        	quadsHCCH[6]= new int[]{6,0,1,3};
        	quadsHCCH[7]= new int[]{6,0,1,5};
        	quadsHCCH[8]= new int[]{6,0,1,7};
        	
           	pIntra.addPotential(p4HCCH, new Atomset4IteratorIndexList(quadsHCCH));
            
        }
        // integrators share a common potentialMaster. just add to one
        sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(0)});///// CALL ONCE IS FINE!
        

        /////////////////////////////////// rotateCH3 ////////////////////////////////
        /////////////////////////////////// rotateCH3 ////////////////////////////////           
        System.out.println("----- add rotateCH3 MC move for all types of alkanes-----");
        rotateCH3Move[0] = new MCMoveClusterRotateCH3(sim, sim.integrators[1].getPotentialMaster(), targetCluster.pointCount(), space);
        sim.integrators[0].getMoveManager().addMCMove(rotateCH3Move[0]);//reference system
        rotateCH3Move[1] = new MCMoveClusterRotateCH3(sim, sim.integrators[1].getPotentialMaster(), targetCluster.pointCount(), space);
        sim.integrators[1].getMoveManager().addMCMove(rotateCH3Move[1]);//target system
        
        if (nSpheres > 2) {
        	//====================== bending (triplets)==========================================//
            System.out.println("n>2, add CCC bending intramolecular potential ");  
        	double bendAngle = Math.PI * 112.7 / 180.0;//degree ->radian
        	double bendEpsilon = Kelvin.UNIT.toSim(58765.0);//58765 in K/rad^2
            P3BondAngle p3 = new P3BondAngle(space);
            p3.setAngle(bendAngle);
            p3.setEpsilon(bendEpsilon);
            int[][] triplets = new int[nSpheres-2][3];// all carbon atoms, exclude hydrogens
            for (int i=0; i<nSpheres-2; i++) {
            	triplets[i][0] = i;
            	triplets[i][1] = i+1;
            	triplets[i][2] = i+2;
            }
            pIntra.addPotential(p3, new Atomset3IteratorIndexList(triplets));

        	//====================== XCCH torsion ========================================== //
            int[][] quadsCCCH  = new int[6][0];// H-C0-C1-H  & H-C(n-2)-C(n-1)-H, 6*2 = 12 sets ===> HCCH potential
            int[][] quadsHCCH1 = new int[6][0];// H-C0-C1-C2 potential, 3*2 = 6 sets   ======> CCCH potential
            int[][] quadsHCCH2 = new int[6][0];// C(n-3)-C(n-2)-C(n-1)-H potential, 3*2 = 6 sets ======> CCCH potential
            
            quadsCCCH[0]= new int[]{nSpheres,  0,1,2};
            quadsCCCH[1]= new int[]{nSpheres*2,0,1,2};
            quadsCCCH[2]= new int[]{nSpheres*3,0,1,2};
            quadsCCCH[3]= new int[]{nSpheres*2-1, nSpheres-1,nSpheres-2, nSpheres-3};
            quadsCCCH[4]= new int[]{nSpheres*3-1, nSpheres-1,nSpheres-2, nSpheres-3};
            quadsCCCH[5]= new int[]{nSpheres*3+1, nSpheres-1,nSpheres-2, nSpheres-3};

            quadsHCCH1[0]= new int[]{nSpheres,   0,1, nSpheres+1};
            quadsHCCH1[1]= new int[]{nSpheres,   0,1, nSpheres*2+1};
            quadsHCCH1[2]= new int[]{nSpheres*2, 0,1, nSpheres+1};
            quadsHCCH1[3]= new int[]{nSpheres*2, 0,1, nSpheres*2+1};
            quadsHCCH1[4]= new int[]{nSpheres*3, 0,1, nSpheres+1};
            quadsHCCH1[5]= new int[]{nSpheres*3, 0,1, nSpheres*2+1};

            quadsHCCH2[0]= new int[]{nSpheres*2-2, nSpheres-2,nSpheres-1, nSpheres*2-1};
            quadsHCCH2[1]= new int[]{nSpheres*2-2, nSpheres-2,nSpheres-1, nSpheres*3-1};
            quadsHCCH2[2]= new int[]{nSpheres*2-2, nSpheres-2,nSpheres-1, nSpheres*3+1};
            quadsHCCH2[3]= new int[]{nSpheres*3-2, nSpheres-2,nSpheres-1, nSpheres*2-1};
            quadsHCCH2[4]= new int[]{nSpheres*3-2, nSpheres-2,nSpheres-1, nSpheres*3-1};
            quadsHCCH2[5]= new int[]{nSpheres*3-2, nSpheres-2,nSpheres-1, nSpheres*3+1};
            
            
            System.out.println("n>2, add HCCH torsion intramolecular potential !");  
            System.out.println("n>2, add CCCH torsion intramolecular potential !");  
            pIntra.addPotential(p4CCCH, new Atomset4IteratorIndexList(quadsCCCH));
            pIntra.addPotential(p4HCCH, new Atomset4IteratorIndexList(quadsHCCH1));
            pIntra.addPotential(p4HCCH, new Atomset4IteratorIndexList(quadsHCCH2));
            
            /////////////////////////////////// wiggle ////////////////////////////////
            /////////////////////////////////// wiggle ////////////////////////////////           
            System.out.println("n>2, add wiggle MC move!");  
            wiggleMove[0] = new MCMoveClusterWiggleAlkaneEH(sim,sim.integrators[1].getPotentialMaster(),  targetCluster.pointCount(), space);
            sim.integrators[0].getMoveManager().addMCMove(wiggleMove[0]);//reference system
            wiggleMove[1] = new MCMoveClusterWiggleAlkaneEH(sim,sim.integrators[1].getPotentialMaster(),  targetCluster.pointCount(), space);
            sim.integrators[1].getMoveManager().addMCMove(wiggleMove[1]);//target system
        }
          
        if (nSpheres > 3) {
            //u(CCCC) = c1*[1+cos(phi)] + c2*[1-cos(2*phi)] + C3 * [1+cos(3phi)], c1/kB = 355.03 K, c2/kB = -68.19K, c3/kB = 791.32K
        	P4BondTorsion p4CCCC = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quadsCCCC = new int[nSpheres-3][4];
            for (int i=0; i<nSpheres-3; i++) {
            	quadsCCCC[i][0] = i;
            	quadsCCCC[i][1] = i+1;
            	quadsCCCC[i][2] = i+2;
            	quadsCCCC[i][3] = i+3;
            }
            System.out.println("n>3, add CCCC torsion intramolecular potential !");  
            pIntra.addPotential(p4CCCC, new Atomset4IteratorIndexList(quadsCCCC));
            
            /////////////////////////////////// torsion ////////////////////////////////
            /////////////////////////////////// torsion ////////////////////////////////           

            System.out.println("n>3, add torsion MC move");
            torsionMove = new MCMoveClusterTorsionAlkaneEH[2];
            torsionMove[0] = new MCMoveClusterTorsionAlkaneEH(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4CCCC, 400);//reference
            torsionMove[0].setTemperature(temperature);
            torsionMove[1] = new MCMoveClusterTorsionAlkaneEH(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4CCCC, 400);//target
            torsionMove[1].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMove[0]);
            sim.integrators[1].getMoveManager().addMCMove(torsionMove[1]);

        }

        if (nSpheres > 4) { 
            System.out.println("n>4, add H(of CH3)--H(of CH3), C(of CH3)--C(of CH3), C(of CH3)--H(of CH3) pair intramolecular potential !");  
        	// ========================= Indices of atoms in CH3 groups ============================//
        	int h01 = nSpheres;// H1 on left C(H3)
            int h02 = nSpheres * 2;// H2 on left C(H3)
            int h03 = nSpheres * 3;// H3 on left C(H3)
            int h11 = nSpheres * 2 - 1;// H1 on right C(H3)
            int h12 = nSpheres * 3 - 1;// H2 on right C(H3)
            int h13 = nSpheres * 3 + 1;// H3 on right C(H3)
            int c0  = 0;// carbon on the left
            int c1  = nSpheres -1;// carbon on the right
            int[][] pairsCH = new int[6][2];
            // 1. H (of CH3) and C (of end CH3)
            pairsCH[0][0] = h01;
            pairsCH[0][1] = c1;
            pairsCH[1][0] = h02;
            pairsCH[1][1] = c1;
            pairsCH[2][0] = h03;
            pairsCH[2][1] = c1;
            // 2. C (of CH3) and H (of end CH3)
            pairsCH[3][0] = c0;
            pairsCH[3][1] = h11;
            pairsCH[4][0] = c0;
            pairsCH[4][1] = h12;
            pairsCH[5][0] = c0;
            pairsCH[5][1] = h13;
            
            // 3. H (of CH3) and H (of end CH3)
            int[][] pairsHH = new int[9][2];
            pairsHH[0][0] = h01;
            pairsHH[0][1] = h11;
            pairsHH[1][0] = h01;
            pairsHH[1][1] = h12;
            pairsHH[2][0] = h01;
            pairsHH[2][1] = h13;
                     
            pairsHH[3][0] = h02;
            pairsHH[3][1] = h11;
            pairsHH[4][0] = h02;
            pairsHH[4][1] = h12;
            pairsHH[5][0] = h02;
            pairsHH[5][1] = h13;            
            
            pairsHH[6][0] = h03;
            pairsHH[6][1] = h11;            
            pairsHH[7][0] = h03;
            pairsHH[7][1] = h12;            
            pairsHH[8][0] = h03;
            pairsHH[8][1] = h13;  
            
            pIntra.addPotential(pCH3CH3,new ApiIndexList(new int[][]{{0,nSpheres-1}}));// C(of CH3)-C(of CH3) potential
            pIntra.addPotential(pCH3H,new ApiIndexList(pairsCH));// H(of CH3)-C(of CH3) potential
            pIntra.addPotential(pHH,new ApiIndexList(pairsHH));// H(of CH3)-H(of CH3) potential       
            
        }

        
        if (nSpheres > 5) {
        	
            System.out.println("n>5, add CH3-CH2, HH, CH2-H, CH3-H pair intramolecular potential !");  
            
            int[][] pairsCH2CH3  = new int[2*(nSpheres-5)][2];// add C(of CH3)-C(of CH2)
            int[][] pairsH_CH2_1 = new int[2*(nSpheres-5)][2];// add (C of CH2)--(H of CH3)
            int[][] pairsH_CH2_2 = new int[2*(nSpheres-5)][2];// add (C of CH2)--(H of CH3)
            int[][] pairsH_CH2_3 = new int[2*(nSpheres-5)][2];// add (C of CH2)--(H of CH3)
            
            int[][] pairsH_CH3_1 = new int[2*(nSpheres-5)][2];// add (C of CH3)--(H of CH2)
            int[][] pairsH_CH3_2 = new int[2*(nSpheres-5)][2];// add (C of CH3)--(H of CH2)
            int[][] pairsH1H1    = new int[2*(nSpheres-5)][2];// add (H1 of CH3)--(H1 of CH2)
            int[][] pairsH1H2    = new int[2*(nSpheres-5)][2];// add (H1 of CH3)--(H2 of CH2)
            int[][] pairsH2H1    = new int[2*(nSpheres-5)][2];// add (H2 of CH3)--(H1 of CH2)
            int[][] pairsH2H2    = new int[2*(nSpheres-5)][2];// add (H2 of CH3)--(H2 of CH2)
            int[][] pairsH3H1    = new int[2*(nSpheres-5)][2];// add (H3 of CH3)--(H1 of CH2)
            int[][] pairsH3H2    = new int[2*(nSpheres-5)][2];// add (H3 of CH3)--(H2 of CH2)

            for (int i=0; i<nSpheres-5; i++) {
            	int cL = 0;// C(H3) on left
            	int cH2L = nSpheres - 2 - i; // C(H2) corresponding to left C(H3) 
            	int cR = nSpheres - 1; // C(H3) on right
            	int cH2R = i + 1;// C(H2) corresponding to right C(H3)
                int hL1 = nSpheres;         // H1 on left CH3
                int hL2 = nSpheres * 2 ;    // H2 on left CH3
                int hL3 = nSpheres * 3 ;    // H3 on left CH3
                int hR1 = nSpheres * 2 - 1 ;// H1 on right CH3
                int hR2 = nSpheres * 3 - 1 ;// H2 on right CH3
                int hR3 = nSpheres * 3 + 1 ;// H3 on right CH3

                int hL1CH2  = cH2L + nSpheres;// H1 on CH2 corresponding to left CH3
                int hL2CH2  = cH2L + nSpheres * 2 ;// H2 on CH2 corresponding to left CH3
                int hR1CH2  = cH2R + nSpheres ;// H1 on CH2 corresponding to right CH3
                int hR2CH2  = cH2R + nSpheres * 2 ;// H2 on CH2 corresponding to right CH3
                
            	/////////////////////////////////// (a) C(of CH3)--C(of CH2) potential ////////////////////////////////////////////
                pairsCH2CH3[2*i][0] = cL;// C(H3) on the left
                pairsCH2CH3[2*i][1] = cH2L;
                pairsCH2CH3[2*i+1][0] = cR;//C(H3) on the right
                pairsCH2CH3[2*i+1][1] = cH2R;
            	/////////////////////////////////// (b) H(on CH3)-C(on CH2) potential ////////////////////////////////////////////
                pairsH_CH2_1[2*i][0] = hL1 ;// H1 on left CH3
                pairsH_CH2_1[2*i][1] = cH2L ;// C(H2)
                pairsH_CH2_1[2*i+1][0] = hR1; // H1 on right CH3
                pairsH_CH2_1[2*i+1][1] = cH2R;// C(H2)

                pairsH_CH2_2[2*i][0] = hL2 ; //H2 on left CH3 
                pairsH_CH2_2[2*i][1] = cH2L ;
                pairsH_CH2_2[2*i+1][0] = hR2 ; //H2 on right CH3
                pairsH_CH2_2[2*i+1][1] = cH2R ;
                
                pairsH_CH2_3[2*i][0]   = hL3 ; // H2 on left CH3
                pairsH_CH2_3[2*i][1]   = cH2L ; 
                pairsH_CH2_3[2*i+1][0] = hR3 ; // H2 on left CH3
                pairsH_CH2_3[2*i+1][1] = cH2R ; 
            	/////////////////////////////////// (c) C(of CH3)-H(on CH2) potential////////////////////////////////////////////
                pairsH_CH3_1[2*i][0] = cL ; // C(H3) on the left
                pairsH_CH3_1[2*i][1] = hL1CH2; // H1 on CH2
                pairsH_CH3_1[2*i+1][0] = cR ; // C(H3) on the right
                pairsH_CH3_1[2*i+1][1] = hR1CH2; // H1 on CH2 
                
                pairsH_CH3_2[2*i][0] = cL ; // C(H3) on the left
                pairsH_CH3_2[2*i][1] = hL2CH2 ; // H2 on CH2
                pairsH_CH3_2[2*i+1][0] = cR ; // C(H3) on the right
                pairsH_CH3_2[2*i+1][1] = hR2CH2 ; // H2 on CH2 
            	/////////////////////////////////// (d) H(of CH3)-H(on CH2) potential ////////////////////////////////////////////
                pairsH1H1[2*i][0] = hL1 ; // H1 on left CH3
                pairsH1H1[2*i][1] = hL1CH2; // H1 on CH2
                pairsH1H1[2*i+1][0] =  hR1; //H1 on the right CH3
                pairsH1H1[2*i+1][1] = hR1CH2; // H1 on CH2 
                
                pairsH1H2[2*i][0] = hL1 ; // H1 on left CH3
                pairsH1H2[2*i][1] = hL2CH2; // H2 on CH2
                pairsH1H2[2*i+1][0] =  hR1; // H1 on the right CH3
                pairsH1H2[2*i+1][1] = hR2CH2; // H2 on CH2 
                
                pairsH2H1[2*i][0] = hL2 ; // H2 on left CH3
                pairsH2H1[2*i][1] = hL1CH2; // H1 on CH2
                pairsH2H1[2*i+1][0] =  hR2; // H2 on the right CH3
                pairsH2H1[2*i+1][1] = hR1CH2; // H1 on CH2 
                
                pairsH2H2[2*i][0] = hL2 ; // H2 on left CH3
                pairsH2H2[2*i][1] = hL2CH2; // H2 on CH2
                pairsH2H2[2*i+1][0] =  hR2; // H2 on the right CH3
                pairsH2H2[2*i+1][1] = hR2CH2; // H2 on CH2 
                
                
                pairsH3H1[2*i][0] = hL3 ; // H3 on left CH3
                pairsH3H1[2*i][1] = hL1CH2; // H1 on CH2
                pairsH3H1[2*i+1][0] =  hR3; // H3 on the right CH3
                pairsH3H1[2*i+1][1] = hR1CH2; // H1 on CH2 
                
                pairsH3H2[2*i][0] = hL3 ; // H3 on left CH3
                pairsH3H2[2*i][1] = hL2CH2; // H2 on CH2
                pairsH3H2[2*i+1][0] =  hR3; // H3 on the right CH3
                pairsH3H2[2*i+1][1] = hR2CH2; // H2 on CH2 
                
            }
         
            pIntra.addPotential(pCH2CH3,new ApiIndexList(pairsCH2CH3));
            pIntra.addPotential(pCH2H,new ApiIndexList(pairsH_CH2_1));
            pIntra.addPotential(pCH2H,new ApiIndexList(pairsH_CH2_2));
            pIntra.addPotential(pCH2H,new ApiIndexList(pairsH_CH2_3));
            pIntra.addPotential(pCH3H,new ApiIndexList(pairsH_CH3_1));
            pIntra.addPotential(pCH3H,new ApiIndexList(pairsH_CH3_2));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH1H1));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH1H2));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH2H1));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH2H2));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH3H1));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH3H2));

        }
        
             
        
        if (nSpheres > 6) {
            System.out.println("n>6, add CH2-CH2, H(of CH2)-H(of CH2) pair intramolecular potential !");  

        	// add C(of CH2)--C(of CH2) potential
           	int[][] pairs = new int[(nSpheres-6)*(nSpheres-5)/2][2];
            // add H(of CH2)--H(of CH2) potential
            int[][] pairsH1H1 = new int[(nSpheres-6)*(nSpheres-5)/2][2];// H1--H1
            int[][] pairsH1H2 = new int[(nSpheres-6)*(nSpheres-5)/2][2];// H1--H2
            int[][] pairsH2H1 = new int[(nSpheres-6)*(nSpheres-5)/2][2];// H2--H1
            int[][] pairsH2H2 = new int[(nSpheres-6)*(nSpheres-5)/2][2];// H2--H2

            int k = 0;
            for (int i=1; i<nSpheres-5; i++) {
                for (int j=i+4; j<nSpheres-1; j++) {
                    pairs[k][0] = i;// C(of CH2)
                    pairs[k][1] = j;// C(of the other CH2)
                
                    pairsH1H1[k][0] = i + nSpheres;
                    pairsH1H1[k][1] = j + nSpheres;
                    
                    pairsH1H2[k][0] = i + nSpheres;
                    pairsH1H2[k][1] = j + nSpheres * 2;
                    
                    pairsH2H1[k][0] = i + nSpheres * 2;
                    pairsH2H1[k][1] = j + nSpheres ;
                    
                    pairsH2H2[k][0] = i + nSpheres * 2;
                    pairsH2H2[k][1] = j + nSpheres * 2;
                   
                    k++;
                }
            }
            pIntra.addPotential(pCH2CH2,new ApiIndexList(pairs));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH1H1));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH1H2));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH2H1));
            pIntra.addPotential(pHH,new ApiIndexList(pairsH2H2));

        }
        
        if (false) {
        	
        	   double size = (nSpheres+5)*1.5;
               sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
               sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
               SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
               DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
               DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
               displayBox0.setPixelUnit(new Pixel(300.0/size));
               displayBox1.setPixelUnit(new Pixel(300.0/size));
               displayBox0.setShowBoundary(false);
               displayBox1.setShowBoundary(false);
               ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
               ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
               
               
               DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
               diameterManager.setDiameter(typeCH2, 0.2*sigmaCH2);
               diameterManager.setDiameter(typeCH3, 0.2*sigmaCH3);
               diameterManager.setDiameter(typeH, 0.1*sigmaH);

               displayBox1.setDiameterHash(diameterManager);
               ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[0], sim.getRandom());
              
               displayBox0.setColorScheme(colorScheme);
               colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box[1], sim.getRandom());
               displayBox1.setColorScheme(colorScheme);
                          
//            double size = (nSpheres+5) * 1.5;
//            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
//            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
//            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
//            simGraphic.getDisplayBox(sim.box[0]).setPixelUnit(new Pixel(300.0/size));
//            simGraphic.getDisplayBox(sim.box[1]).setPixelUnit(new Pixel(300.0/size));
//            simGraphic.getDisplayBox(sim.box[0]).setShowBoundary(false);
//            simGraphic.getDisplayBox(sim.box[1]).setShowBoundary(false);
// 
//            DiameterHashByType diameter = new DiameterHashByType(sim); 
//            diameter.setDiameter(species.getAtomType(0), 0.3*sigmaCH3);// C(H3)
//            diameter.setDiameter(species.getAtomType(1), 0.3*sigmaCH2);// C(H2)
//            diameter.setDiameter(species.getAtomType(2), 0.2*sigmaH);// H  
//            simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
//            simGraphic.getDisplayBox(sim.box[0]).setDiameterHash(diameter);
//        
//            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
//            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.red);// C(H3)
//            colorScheme.setColor(sim.getSpecies(0).getAtomType(1), Color.blue);// C(H2)          
//            colorScheme.setColor(sim.getSpecies(0).getAtomType(2), Color.cyan);// H
//            
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

            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
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
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        
        MeterVirial meterDiagrams = new MeterVirial(targetDiagrams);
        meterDiagrams.setBox(sim.box[1]);
        AccumulatorAverageCovariance accumulatorDiagrams = null;
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/20);

        
        System.out.println("equilibration finished");   
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
      
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize()+" "+sim.mcMoveRotate[0].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[0].getStepSize())));
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize()+" "+sim.mcMoveRotate[1].getStepSize()+" "
                +(sim.mcMoveWiggle==null ? "" : (""+sim.mcMoveWiggle[1].getStepSize())));
        
        System.out.println("RotateCH3 move acceptance "+ rotateCH3Move[0].getTracker().acceptanceRatio()+" "+
        		rotateCH3Move[1].getTracker().acceptanceRatio());
        
        if (nSpheres > 2) {
            System.out.println("Wiggle move acceptance "+ wiggleMove[0].getTracker().acceptanceRatio()+" "+
     	       		wiggleMove[1].getTracker().acceptanceRatio());
        	System.out.println("Wiggle move step sizes " + wiggleMove[0].getStepSize() + " "+	wiggleMove[1].getStepSize());
        }
        
        if (nSpheres > 3) {
        	System.out.println("Torsion move acceptance "+torsionMove[0].getTracker().acceptanceRatio()+" "+
    	                torsionMove[1].getTracker().acceptanceRatio());
        	System.out.println("Torsion move step size "+torsionMove[0].getStepSize() + " "+   torsionMove[1].getStepSize());
           
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
        sim.getController().actionPerformed();
        
        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());
        sim.printResults(refIntegral);
 
        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(AccumulatorAverage.AVERAGE.index);
        IData dataErr = allData.getData(AccumulatorAverage.ERROR.index);
        IData dataCov = allData.getData(AccumulatorAverageCovariance.BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        int nTotal = (targetDiagrams.length+2);
        double oVar = dataCov.getValue(nTotal*nTotal-1);
        for (int i=0; i<targetDiagrams.length; i++) {
            if (targetDiagramNumbers[i]<0) {
                System.out.print("diagram "+(-targetDiagramNumbers[i])+("bc "));
            }
            else {
                System.out.print("diagram "+targetDiagramNumbers[i]);
//                    if (fTargetDiagramNumbers[i] != 0) {
//                        System.out.print(fTargetDiagramNumbers[i]);
//                    }

                if (diagramFlexCorrection[i]) {
                    System.out.print("c");
                }
                System.out.print(" ");
            }
            // average is vi/|v| average, error is the uncertainty on that average
            // ocor is the correlation coefficient for the average and overlap values (vi/|v| and o/|v|)
            double ivar = dataCov.getValue((i+1)*nTotal+(i+1));
            double ocor = ivar*oVar == 0 ? 0 : dataCov.getValue(nTotal*(i+1)+nTotal-1)/Math.sqrt(ivar*oVar);
            System.out.print(String.format("average: %20.15e  error: %10.15e  ocor: %7.5f", dataAvg.getValue(i+1), dataErr.getValue(i+1), ocor));
            if (targetDiagrams.length > 1) {
                System.out.print("  dcor:");
                for (int j=0; j<targetDiagrams.length; j++) {
                    if (i==j) continue;
                    double jvar = dataCov.getValue((j+1)*nTotal+(j+1));
                    double dcor = ivar*jvar == 0 ? 0 : dataCov.getValue((i+1)*nTotal+(j+1))/Math.sqrt(ivar*jvar);
                    System.out.print(String.format(" %8.6f", dcor));
                }
            }
            System.out.println();
        }
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
    public static class VirialParam extends ParameterBase {
        public int nPoints = 3;
        public int nSpheres = 6;
        public double temperature = 1500;   // Kelvin
        public long numSteps = 1000000;
        public double refFreq = -1;
    }

}
