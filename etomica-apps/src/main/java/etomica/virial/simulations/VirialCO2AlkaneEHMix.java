/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import java.awt.Color;
import java.util.Map;
import java.util.Set;

import etomica.AlkaneEH.SpeciesAlkaneEH;
import etomica.action.IAction;
import etomica.atom.IAtomType;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.ISpecies;
import etomica.atom.AtomPositionGeometricCenterAlkaneEH;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.data.IData;
import etomica.data.types.DataGroup;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.potential.P2CO2TraPPE;
import etomica.potential.P2LennardJones;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.P4BondTorsionAlkaneXCCH;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.ClusterAbstract;
import etomica.virial.ClusterSum;
import etomica.virial.ClusterSumShell;
import etomica.virial.ClusterWeight;
import etomica.virial.ClusterWeightAbs;
import etomica.virial.CoordinatePairMoleculeSet;
import etomica.virial.MCMoveClusterMoleculeMulti;
import etomica.virial.MCMoveClusterRotateCH3;
import etomica.virial.MCMoveClusterRotateMoleculeMulti;
import etomica.virial.MCMoveClusterTorsionAlkaneEH;
import etomica.virial.MCMoveClusterWiggleAlkaneEH;
import etomica.virial.MayerFunction;
import etomica.virial.MayerGeneral;
import etomica.virial.MayerHardSphere;
import etomica.virial.SpeciesTraPPECO2;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.cluster.VirialDiagramsMix2;

/**
 *   Mayer sampling simulation for CO2(rigid, TraPPE)-alkanes(TraPPE-EH) mixture
 *   cross virial coefficients
 *   Using VirialDiagramMix2 to generate diagrams
 *   rigid diagrams and flex diagrams are calculated separately (Bij = Bij[flexID=-1] + Bij[flexID=1])
 *   
 *   @author shu
 *   May 2013
 * 
 */
public class VirialCO2AlkaneEHMix {
	
    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagramsMix2 flexDiagrams, boolean correction) {
        DeleteEdge edgeDeleter = new DeleteEdge();
        DeleteEdgeParameters ede = new DeleteEdgeParameters(flexDiagrams.eBond);
        boolean first = true;
        String str = "";
        for (Graph gs : gSplit) {
            // str +="_";//add a space between number and letters
            //loop over every node in the graph to get the node color,append this char to str
        	
////            for ( Node node : gs.nodes()){
////            	str += node.getColor();
////            }
            byte nc = gs.nodeCount();
            if (VirialDiagrams.graphHasEdgeColor(gs, flexDiagrams.efbcBond)) {
                str += " "+gs.nodeCount()+"bc";
            }
            else {
                str += " "+gs.getStore().toNumberString();
                if (flexDiagrams.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
//            System.out.println("first is:"+first);
//
//            System.out.println("correction is:"+correction);
            first = false;
            for ( Node node : gs.nodes()){
            	str += node.getColor();
            }
        }
        return str;
    }
	
	public static void main(String[] args) {
        VirialMixParam params = new VirialMixParam();
        if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
			
		}
        final int nPoints = params.nPoints;
        double temperature = params.temperature;
        long steps = params.numSteps;
        int nSpheres = params.nSpheres;// number of carbons
        int flexID = params.flexID;
        int[] nTypes = params.nTypes;// composition
        double refFrac = params.refFrac;
        
        double sigmaH = 3.31; 
        double sigmaCH2 = 3.65;//C in CH2
        double sigmaCH3 = 3.30;//C in CH3
        double sigmaCH2H = 0.5 * (sigmaH + sigmaCH2);
        double sigmaCH3H = 0.5 * (sigmaH + sigmaCH3);
        double sigmaCH2CH3 = 0.5 * (sigmaCH2 + sigmaCH3);

        double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;
       
        boolean[] flex = new boolean[]{false, nPoints > 2};//CO2:always rigid;alkane depends on chain length and Bn order

        if (nSpheres==1){
        	throw new RuntimeException("I cannot handle methane!");
        }
        if ( nTypes[0]==0 || nTypes[1]==0){
        	throw new RuntimeException("refer to pure component virial coefficient calculation!");
        }
        if ( (nTypes[0] + nTypes[1])!= nPoints ){
        	throw new RuntimeException("wrong composition!");
        }
        if (!flex[1]){// alkane is rigid
        	if (flexID != -1){
        		throw new RuntimeException("both are rigid, flexID can only be -1!");
        	}
        	System.out.println("both components are rigid, rigid diagram calculation only");
        }
        if (flex[1]){
        	System.out.println("diagrams with letters, alkane is flexible");
        	if (flexID == -1){
        		System.out.println("rigid diagram calculation only");
        	}
        	else if (flexID == 1){
        		System.out.println("flex B diagram calculation only");
        	}
        	else {
        		throw new RuntimeException("wrong flexID!");
        	}
        }
        
        System.out.println("CO2(TraPPE,rigid)+"+nSpheres+"smer n-alkane(TraPPE-EH) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
        temperature = Kelvin.UNIT.toSim(temperature);
        System.out.println("flexID:"+flexID);
        
        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
		
        Space space = Space3D.getInstance();
        // ------------ ref cluster ------------------------- //
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        VirialDiagrams rigidDiagrams = new VirialDiagrams(nPoints, false, false);
        ClusterSum refCluster = rigidDiagrams.makeVirialCluster(fRef);
        
        //-------------- CO2 potential & CO2 mayer function-------------//
        P2CO2TraPPE pCO2 = new P2CO2TraPPE(space);
        MayerGeneral fCO2 = new MayerGeneral(pCO2);
        
        // ------------- Alkane potential & alkane mayer function---------- //
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
        PotentialGroup pAlkaneEH = new PotentialGroup(2);
        MayerGeneral fAlkaneEH = new MayerGeneral(pAlkaneEH);
        
        // ------------ CO2(1st)-alkane(2nd) potential & mayer function ------------//
        P2LennardJones pC_H = new P2LennardJones(space, 0.5*(sigmaH+pCO2.getSigmaC()), Math.sqrt(epsilonH*pCO2.getEpsilonC()));
        P2LennardJones pC_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pCO2.getSigmaC()), Math.sqrt(epsilonCH2*pCO2.getEpsilonC()));
        P2LennardJones pC_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pCO2.getSigmaC()), Math.sqrt(epsilonCH3*pCO2.getEpsilonC()));
        P2LennardJones pO_H = new P2LennardJones(space, 0.5*(sigmaH+pCO2.getSigmaO()), Math.sqrt(epsilonH*pCO2.getEpsilonO()));
        P2LennardJones pO_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pCO2.getSigmaO()), Math.sqrt(epsilonCH2*pCO2.getEpsilonO()));
        P2LennardJones pO_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pCO2.getSigmaO()), Math.sqrt(epsilonCH3*pCO2.getEpsilonO()));
        PotentialGroup pCO2AlkaneEH = new PotentialGroup(2);
        MayerGeneral fCO2AlkaneEH= new MayerGeneral(pCO2AlkaneEH);
        
        // use VirialDiagram to construct target cluster
        VirialDiagramsMix2 diagrams = new VirialDiagramsMix2(nPoints,flex);
        diagrams.setDoReeHoover(false);
        MayerFunction[][] f = new MayerFunction[][]{{fCO2,fCO2AlkaneEH},{fCO2AlkaneEH,fAlkaneEH}};
        ClusterSum targetCluster = diagrams.makeVirialCluster(f,flexID, nTypes);
   
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        // add an array of string
        String[] targetDiagramStrings = new String[0];///////////////////////////
        boolean[] diagramFlexCorrection = new boolean[targetDiagrams.length];
 
        targetDiagrams = diagrams.makeSingleVirialClusters(targetCluster, f, flexID, nTypes);
        targetDiagramNumbers = new int[targetDiagrams.length];
        targetDiagramStrings = new String[targetDiagrams.length];//////////////////////////////

        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = diagrams.getMSMCGraphs(true,flexID,nTypes);
        Map<Graph,Graph> cancelMap = diagrams.getCancelMap();
        int iGraph = 0;
        diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
        	String gString = g.getStore().toNumberString();
        	for (Node node:g.nodes()){
        		gString += node.getColor();
        	}
        	targetDiagramStrings[iGraph] = gString;///////////////////////////////// 
        	System.out.print(iGraph+" ("+g.coefficient()+") "+gString);
        	targetDiagramNumbers[iGraph] = Integer.parseInt(g.getStore().toNumberString());
        	
        	Graph cancelGraph = cancelMap.get(g);
        	if (cancelGraph != null) {
        		diagramFlexCorrection[iGraph] = true;
        		String gnStr = cancelGraph.getStore().toNumberString();// toNumberString: its corresponding diagram number
        		Set<Graph> gSplit = diagrams.getSplitDisconnectedVirialGraphs(cancelGraph);
        		System.out.print(" - "+getSplitGraphString(gSplit, diagrams, false));
        	}
        	System.out.println();
        	iGraph++;
        }
        System.out.println();
        
        Set<Graph> disconnectedGraphs = diagrams.getExtraDisconnectedVirialGraphs(nTypes);
        if ( flexID==-1 && disconnectedGraphs.size() > 0) {
        	System.out.println("shown only when flexID = -1\nextra clusters:");
        	for (Graph g : disconnectedGraphs) {
        		Set<Graph> gSplit = diagrams.getSplitDisconnectedVirialGraphs(g);
        		System.out.print(g.coefficient()+" ");
        		boolean first = true;
        		for (Graph gs : gSplit) {
        			if (!first) {
        				System.out.print(" ");
        			}
        			String gsString = gs.getStore().toNumberString();
        	//		getSplitGraphString(gSplit, diagrams, true);
                	for (Node node:gs.nodes()){
                		gsString += node.getColor();
                	}
    //            	gsString += getSplitGraphString(gSplit, diagrams, true);///?????????????????????????????????????????????????

                	System.out.print(gsString);
        			first = false;
        		}
        		System.out.println();
        		//System.out.println(g);
        	}
        	System.out.println();
        }
        
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        for (int i=0; i<targetDiagrams.length; i++) {
            targetDiagrams[i].setTemperature(temperature);
        }
        
        SpeciesTraPPECO2 speciesCO2 = new SpeciesTraPPECO2(space);// CO2 
        SpeciesAlkaneEH speciesAlkaneEH = new SpeciesAlkaneEH(space, nSpheres);// alkaneEH
        
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), 
        		                                             ClusterWeightAbs.makeWeightCluster(targetCluster)};

        if(flexID != -1){
        	nTypes[flexID]++;// flexID:0=> add a nTypes[0] point; flexID:1=> add a nTypes[1] point
        } 
        
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesCO2,speciesAlkaneEH}, 
        		nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},targetDiagrams,sampleClusters,false);

        if (flexID == 1) {// flex B only
            int[] constraintMap = new int[nPoints+1];
            for (int i=0; i<nPoints; i++) {
                constraintMap[i] = i;
            }
            constraintMap[nPoints] = nTypes[0];// add the duplicate point(of alkane species) to superimpose the 1st alkane point in the diagram
            
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterMoleculeMulti)sim.mcMoveTranslate[1]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[0]).setConstraintMap(constraintMap);
            ((MCMoveClusterRotateMoleculeMulti)sim.mcMoveRotate[1]).setConstraintMap(constraintMap);
            
        }
        
        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;
        
        IAtomType typeCH3 = speciesAlkaneEH.getC_3Type();// C in CH3
        IAtomType typeCH2 = speciesAlkaneEH.getC_2Type();// C in CH2
        IAtomType typeH   = speciesAlkaneEH.getHType();  // H
        IAtomType typeC = speciesCO2.getAtomType(0);
        IAtomType typeO = speciesCO2.getAtomType(1);
        
        // alkane potential
        pAlkaneEH.addPotential(pHH, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeH, typeH}));//H-H
        pAlkaneEH.addPotential(pCH2H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeH, typeCH2}));//H-CH2
        pAlkaneEH.addPotential(pCH3H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeH, typeCH3}));//H-CH3
        pAlkaneEH.addPotential(pCH2H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeH, typeCH2}));//CH2-H
        pAlkaneEH.addPotential(pCH2CH2, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH2}));//CH2-CH2
        pAlkaneEH.addPotential(pCH2CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH3}));//CH2-CH3
        pAlkaneEH.addPotential(pCH3H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeH, typeCH3}));//CH3-H
        pAlkaneEH.addPotential(pCH2CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH2, typeCH3}));//CH3-CH2
        pAlkaneEH.addPotential(pCH3CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeCH3, typeCH3}));//CH3-CH3
        // CO2-alkane potential
        pCO2AlkaneEH.addPotential(pC_CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeC, typeCH3}));
        pCO2AlkaneEH.addPotential(pC_CH2, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeC, typeCH2}));
        pCO2AlkaneEH.addPotential(pC_H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeC, typeH}));
        pCO2AlkaneEH.addPotential(pO_CH3, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeO, typeCH3}));
        pCO2AlkaneEH.addPotential(pO_CH2, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeO, typeCH2}));
        pCO2AlkaneEH.addPotential(pO_H, ApiBuilder.makeIntergroupTypeIterator(new IAtomType[]{typeO, typeH}));
        
        // create the intramolecular potential here, add to it to the "potential master" if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        // *************************** intramolecular potential *************************** //
        //*********************************** torsional potential (XCCH) *********************************//
        double torsionHCCH_c3 = 717; //u(HCCH) = C3 * [1 - cos(3phi)], c3/kB = 717K
    	double torsionCCCH_c3 = 854;//u(CCCH) = C3 * [1 - cos(3phi)], c3/kB = 854K
    	P4BondTorsionAlkaneXCCH p4HCCH = new P4BondTorsionAlkaneXCCH(space, 0.0, 0.0, 0.0, Kelvin.UNIT.toSim(torsionHCCH_c3));
    	P4BondTorsionAlkaneXCCH p4CCCH = new P4BondTorsionAlkaneXCCH(space, 0.0, 0.0, 0.0, Kelvin.UNIT.toSim(torsionCCCH_c3));
    	
        //*********************************** Set Geometric center *********************************//
    	// geometric center is based on all carbons, hydrogens not included
    	AtomPositionGeometricCenterAlkaneEH center = new AtomPositionGeometricCenterAlkaneEH(space,speciesAlkaneEH);
    	((MCMoveRotateMolecule3D)sim.mcMoveRotate[0]).setPositionDefinition(center);
    	((MCMoveRotateMolecule3D)sim.mcMoveRotate[1]).setPositionDefinition(center);
    	((CoordinatePairMoleculeSet)sim.box[0].getCPairSet()).setPositionDefinition(center);
    	((CoordinatePairMoleculeSet)sim.box[1].getCPairSet()).setPositionDefinition(center);
    	
        //*********************************** Add moves *********************************//
        MCMoveClusterRotateCH3[] rotateCH3Move = new MCMoveClusterRotateCH3[2];
        MCMoveClusterWiggleAlkaneEH[] wiggleMove = new MCMoveClusterWiggleAlkaneEH[2];
        MCMoveClusterTorsionAlkaneEH[] torsionMove = null;
            
        if (nSpheres==2){// C2H6, only H-C-C-H torsion potential, 9 sets in total 
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
        sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(1)});///// CALL ONCE IS FINE!,get alkane
        
        /////////////////////////////////// rotateCH3 ////////////////////////////////
        /////////////////////////////////// rotateCH3 ////////////////////////////////           
        System.out.println("----- add rotateCH3 MC move for all types of alkanes-----");
        
        for ( int j = 0; j<2;j++){//reference system and target system
            rotateCH3Move[j] = new MCMoveClusterRotateCH3(sim, sim.integrators[1].getPotentialMaster(), targetCluster.pointCount(), space);
            rotateCH3Move[j].setSpecies(sim.getSpecies(1));// always add to alkane
            sim.integrators[j].getMoveManager().addMCMove(rotateCH3Move[j]);
        }
        	
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
            for ( int j = 0; j<2;j++){//reference system and target system
            	wiggleMove[j] = new MCMoveClusterWiggleAlkaneEH(sim,sim.integrators[1].getPotentialMaster(), targetCluster.pointCount(), space);
            	wiggleMove[j].setSpecies(sim.getSpecies(1));// always add to alkane
            	sim.integrators[j].getMoveManager().addMCMove(wiggleMove[j]);
            }

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
            for ( int j = 0; j<2;j++){//reference system and target system
                torsionMove[j] = new MCMoveClusterTorsionAlkaneEH(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4CCCC, 400);
                torsionMove[j].setTemperature(temperature);
                torsionMove[j].setSpecies(sim.getSpecies(1));// always add to alkane
                sim.integrators[j].getMoveManager().addMCMove(torsionMove[j]);
            }
        }
    	
        // ======== add more vdW intra molecular potential for longer alkane chain ====== //
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
        
        // find a proper configuration (for nominal point and duplicate point)?????????????????
    	double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        for (int j=0; j<10000 && (pi < 1e-10 || Double.isNaN(pi)); j++) {
        	sim.integrators[1].doStep();
        	pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        }
        if ( pi == 0 || Double.isNaN(pi) ) {
            throw new RuntimeException("could not find a configuration for target system:"+pi);
        }
        pi = sim.box[0].getSampleCluster().value(sim.box[0]);
        if ( pi == 0 || Double.isNaN(pi) ) {
            throw new RuntimeException("could not find a configuration for target system:"+pi);
        }
        sim.accumulators[1].reset();// don't want to collect these data!!!!
        
        if (false) {
        	  double size = 10;
              sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
              sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{size,size,size}));
              SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
              DisplayBox dBox0 = simGraphic.getDisplayBox(sim.box[0]);
              DisplayBox dBox1 = simGraphic.getDisplayBox(sim.box[1]);
              dBox0.setPixelUnit(new Pixel(300.0/size));
              dBox1.setPixelUnit(new Pixel(300.0/size));
              dBox0.setShowBoundary(false);
              dBox1.setShowBoundary(false);
              
              //set diameters
              DiameterHashByType diameter = new DiameterHashByType(sim); 
              diameter.setDiameter(speciesCO2.getAtomType(0),0.2);
              diameter.setDiameter(speciesCO2.getAtomType(1),0.3);
              diameter.setDiameter(speciesAlkaneEH.getC_2Type(), 0.3);
              diameter.setDiameter(speciesAlkaneEH.getC_3Type(), 0.4);
              diameter.setDiameter(speciesAlkaneEH.getHType(), 0.2);

              simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
              ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
              colorScheme.setColor(speciesCO2.getAtomType(0), Color.blue);
              colorScheme.setColor(speciesCO2.getAtomType(1), Color.red);
              colorScheme.setColor(speciesAlkaneEH.getC_2Type(), Color.green);
              colorScheme.setColor(speciesAlkaneEH.getC_3Type(), Color.yellow);          
              colorScheme.setColor(speciesAlkaneEH.getHType(), Color.cyan);          

              ((DisplayBoxCanvasG3DSys)dBox1.canvas).setBackgroundColor(Color.WHITE);
              
              simGraphic.makeAndDisplayFrame();

              sim.integratorOS.setNumSubSteps(1000);
              sim.setAccumulatorBlockSize(100);
                  
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
              if (Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0) {
              	throw new RuntimeException("Oops");
              }

              return;
          }
          

        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/100);
        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, steps/40);
        if (sim.refPref == 0 || Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref)) {
            throw new RuntimeException("oops");
        }

        System.out.println("equilibration finished");
        
        sim.setAccumulatorBlockSize(steps);
        sim.integratorOS.setNumSubSteps((int)steps);
        sim.ai.setMaxSteps(1000);
        sim.integratorOS.getMoveManager().setEquilibrating(false);
        
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
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        if (false) {
            IIntegratorListener progressReport = new IIntegratorListener() {
                public void integratorInitialized(IIntegratorEvent e) {}
                public void integratorStepStarted(IIntegratorEvent e) {}
                public void integratorStepFinished(IIntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % sim.ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
        }

        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step frequency "+sim.integratorOS.getRefStepFraction());

        sim.printResults(HSB[nPoints]);
        
        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(sim.accumulators[1].AVERAGE.index);
        IData dataErr = allData.getData(sim.accumulators[1].ERROR.index);
        IData dataCov = allData.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        int nTotal = (targetDiagrams.length+2);
        double oVar = dataCov.getValue(nTotal*nTotal-1);
        for (int i=0; i<targetDiagrams.length; i++) {
            if (targetDiagramNumbers[i]<0) {
                System.out.print("diagram "+(-targetDiagramNumbers[i])+("bc "));
            }
            else {
                System.out.print("diagram "+targetDiagramStrings[i]);//////////////// instead of printing targetDiagramsNumbers, print String
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

    /*
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 2;
        public double temperature = 298;
        public long numSteps = 1000000;
        public int nSpheres = 2;
        public int flexID = -1;
        public int[] nTypes = new int[]{1,1};
        public double refFrac = -1;
    }
}
