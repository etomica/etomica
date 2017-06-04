/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.api.IElement;
import etomica.api.IIntegratorEvent;
import etomica.api.IIntegratorListener;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.config.IConformation;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
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
import etomica.potential.*;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;
import etomica.virial.cluster.VirialDiagrams;
import etomica.virial.cluster.VirialDiagramsMix2;

import java.awt.*;
import java.util.Map;
import java.util.Set;

/**
 *   Mayer sampling simulation for N2(rigid, TraPPE)-alkanes(TraPPE-UA) mixture
 *   Martin and Siepmann, JPC(B) 102, 2569-2577 (1998)
 *   eovererr can be used
 *   cross virial coefficients
 *   Using VirialDiagramMix2 to generate diagrams
 *   see "VirialN2AlkaneUARigid" B11 code(use rigid mix diagram)
 *   rigid diagrams and flex diagrams are calculated respectively(Bij = Bij[flexID=-1] + Bij[flexID=1])
 *   
 *   @author shu
 *   April 2013
 * 
 */
public class VirialN2AlkaneUAMix {
	
    public static String getSplitGraphString(Set<Graph> gSplit, VirialDiagramsMix2 flexDiagrams, boolean correction) {
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
                if (VirialDiagramsMix2.graphHasEdgeColor(gs, flexDiagrams.eBond)) {
                    str += "p" + edgeDeleter.apply(gs, ede).getStore().toNumberString();
                }
            }
            if (first && correction) str += "c";
            first = false;
            
           // str +="_";//add a space between number and letters
            //loop over every node in the graph to get the node color,append this char to str
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
        double sigmaHSRef = params.sigmaHSRef;
        double sigmaCH2= params.sigmaCH2;
        double sigmaCH3= params.sigmaCH3;
        int nSpheres = params.nSpheres;// number of carbons
        int flexID = params.flexID;
        int[] nTypes = params.nTypes;// composition
        double refFrac = params.refFrac;
        
        if (nSpheres==1){
        	throw new RuntimeException("I cannot handle methane!");
        }
        if ( nTypes[0]==0 || nTypes[1]==0){
        	throw new RuntimeException("refer to pure component virial coefficient calculation!");
        }
        if ( (nTypes[0]+nTypes[1])!= nPoints ){
        	throw new RuntimeException("wrong composition!");
        }
        System.out.println("diagram with letters");
        System.out.println("N2(TraPPE,rigid)+"+nSpheres+"smer n-alkane(TraPPE-UA) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
        temperature = Kelvin.UNIT.toSim(temperature);
        System.out.println("flexID:"+flexID);
        boolean[] flex = new boolean[]{false,(nSpheres > 2 && nPoints > 2)};//N2:always rigid;alkane depends on chain length and Bn order
        if (!flex[1]){// alkane is rigid
        	if (flexID != -1){
        		throw new RuntimeException("both are rigid, flexID can only be -1!");
        	}
        	System.out.println("both components are rigid, rigid diagram calculation only");
        }
        if (flex[1]){
        	System.out.println("alkane is flexible");
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
        // ------------ potential parameters ------------------------- //
        double sigmaA = 1.0;// set arbitrarily
        double sigmaN = 3.31;
        double epsilonA = 0.0;
        double epsilonN = 36.0;
        double charge = 0.964;// in the center of mass
        double sigmaAN= (sigmaA + sigmaN) * 0.5 ;
        double epsilonAN = 0.0;
        //-------------- N2 potential & N2 mayer function-------------//
        P2CO2EMP pN2 = new P2CO2EMP(space, sigmaA, sigmaAN, sigmaN,epsilonA ,epsilonAN, Kelvin.UNIT.toSim(epsilonN), Electron.UNIT.toSim(charge));
        MayerGeneral fN2 = new MayerGeneral(pN2);
        // ------------- Alkane potential & alkane mayer function---------- //
        double epsilonCH2 = Kelvin.UNIT.toSim(46.0);
        double epsilonCH3 = Kelvin.UNIT.toSim(98.0);
        double epsilonCH2CH3 = Math.sqrt(epsilonCH2*epsilonCH3);
        P2LennardJones p2CH2 = new P2LennardJones(space, sigmaCH2, epsilonCH2);
        P2LennardJones p2CH3 = new P2LennardJones(space, sigmaCH3, epsilonCH3);
        P2LennardJones p2CH2CH3 = new P2LennardJones(space, 0.5*(sigmaCH2+sigmaCH3), epsilonCH2CH3);
        PotentialGroup pAlkane = new PotentialGroup(2);
        MayerGeneral fAlkane= new MayerGeneral(pAlkane);
        // ------------ N2(1st)-alkane(2nd) potential & mayer function ------------//
        P2LennardJones pN_CH2 = new P2LennardJones(space, 0.5*(sigmaCH2+pN2.getSigmaO()), Math.sqrt(epsilonCH2*pN2.getEpsilonO()));
        P2LennardJones pN_CH3 = new P2LennardJones(space, 0.5*(sigmaCH3+pN2.getSigmaO()), Math.sqrt(epsilonCH3*pN2.getEpsilonO()));
        PotentialGroup pN2Alkane = new PotentialGroup(2);
        MayerGeneral fN2Alkane= new MayerGeneral(pN2Alkane);
        
        // use VirialDiagram to construct target cluster
        VirialDiagramsMix2 diagrams = new VirialDiagramsMix2(nPoints,flex);
        diagrams.setDoReeHoover(false);
        MayerFunction[][] f = new MayerFunction[][]{{fN2,fN2Alkane},{fN2Alkane,fAlkane}};
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
        		Set<Graph> gSplit = diagrams.getSplitDisconnectedVirialGraphs(g);//use happy ap ???????
        		System.out.print(g.coefficient()+" ");
        		boolean first = true;
        		for (Graph gs : gSplit) {
        			if (!first) {
        				System.out.print(" ");
        			}
        			String gsString = gs.getStore().toNumberString();
                	for (Node node:gs.nodes()){
                		gsString += node.getColor();
                	}
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
        // nitrogen
        final IConformation conformation = new IConformation() {
        	public void initializePositions(IAtomList atomList) {
                // atoms are C, O and O, so we arrange them as 1-0-2
                double bondL = 1.10;
                atomList.getAtom(0).getPosition().E(0);
                atomList.getAtom(1).getPosition().E(0);
                atomList.getAtom(1).getPosition().setX(0, -0.5 * bondL);
                atomList.getAtom(2).getPosition().E(0);
                atomList.getAtom(2).getPosition().setX(0, + 0.5 * bondL);
            }
        };
        SpeciesSpheresHetero speciesN2 = new SpeciesSpheresHetero(space, new IElement[]{new ElementSimple("A"), Nitrogen.INSTANCE});
        speciesN2.setChildCount(new int[]{1,2});
        speciesN2.setConformation(conformation);
        // alkane
        SpeciesAlkane speciesAlkane = new SpeciesAlkane(space, nSpheres);
        
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), 
        		                                             ClusterWeightAbs.makeWeightCluster(targetCluster)};

        if(flexID != -1){
        	nTypes[flexID]++;// flexID:0=> add a nTypes[0] point; flexID:1=> add a nTypes[1] point
        } 
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesN2,speciesAlkane}, 
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
        
        MCMoveClusterWiggleMulti[] wiggleMove = new MCMoveClusterWiggleMulti[2];
        if (nSpheres==2){// if ethane, then do not do wiggle
        	System.out.println("n=2, no wiggles !");
        }
        else {
        	System.out.println("not ethane, add wiggle move!");
        	for (int i=0; i<2; i++) {
        		wiggleMove[i] = new MCMoveClusterWiggleMulti(sim, sim.integrators[i].getPotentialMaster(), targetCluster.pointCount(), space);
        		wiggleMove[i].setSpecies(sim.getSpecies(1));
        		sim.integrators[i].getMoveManager().addMCMove(wiggleMove[i]);
        	}
        }
        
        sim.integratorOS.setNumSubSteps(1000);

        AtomType typeCH3 = speciesAlkane.getAtomType(0);
        AtomType typeCH2 = speciesAlkane.getAtomType(1);
        AtomType typeA = speciesN2.getAtomType(0);
        AtomType typeN = speciesN2.getAtomType(1);
        
        // alkane potential
        pAlkane.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH2}));
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2, typeCH3}));// CH2 on molecule1 to CH3 on molecule2
        pAlkane.addPotential(p2CH2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH2}));// CH3 on molecule1 to CH2 on molecule2
        pAlkane.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3, typeCH3}));
        // N2-alkane potential
        pN2Alkane.addPotential(pN_CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeCH3}));
        pN2Alkane.addPotential(pN_CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeN, typeCH2}));
        
        sim.integratorOS.setNumSubSteps(1000);

        // create the intramolecular potential here, add to it to the potential master if needed
        PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        //*********************************************************************************************************//
        // *************************** intramolecular potential *************************** //
        //*********************************************************************************************************//
        if (nSpheres > 2) {
        	System.out.println("n>2, add bending potential to intramolecular potential!");
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
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(1)});
        }
        MCMoveClusterTorsionMulti[] torsionMoves = null;
        if (nSpheres > 3) {
        	System.out.println("n>3, add torsion potential to intramolecular potential!");
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres-3][4];
            for (int i=0; i<nSpheres-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            
            pIntra.addPotential(p4, new Atomset4IteratorIndexList(quads));
            
            System.out.println("n>3, add torsion mc move!");
            torsionMoves = new MCMoveClusterTorsionMulti[2];
            for (int j = 0; j<2; j++){   
            	torsionMoves[j] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            	torsionMoves[j].setTemperature(temperature);
                torsionMoves[j].setSpecies(sim.getSpecies(1));
                sim.integrators[j].getMoveManager().addMCMove(torsionMoves[j]);

            }
            
        }
 
        // =====================add move 2 and change 1&3============================
//        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
//        IVectorRandom groupTranslationVector = (IVectorRandom)translator.getTranslationVector();
//        MoleculeChildAtomAction moveMoleculeAction = new MoleculeChildAtomAction(translator);
//        IMoleculeList molecules = sim.box[1].getMoleculeList();
//        double r = 4;
//        // put the molecules in a ring around the origin, with one atom
//        // from each scaled in toward the origin
//        int totalPoints = 0;
//        for (int t=0; t<nTypes.length;t++){
//        	totalPoints +=nTypes[t];
//        }
//    	for (int i=1; i<totalPoints; i++) {
//        	groupTranslationVector.setX(0, r*Math.cos(2*(i-1)*Math.PI/(totalPoints-1)));
//            groupTranslationVector.setX(1, r*Math.sin(2*(i-1)*Math.PI/(totalPoints-1)));
//            moveMoleculeAction.actionPerformed(molecules.getMolecule(i));
//        }
    	
    	// ======== add more vdW intra molecular potential for longer alkane chain ====== //
        if (nSpheres > 4) {
        	System.out.println("n>4, add CH3-CH3 pair intramolecular potential!");
            pIntra.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres-1}}));
        }
        if (nSpheres > 5) {
        	System.out.println("n>5, add CH3-CH2 pair intramolecular potential!");
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
        	System.out.println("n>6, add CH2-CH2 pair intramolecular potential!");
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
        
        // find a proper configuration (for nominal point and duplicate point)?????????????????
        double pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        for (int j=0; j<10000 && (pi < 1e-10 || Double.isNaN(pi)); j++) {
            sim.integrators[1].doStep();
            pi = sim.box[1].getSampleCluster().value(sim.box[1]);
        }
        if ( pi == 0) {
            throw new RuntimeException("could not find a configuration for target system");
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
              diameter.setDiameter(speciesN2.getAtomType(0),0.2);
              diameter.setDiameter(speciesN2.getAtomType(1),0.3);
              diameter.setDiameter(speciesAlkane.getCH2Type(), 0.3);
              diameter.setDiameter(speciesAlkane.getCH3Type(), 0.4);

              simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
              ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
              colorScheme.setColor(speciesN2.getAtomType(0), Color.blue);
              colorScheme.setColor(speciesN2.getAtomType(1), Color.red);
              colorScheme.setColor(speciesAlkane.getCH2Type(), Color.green);
              colorScheme.setColor(speciesAlkane.getCH3Type(), Color.yellow);          

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
        
        for (int i=0; i<2; i++) {
            System.out.println("MC Move step sizes "+sim.mcMoveTranslate[i].getStepSize()+" "+sim.mcMoveRotate[i].getStepSize());
        }
        if (nSpheres > 2) {
            System.out.println("Wiggle move acceptance "+ wiggleMove[0].getTracker().acceptanceRatio()+" "+
     	       		wiggleMove[1].getTracker().acceptanceRatio());
        	System.out.println("Wiggle move step sizes " + wiggleMove[0].getStepSize() + " "+	wiggleMove[1].getStepSize());
        }
        
        if (nSpheres > 3) {
        	System.out.println("Torsion move acceptance "+torsionMoves[0].getTracker().acceptanceRatio()+" "+
    	                torsionMoves[1].getTracker().acceptanceRatio());
        	System.out.println("Torsion move step size "+torsionMoves[0].getStepSize() + " "+   torsionMoves[1].getStepSize());
           
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

    /**
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 4;
        public double temperature = 320;
        public long numSteps = 1000000;
        public int nSpheres = 20;
        public double sigmaCH2 = 3.95;// in alkane
        public double sigmaCH3 = 3.75;// in alkane
        public double sigmaHSRef = sigmaCH3 + 0.5*nSpheres;
        public int flexID = -1;
        public int[] nTypes = new int[]{2,2};
        public double refFrac = -1;
    }
}
