/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ApiBuilder;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Nitrogen;
import etomica.config.IConformation;
import etomica.graph.model.Graph;
import etomica.graph.model.Node;
import etomica.graph.operations.DeleteEdge;
import etomica.graph.operations.DeleteEdgeParameters;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.models.co2.SpeciesTraPPECO2;
import etomica.potential.P2CO2EMP;
import etomica.potential.P2CO2TraPPE;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
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
 *   Mayer sampling simulation for CO2(rigid, TraPPE)-N2(rigid, TraPPE) mixture
 *   cross virial coefficients
 *   Using VirialDiagramMix2 to generate diagrams
 *   rigid diagrams only since both of the components are rigid
 *   
 *   @author shu
 *   May 2013
 * 
 */
public class VirialCO2N2Mix {
	
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
        int[] nTypes = params.nTypes;// composition
        double refFrac = params.refFrac;
        
        if ( nTypes[0]==0 || nTypes[1]==0){
        	throw new RuntimeException("refer to pure component virial coefficient calculation!");
        }
        if ( (nTypes[0]+nTypes[1])!= nPoints ){
        	throw new RuntimeException("wrong composition!");
        }
        System.out.println("\nCO2(TraPPE)+ N2(TraPPE) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
        System.out.println("both components are rigid, rigid diagram only");
        temperature = Kelvin.UNIT.toSim(temperature);
        boolean[] flex = new boolean[]{false,false};// both are rigid
                
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
        // ------------ N2 potential parameters ------------------------- //
        double sigmaA = 1.0;// set arbitrarily
        double sigmaN = 3.31;
        double epsilonN = 36.0;
        double charge = 0.964;// at the center of mass
        double sigmaAN= (sigmaA + sigmaN) * 0.5 ;
        //-------------- N2 potential & N2 mayer function-------------//
        P2CO2EMP pN2 = new P2CO2EMP(space, sigmaA, sigmaAN, sigmaN, 0.0, 0.0, Kelvin.UNIT.toSim(epsilonN), Electron.UNIT.toSim(charge));
        MayerGeneral fN2 = new MayerGeneral(pN2);
        // ------------ CO2(1st)-N2(2nd) potential & mayer function ------------//
        P2LennardJones pC_N = new P2LennardJones(space, 0.5*(sigmaN + pCO2.getSigmaC()), Math.sqrt(epsilonN * pCO2.getEpsilonC()));
        P2LennardJones pO_N = new P2LennardJones(space, 0.5*(sigmaN + pCO2.getSigmaO()), Math.sqrt(epsilonN * pCO2.getEpsilonO()));
        PotentialGroup pCO2N2 = new PotentialGroup(2);
        MayerGeneral fCO2N2= new MayerGeneral(pCO2N2); 
        
        // use VirialDiagram to construct target cluster
        VirialDiagramsMix2 diagrams = new VirialDiagramsMix2(nPoints,flex);
        diagrams.setDoReeHoover(false);
        MayerFunction[][] f = new MayerFunction[][]{{fCO2,fCO2N2},{fCO2N2,fN2}};
        ClusterSum targetCluster = diagrams.makeVirialCluster(f,-1, nTypes);//flexID always -1 since both are rigid
        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
        int[] targetDiagramNumbers = new int[0];
        boolean[] diagramFlexCorrection = new boolean[targetDiagrams.length];
        targetDiagrams = diagrams.makeSingleVirialClusters(targetCluster, f, -1, nTypes);//flexID always -1 since both are rigid
        targetDiagramNumbers = new int[targetDiagrams.length];

        System.out.println("individual clusters:");
        Set<Graph> singleGraphs = diagrams.getMSMCGraphs(true,-1,nTypes);//flexID always -1 since both are rigid
        Map<Graph,Graph> cancelMap = diagrams.getCancelMap();
        int iGraph = 0;
        diagramFlexCorrection = new boolean[targetDiagrams.length];
        for (Graph g : singleGraphs) {
        	String gString = g.getStore().toNumberString();
        	for (Node node:g.nodes()){
        		gString += node.getColor();
        	}
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
        if ( disconnectedGraphs.size() > 0) {
        	System.out.println("shown only when flexID = -1\nextra clusters:");
                
        	for (Graph g : disconnectedGraphs) {
        		Set<Graph> gSplit = diagrams.getSplitDisconnectedVirialGraphs(g);
        		System.out.print(g.coefficient()+" ");
        		boolean first = true;
        		for (Graph gs : gSplit) {
        			if (!first) {
        				System.out.print("*");
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

        SpeciesTraPPECO2 speciesCO2 = new SpeciesTraPPECO2(space);// CO2
        // N2
        final IConformation conformation = new IConformation() {
        	public void initializePositions(IAtomList atomList) {
                // atoms are C, O and O, so we arrange them as 1-0-2
                double bondL = 1.10;
                atomList.get(0).getPosition().E(0);
                atomList.get(1).getPosition().E(0);
                atomList.get(1).getPosition().setX(0, -0.5 * bondL);
                atomList.get(2).getPosition().E(0);
                atomList.get(2).getPosition().setX(0, + 0.5 * bondL);
            }
        };
        SpeciesSpheresHetero speciesN2 = new SpeciesSpheresHetero(space, new IElement[]{new ElementSimple("A"), Nitrogen.INSTANCE});
        speciesN2.setChildCount(new int[]{1,2});
        speciesN2.setConformation(conformation);        
        
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), 
        		                                             ClusterWeightAbs.makeWeightCluster(targetCluster)};

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesCO2,speciesN2}, 
        		nTypes, temperature,new ClusterAbstract[]{refCluster,targetCluster},targetDiagrams,sampleClusters,false);

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        AtomType typeC = speciesCO2.getAtomType(0);//  C in CO2
        AtomType typeO = speciesCO2.getAtomType(1);// O in CO2
        AtomType typeN = speciesN2.getAtomType(1);//N in N2

   
        // CO2-N2 potential
        pCO2N2.addPotential(pC_N, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC, typeN}));
        pCO2N2.addPotential(pO_N, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO, typeN}));
                
        sim.integratorOS.setNumSubSteps(1000);

        // find a proper configuration 
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
              sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
              sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{size, size, size}));
              SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
              DisplayBox dBox0 = simGraphic.getDisplayBox(sim.box[0]);
              DisplayBox dBox1 = simGraphic.getDisplayBox(sim.box[1]);
              dBox0.setPixelUnit(new Pixel(300.0/size));
              dBox1.setPixelUnit(new Pixel(300.0/size));
              dBox0.setShowBoundary(false);
              dBox1.setShowBoundary(false);
              
              //set diameters
              DiameterHashByType diameter = new DiameterHashByType();
              diameter.setDiameter(speciesCO2.getAtomType(0),0.2);
              diameter.setDiameter(speciesCO2.getAtomType(1),0.3);
              diameter.setDiameter(speciesN2.getAtomType(0), 0.3);
              diameter.setDiameter(speciesN2.getAtomType(1), 0.4);

              simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
              ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
              colorScheme.setColor(speciesCO2.getAtomType(0), Color.blue);
              colorScheme.setColor(speciesCO2.getAtomType(1), Color.red);
              colorScheme.setColor(speciesN2.getAtomType(0), Color.green);
              colorScheme.setColor(speciesN2.getAtomType(1), Color.yellow);          

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
                   
        if (refFrac >= 0) {
            sim.integratorOS.setRefStepFraction(refFrac);
            sim.integratorOS.setAdjustStepFraction(false);
        }

        if (false) {
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
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
        
        
	}

    /**
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 5;
        public double temperature = 298;
        public long numSteps = 1000000;

        public double sigmaHSRef = 4.5;
        public int[] nTypes = new int[]{1,4};
        public double refFrac = -1;
    }
}
