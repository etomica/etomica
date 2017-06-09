package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.integrator.IntegratorListener;
import etomica.integrator.IntegratorEvent;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.chem.elements.ElementSimple;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.math.SpecialFunctions;
import etomica.potential.P2CO2TraPPE;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

import java.awt.*;
import java.util.Arrays;

/**
 *   Mayer sampling simulation for CO2(rigid, TraPPE)-CH4(rigid, TraPPE-UA, single LJ site) mixture
 *   cross virial coefficients
 *   Using VirialDiagramMix2 to generate diagrams
 *   both of the components are rigid
 *   Implement Wheatley method
 *   
 *   @author shu
 *   November 2014
 * 
 */
public class VirialCO2CH4UAMixWheatley {
	
	public static void main(String[] args) {
        VirialMixParam params = new VirialMixParam();
        if (args.length > 0) {
			ParseArgs.doParseArgs(params, args);
		} else {
			params.nPoints = 7;
			params.nTypes = new int[]{3,4};
			params.numSteps = 100000000;
			params.temperature = 323.15;
		}
        long t1 = System.currentTimeMillis();
        
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
        System.out.println("wheatley approach");
        System.out.println("\nCO2(TraPPE)+ CH4(TraPPE-UA) overlap sampling B"+nTypes[0]+""+nTypes[1]+" at T="+temperature+" Kelvin");
        System.out.println("both components are rigid, rigid diagram only");
        temperature = Kelvin.UNIT.toSim(temperature);
                
        final double[] HSB = new double[9];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        HSB[8] = Standard.B8HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
       
        double coefficient = SpecialFunctions.factorial(nPoints)/SpecialFunctions.factorial(nTypes[0])/SpecialFunctions.factorial(nTypes[1]) ;
        
//        HSB[nPoints] = HSB[nPoints] * coefficient;
        System.out.println("B"+nPoints+"HS: "+HSB[nPoints]);
 
        Space space = Space3D.getInstance();
        // ------------ ref cluster ------------------------- //
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = new ClusterWheatleyHS(nPoints, fRef);
        refCluster.setTemperature(temperature);
        //-------------- CO2 potential & CO2 mayer function-------------//
        P2CO2TraPPE pCO2 = new P2CO2TraPPE(space);
        MayerGeneral fCO2 = new MayerGeneral(pCO2);        
        // CH4(TraPPE-UA) potential, single site, apply P2LennardJones
        double sigmaCH4 = 3.73;
        double epsilonCH4 = Kelvin.UNIT.toSim(148.0);
        P2LennardJones pCH4 = new P2LennardJones(space, sigmaCH4, epsilonCH4);
        MayerGeneralSpherical fCH4 = new MayerGeneralSpherical(pCH4);
        // ------------ CO2(1st)-CH4(2nd) potential & mayer function ------------//
        P2LennardJones pC_CH4 = new P2LennardJones(space, 0.5*(sigmaCH4 + pCO2.getSigmaC()), Math.sqrt(epsilonCH4 * pCO2.getEpsilonC()));
        P2LennardJones pO_CH4 = new P2LennardJones(space, 0.5*(sigmaCH4 + pCO2.getSigmaO()), Math.sqrt(epsilonCH4 * pCO2.getEpsilonO()));
        PotentialGroup pCO2CH4 = new PotentialGroup(2);
        MayerGeneral fCO2CH4= new MayerGeneral(pCO2CH4); 
        
        // use VirialDiagram to construct target cluster
//        VirialDiagramsMix2 diagrams = new VirialDiagramsMix2(nPoints,flex);
//        diagrams.setDoReeHoover(false);
        MayerFunction[][] f = new MayerFunction[][]{{fCO2,fCO2CH4},{fCO2CH4,fCH4}};
//        ClusterSum targetCluster = diagrams.makeVirialCluster(f,-1, nTypes);//flexID always -1 since both are rigid
//        ClusterSumShell[] targetDiagrams = new ClusterSumShell[0];
//        targetDiagrams = diagrams.makeSingleVirialClusters(targetCluster, f, -1, nTypes);//flexID always -1 since both are rigid
//       
//        
        ClusterAbstract targetCluster = new ClusterWheatleySoftMix(nPoints,nTypes, f, 1e-12);	

        targetCluster.setTemperature(temperature);

        
        SpeciesTraPPECO2 speciesCO2 = new SpeciesTraPPECO2(space);// CO2 
        SpeciesSpheresMono speciesCH4 = new SpeciesSpheresMono(space, new ElementSimple("A"));//CH4
       
        ClusterWeight[] sampleClusters = new ClusterWeight[]{ClusterWeightAbs.makeWeightCluster(refCluster), 
        		                                             ClusterWeightAbs.makeWeightCluster(targetCluster)};

        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space,new ISpecies[]{speciesCO2,speciesCH4}, 
        		nTypes, temperature,refCluster, targetCluster);
        sim.setDoWiggle(false);
//        sim.setRandom(new RandomMersenneTwister(new int[]{-1618344162, 827873331, -207241658, 1776043175}));
        sim.init();
        int[] seeds = sim.getRandomSeeds();
        System.out.println("random seeds: "+Arrays.toString(seeds));

        sim.integratorOS.setNumSubSteps(1000);
        sim.integratorOS.setAggressiveAdjustStepFraction(true);
        System.out.println(steps+" steps (1000 blocks of "+steps/1000+")");
        steps /= 1000;

        AtomType typeCH4 = speciesCH4.getAtomType(0);//C in CH4
        AtomType typeC_CO2 = speciesCO2.getAtomType(0);//  C in CO2
        AtomType typeO_CO2 = speciesCO2.getAtomType(1);// O in CO2
      
        // CO2-CH4 potential
        pCO2CH4.addPotential(pC_CH4, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeC_CO2, typeCH4}));
        pCO2CH4.addPotential(pO_CH4, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeO_CO2, typeCH4}));
                
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
              diameter.setDiameter(speciesCH4.getAtomType(0), 0.3);

              simGraphic.getDisplayBox(sim.box[1]).setDiameterHash(diameter);
              ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box[1]).getColorScheme();
              colorScheme.setColor(speciesCO2.getAtomType(0), Color.blue);
              colorScheme.setColor(speciesCO2.getAtomType(1), Color.red);
              colorScheme.setColor(speciesCH4.getAtomType(0), Color.green);

              ((DisplayBoxCanvasG3DSys)dBox1.canvas).setBackgroundColor(Color.WHITE);
              
              simGraphic.makeAndDisplayFrame();

              sim.integratorOS.setNumSubSteps(1000);
              sim.setAccumulatorBlockSize(100);
                  
              // if running interactively, set filename to null so that it doesn't read
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
		long t2 = System.currentTimeMillis();
		System.out.println("simulation time is:    "+ (t2- t1)*0.001);
		System.out.println("conv is "+Math.pow(6.0221415e23/1e24,nPoints-1));
	}

    /**
     * Inner class for parameters
     */
    public static class VirialMixParam extends ParameterBase {
        public int nPoints = 7;
        public double temperature = 323.15;
        public long numSteps = 1000000000;
        public double sigmaHSRef = 4.5;
        public int[] nTypes = new int[]{4,3};
        public double refFrac = -1;
    }
}
