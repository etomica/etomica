/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.potential.P2Exp6Buckingham;
import etomica.potential.P3BondAngle;
import etomica.potential.P4BondTorsion;
import etomica.potential.PotentialGroup;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.util.ParameterBase;
import etomica.virial.*;
import etomica.virial.cluster.Standard;

public class VirialAlkaneMix2_Kong {


    public static void main(String[] args) {
        VirialAlkaneMixParam params = new VirialAlkaneMixParam();
        final int nPoints;
        int nSpheres1;
        int nSpheres2;
        double temperature;
        long steps;
        int nComp1;
        int nComp2;
        
        
        if (args.length == 0) {
        	nPoints = params.nPoints;
            nSpheres1 = params.nSpheres1;//the number of unit of C of component1 
            nSpheres2 = params.nSpheres2;//the number of unit of C of component2
            temperature = params.temperature;
            steps = params.numSteps;
            nComp1 = params.nComp1;//the number of component1
            nComp2 = params.nComp2;//the number of component2
        	
        } else if (args.length == 7) {
        	nPoints = Integer.parseInt(args[0]);
            nSpheres1 = Integer.parseInt(args[1]);
            nSpheres2 = Integer.parseInt(args[2]);
            temperature = Double.parseDouble(args[3]);
            steps = Integer.parseInt(args[4]);
            nComp1 = Integer.parseInt(args[5]);//the number of component1
            nComp2 = Integer.parseInt(args[6]);//the number of component2
        } else {
        	throw new IllegalArgumentException("Wrong number of arguments");
        }
        
        double sigmaCH3 = 3.679;
        double sigmaCH2 = 4.000;
        double bondL1 = 1.839;//bond length of CH3-CH3
        double bondL2 = 1.687;//bond length of CH3-CH2
        double bondL3 = 1.535;//bond length of CH2-CH2
        double sigmaHSRef = 0.5*((sigmaCH3+0.5*nSpheres1)+(sigmaCH3+0.5*nSpheres2));//we need to modify by compare relative uncertainty
        final double[] HSB = new double[8];
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        System.out.println("B2HS: "+HSB[2]);
        System.out.println("B3HS: "+HSB[3]+" = "+(HSB[3]/(HSB[2]*HSB[2]))+" B2HS^2");
        System.out.println("B4HS: "+HSB[4]+" = "+(HSB[4]/(HSB[2]*HSB[2]*HSB[2]))+" B2HS^3");
        System.out.println("B5HS: "+HSB[5]+" = 0.110252 B2HS^4");
        System.out.println("B6HS: "+HSB[6]+" = 0.03881 B2HS^5");
        System.out.println("B7HS: "+HSB[7]+" = 0.013046 B2HS^6");
       
		
        Space space = Space3D.getInstance();
        
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);//Mayer f-function of HS potential
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);//Mayer e-function of HS potential
        //PotentialGroup pMethaneMethaneGroup = new PotentialGroup(2, space);//CH4-CH4 potential group
        //PotentialGroup pMethaneMethylGroup = new PotentialGroup(2, space);//CH4-CH3 potential group
        PotentialGroup pComp1Comp1Group = new PotentialGroup(2);//comp1-comp1 potential group
        PotentialGroup pComp1Comp2Group = new PotentialGroup(2);//comp1-comp2 potential group
        PotentialGroup pComp2Comp2Group = new PotentialGroup(2);//comp2-comp2 potential group
        
        System.out.println("Mixture:"+"C"+nSpheres1+ "+ C"+nSpheres2+"  B"+nPoints+" at "+temperature+"K");
        temperature = Kelvin.UNIT.toSim(temperature);
        //double epsilonCH4 = Kelvin.UNIT.toSim(160.3);// the minimum potential energy
        double epsilonCH3 = Kelvin.UNIT.toSim(129.6);
        double epsilonCH2 = Kelvin.UNIT.toSim(73.5);
        //double epsilonCH4CH3 = Math.sqrt(epsilonCH4*epsilonCH3);
        double epsilonCH3CH2 = 95.25797370;
        //double alphaCH4 = 15;// repulsive steepness of the potential
        double alphaCH3 = 16;
        double alphaCH2 = 22;
        //double alphaCH4CH3 = Math.sqrt(alphaCH4*alphaCH3);
        double alphaCH3CH2 = 18.73414913;
        //double rmCH4 = 4.183766451;//the distance at which the potential is minimum
        double rmCH3 = 4.094113682;
        double rmCH2 = 4.359949403;
        //double rmCH4CH3 = 0.5*(rmCH4 + rmCH3);
        double rmCH3CH2 = 4.245468200;
        //double rmaxCH4 = 0.704;//the potential reaches an unphysical maximum
        double rmaxCH3 = 0.575;
        double rmaxCH2 = 0.22;
        //double rmaxCH4CH3 = 0.5*(rmaxCH4 + rmaxCH3);
        double rmaxCH3CH2 = 0.368577;
        //P2Exp6Buckingham p2CH4 = new P2Exp6Buckingham(space, epsilonCH4, alphaCH4, rmCH4, rmaxCH4);
       /* for (int i=5; i<100;i++){
        	System.out.println(i/10.0+" "+p2CH4.u(i*i/100.0));
            System.out.println(rmCH4+" "+p2CH4.u(rmCH4*rmCH4));
          }
          System.exit(1);*/
        P2Exp6Buckingham p2CH3 = new P2Exp6Buckingham(space, epsilonCH3, alphaCH3, rmCH3, rmaxCH3);
        P2Exp6Buckingham p2CH2 = new P2Exp6Buckingham(space, epsilonCH2, alphaCH2, rmCH2, rmaxCH2);
        //P2Exp6Buckingham p2CH4CH3 = new P2Exp6Buckingham(space, epsilonCH4CH3, alphaCH4CH3, rmCH4CH3, rmaxCH4CH3);
        P2Exp6Buckingham p2CH3CH2 = new P2Exp6Buckingham(space, epsilonCH3CH2, alphaCH3CH2, rmCH3CH2, rmaxCH3CH2);
        
        //MayerGeneral fMethaneMethaneTarget = new MayerGeneral(pMethaneMethaneGroup);//Mayer function of CH4-CH4
        //MayerGeneral fMethaneMethylTarget = new MayerGeneral(pMethaneMethylGroup);
        MayerGeneral fComp1Comp1Target = new MayerGeneral(pComp1Comp1Group);
        MayerGeneral fComp1Comp2Target = new MayerGeneral(pComp1Comp2Group);
        MayerGeneral fComp2Comp2Target = new MayerGeneral(pComp2Comp2Group);
        //MayerEGeneral eMethaneMethaneTarget = new MayerEGeneral(pMethaneMethaneGroup);
        //MayerEGeneral eMethaneMethylTarget = new MayerEGeneral(pMethaneMethylGroup);
        MayerEGeneral eComp1Comp1Target = new MayerEGeneral(pComp1Comp1Group);
        MayerEGeneral eComp1Comp2Target = new MayerEGeneral(pComp1Comp2Group);
        MayerEGeneral eComp2Comp2Target = new MayerEGeneral(pComp2Comp2Group);
        int[] nTypes = new int[]{nComp1,nComp2};      
        ClusterAbstract targetCluster = Standard.virialClusterMixture(nPoints, new MayerFunction[][]{{fComp1Comp1Target,fComp1Comp2Target},{fComp1Comp2Target,fComp2Comp2Target}},
                new MayerFunction[][]{{eComp1Comp1Target,eComp1Comp2Target},{eComp1Comp2Target,eComp2Comp2Target}}, nTypes);
        targetCluster.setTemperature(temperature);
        
        ClusterAbstract refCluster = Standard.virialCluster(nPoints, fRef, nPoints>3, eRef, true);
        refCluster.setTemperature(temperature);

        System.out.println((steps*1000)+" steps ("+steps+" blocks of 1000)");
        
        ElementSimple CH3element = new ElementSimple("CH3",15);
        ElementSimple CH2element = new ElementSimple("CH2",14);
        SpeciesFactorySpheres2 speciesFactoryComp1 = new SpeciesFactorySpheres2(space,nSpheres1, CH3element, CH2element);//make alkane model-component1
        if (nSpheres1 == 2) {
        	speciesFactoryComp1.setBondL(bondL1);
        }
        if (nSpheres1 ==3) {
        	speciesFactoryComp1.setBondL(bondL2);
        } else {
        	speciesFactoryComp1.setBondL(bondL3);
        }
        SpeciesFactorySpheres2 speciesFactoryComp2 = new SpeciesFactorySpheres2(space,nSpheres2,CH3element, CH2element);//make alkane model-component2
        if (nSpheres2 == 2) {
        	speciesFactoryComp2.setBondL(bondL1);
        }
        if (nSpheres2 ==3) {
        	speciesFactoryComp2.setBondL(bondL2);
        } else {
        	speciesFactoryComp1.setBondL(bondL3);
        }

        SpeciesFactory[] speciesFactory = new SpeciesFactory[2];
        
        speciesFactory[0] = speciesFactoryComp1;
        speciesFactory[1] = speciesFactoryComp2;
        
        
        final SimulationVirialMultiOverlap sim = new SimulationVirialMultiOverlap(space, speciesFactory,
                          temperature,refCluster,targetCluster, nSpheres1 > 2 || nSpheres2 > 2, new int[]{nComp1,nComp2} );//overlap-sampling approach to evaluating a cluster diagram.
        //        sim.integratorOS.setAdjustStepFreq(false);
//        sim.integratorOS.setStepFreq0(1);


        SpeciesAlkane species1 = (SpeciesAlkane)sim.species[0];
        SpeciesAlkane species2 = (SpeciesAlkane)sim.species[1];
        AtomType typeCH3A = species1.getAtomType(0);
        AtomType typeCH3B = species2.getAtomType(0);
        AtomType typeCH2A = species1.getAtomType(1);
        AtomType typeCH2B = species2.getAtomType(1);
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3A, typeCH2A}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3A, typeCH2B}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3B, typeCH2A}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3B, typeCH2B}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2A, typeCH3A}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2A, typeCH3B}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2B, typeCH3A}));
        pComp1Comp2Group.addPotential(p2CH3CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2B, typeCH3B}));
        pComp1Comp1Group.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3A, typeCH3A}));
        pComp1Comp1Group.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3A, typeCH3B}));
        pComp1Comp1Group.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3B, typeCH3A}));
        pComp1Comp1Group.addPotential(p2CH3, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH3B, typeCH3B}));
        pComp2Comp2Group.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2A, typeCH2A}));
        pComp2Comp2Group.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2A, typeCH2B}));
        pComp2Comp2Group.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2B, typeCH2A}));
        pComp2Comp2Group.addPotential(p2CH2, ApiBuilder.makeIntergroupTypeIterator(new AtomType[]{typeCH2B, typeCH2B}));
        
        sim.integratorOS.setNumSubSteps(1000);
        PotentialGroup pIntra1 = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        PotentialGroup pIntra2 = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        
        
        if (nSpheres1 > 2) {
            P3BondAngle p3 = new P3BondAngle(space);
            p3.setAngle(Math.PI*114.0/180.0);
            p3.setEpsilon(Kelvin.UNIT.toSim(62500));
            int[][] triplets = new int[nSpheres1-2][3];
            for (int i=0; i<nSpheres1-2; i++) {
                triplets[i][0] = i;
                triplets[i][1] = i+1;
                triplets[i][2] = i+2;
            }
            pIntra1.addPotential(p3, new Atomset3IteratorIndexList(triplets)); 
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra1,new ISpecies[]{sim.species[0]});
        }
        MCMoveClusterTorsionMulti[] torsionMoves1 = null;
        if (nSpheres1 > 3) {
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres1-3][4];
            for (int i=0; i<nSpheres1-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            pIntra1.addPotential(p4, new Atomset4IteratorIndexList(quads));
            torsionMoves1 = new MCMoveClusterTorsionMulti[2];
            torsionMoves1[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves1[0].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMoves1[0]);
            torsionMoves1[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves1[1].setTemperature(temperature);
            sim.integrators[1].getMoveManager().addMCMove(torsionMoves1[1]);
        }
        if (nSpheres1 > 4) {
            pIntra1.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres1-1}}));
        }
        if (nSpheres1 > 5) {
            int[][] pairs = new int[2*(nSpheres1-5)][2];
            for (int i=0; i<nSpheres1-5; i++) {
                pairs[2*i][0] = 0;
                pairs[2*i][1] = nSpheres1-2-i;
                pairs[2*i+1][0] = nSpheres1-1;
                pairs[2*i+1][1] = i+1;
            }
            pIntra1.addPotential(p2CH3CH2,new ApiIndexList(pairs));
        }
        if (nSpheres1 > 6) {
            int[][] pairs = new int[(nSpheres1-6)*(nSpheres1-5)/2][2];
            int k = 0;
            for (int i=1; i<nSpheres1-5; i++) {
                for (int j=i+4; j<nSpheres1-1; j++) {
                    pairs[k][0] = i;
                    pairs[k][1] = j;
                    k++;
                }
            }
            pIntra1.addPotential(p2CH2,new ApiIndexList(pairs));
        }
        
     
        if (nSpheres2 > 2) {
            P3BondAngle p3 = new P3BondAngle(space);
            p3.setAngle(Math.PI*114.0/180.0);
            p3.setEpsilon(Kelvin.UNIT.toSim(62500));
            int[][] triplets = new int[nSpheres2-2][3];
            for (int i=0; i<nSpheres2-2; i++) {
                triplets[i][0] = i;
                triplets[i][1] = i+1;
                triplets[i][2] = i+2;
            }
            pIntra2.addPotential(p3, new Atomset3IteratorIndexList(triplets));
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra2,new ISpecies[]{sim.species[1]});
        }
        MCMoveClusterTorsionMulti[] torsionMoves2 = null;
        if (nSpheres2 > 3) {
            P4BondTorsion p4 = new P4BondTorsion(space, 0, Kelvin.UNIT.toSim(355.03), Kelvin.UNIT.toSim(-68.19), Kelvin.UNIT.toSim(791.32));
            int[][] quads = new int[nSpheres2-3][4];
            for (int i=0; i<nSpheres2-3; i++) {
                quads[i][0] = i;
                quads[i][1] = i+1;
                quads[i][2] = i+2;
                quads[i][3] = i+3;
            }
            pIntra2.addPotential(p4, new Atomset4IteratorIndexList(quads));
            torsionMoves2 = new MCMoveClusterTorsionMulti[2];
            torsionMoves2[0] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves2[0].setTemperature(temperature);
            sim.integrators[0].getMoveManager().addMCMove(torsionMoves2[0]);
            torsionMoves2[1] = new MCMoveClusterTorsionMulti(sim.integrators[1].getPotentialMaster(), space, sim.getRandom(), 1.0, p4, 40);
            torsionMoves2[1].setTemperature(temperature);
            sim.integrators[1].getMoveManager().addMCMove(torsionMoves2[1]);
        }
        if (nSpheres2 > 4) {
            pIntra2.addPotential(p2CH3,new ApiIndexList(new int[][]{{0,nSpheres2-1}}));
        }
        if (nSpheres2 > 5) {
            int[][] pairs = new int[2*(nSpheres2-5)][2];
            for (int i=0; i<nSpheres2-5; i++) {
                pairs[2*i][0] = 0;
                pairs[2*i][1] = nSpheres2-2-i;
                pairs[2*i+1][0] = nSpheres2-1;
                pairs[2*i+1][1] = i+1;
            }
            pIntra2.addPotential(p2CH3CH2,new ApiIndexList(pairs));
        }
        if (nSpheres2 > 6) {
            int[][] pairs = new int[(nSpheres2-6)*(nSpheres2-5)/2][2];
            int k = 0;
            for (int i=1; i<nSpheres2-5; i++) {
                for (int j=i+4; j<nSpheres2-1; j++) {
                    pairs[k][0] = i;
                    pairs[k][1] = j;
                    k++;
                }
            }
            pIntra2.addPotential(p2CH2,new ApiIndexList(pairs));
        }
                               
        if (false) {
            sim.box[0].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            sim.box[1].getBoundary().setBoxSize(Vector.of(new double[]{10, 10, 10}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            displayBox0.setShowBoundary(false);
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
            //diameterManager.setDiameter(typeCH4, sigmaCH4);
            diameterManager.setDiameter(typeCH3A, sigmaCH3);
            diameterManager.setDiameter(typeCH3B, sigmaCH3);
            diameterManager.setDiameter(typeCH2A, sigmaCH2);
            diameterManager.setDiameter(typeCH2B, sigmaCH2);
            displayBox1.setDiameterHash(diameterManager);
            simGraphic.makeAndDisplayFrame();

            sim.integratorOS.setNumSubSteps(1000);
            sim.setAccumulatorBlockSize(1000);
                
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            sim.getController().removeAction(sim.ai);
            sim.getController().addAction(new IAction() {
                public void actionPerformed() {
                    sim.initRefPref(null, 100);
                    sim.equilibrate(null, 200);
                    sim.ai.setMaxSteps(Long.MAX_VALUE);
                }
            });
            sim.getController().addAction(sim.ai);
            if ((Double.isNaN(sim.refPref) || Double.isInfinite(sim.refPref) || sim.refPref == 0)) {
                throw new RuntimeException("Oops");
            }

            return;
        }
        
        // if running interactively, don't use the file
        String refFileName = args.length > 0 ? "refpref"+nPoints+"_"+temperature : null;
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
        
        IAction progressReport = new IAction() {
            public void actionPerformed() {
                System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
                System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
            }
        };
        IntegratorListenerAction progressReportListener = new IntegratorListenerAction(progressReport);
        progressReportListener.setInterval((int)(steps/10));
        sim.integratorOS.getEventManager().addListener(progressReportListener);

        sim.integratorOS.getMoveManager().setEquilibrating(false);
        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();

        System.out.println("final reference step frequency "+sim.integratorOS.getStepFreq0());
        
        double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
        System.out.println("ratio average: "+ratioAndError[0]+", error: "+ratioAndError[1]);
        System.out.println("abs average: "+ratioAndError[0]*HSB[nPoints]+", error: "+ratioAndError[1]*HSB[nPoints]);
        DataGroup allYourBase = (DataGroup)sim.accumulators[0].getData(sim.dsvo.minDiffLocation());
        System.out.println("hard sphere ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);
        System.out.println("hard sphere   average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
        System.out.println("hard sphere overlap average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
        
        allYourBase = (DataGroup)sim.accumulators[1].getData(sim.accumulators[1].getNBennetPoints()-sim.dsvo.minDiffLocation()-1);
        System.out.println("chain ratio average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1]);
        System.out.println("chain average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[0]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[0]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[0]);
        System.out.println("chain overlap average: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(AccumulatorAverage.ERROR.index)).getData()[1]);
    }

    /**
     * Inner class for parameters
     */
    public static class VirialAlkaneMixParam extends ParameterBase {
        public int nPoints = 2; //The number of total components = nComp1+nComp2
        public int nSpheres1 = 2;   // The number of C for component1
        public int nSpheres2 = 12;   // The number of C for component2
        public double temperature = 300;   // Kelvin
        public long numSteps = 1000;
        public int nComp1 = 1;//The number of component1
        public int nComp2 = 1;//The number of component2
    }


}
