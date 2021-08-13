/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations.helium;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ANIntergroupExchange;
import etomica.atom.iterator.ANIntragroupExchange;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.ApiIntergroupExchange;
import etomica.config.ConformationLinear;
import etomica.data.IData;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorEvent;
import etomica.integrator.IntegratorListener;
import etomica.math.DoubleRange;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;
import etomica.units.Kelvin;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;
import etomica.virial.cluster.*;
import etomica.virial.mcmove.MCMoveClusterRingPartialRegrow;
import etomica.virial.mcmove.MCMoveClusterRingScale;
import etomica.virial.simulations.SimulationVirialOverlap2;

/**
 * MSMC calculation for pieces of the "odd" exchange component of helium-4 virial coefficients.
 * These pieces are interacting rings containing paths of one or more helium-4 atoms.
 * Length of the rings array indicates the number of rings (and hard spheres in the reference system). 
 * Elements of the rings array indicate the number of atoms in each ring. 
 * 
 * Adapted from VirialHePIXC by Shaul and Schultz
 */
public class VirialHePIXCOdd {

    public static void main(String[] args) {
        VirialHePIParam params = new VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        final int[] rings = params.rings;
        final double temperatureK = params.temperature;
        long blocks = params.blocks;
        int stepsPerBlock = params.stepsPerBlock;
        long blocksEq = params.blocksEq;
        double refFreq = params.refFrac;
        double sigmaHSRef = params.sigmaHSRef;
        final int nBeads = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        final boolean pairOnly = params.pairOnly;
        final boolean doTotal = params.doTotal;
        boolean subtractHalf = params.subtractHalf;
        final int beadFac = subtractHalf ? 2 : 1;
        final boolean calcApprox = params.calcApprox;
        final boolean subtractApprox = !calcApprox && !subtractHalf && params.subtractApprox;
        
        final int nRings = rings.length;
        

        Space space = Space3D.getInstance();

        double heMass = 4.002602;
        
        int nPoints = 0;
        for (int i=0; i<nRings; i++) {
        	nPoints = nPoints + rings[i];
        	System.out.println(rings[i]+"-molecule ring");
        	if  (i>0) {
        		if (rings[i]>rings[i-1]) {
        			throw new RuntimeException("list rings in descending order.");
        		}
        	}
        }
        
        if (sigmaHSRef == -1) {
            // untested 
            sigmaHSRef = 4 + 20/(10+temperatureK);
            sigmaHSRef = sigmaHSRef*nPoints/nRings;
        }
        
        final double[] HSB = new double[8];
        
        
        System.out.println("He path integral ("+nBeads+"-mer chains) B"+nPoints+"Odd component at "+temperatureK+" K");
        
        
        if (subtractHalf) {
            System.out.println("Calculating difference between "+nBeads/beadFac+" and "+nBeads+" beads");
        } else {
            if (subtractApprox) {
                 System.out.println("computing difference from approximate He");
            }
        }
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (pairOnly) {
            System.out.println("computing pairwise contribution");
        } else {
            System.out.println("computing non-additive contribution");
        }
        if (doTotal) {
        	throw new RuntimeException("doTotal not available yet.");
        }
        final double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        HSB[2] = Standard.B2HS(sigmaHSRef);
        HSB[3] = Standard.B3HS(sigmaHSRef);
        HSB[4] = Standard.B4HS(sigmaHSRef);
        HSB[5] = Standard.B5HS(sigmaHSRef);
        HSB[6] = Standard.B6HS(sigmaHSRef);
        HSB[7] = Standard.B7HS(sigmaHSRef);
        MayerHardSphere fRef = new MayerHardSphere(sigmaHSRef);
        MayerEHardSphere eRef = new MayerEHardSphere(sigmaHSRef);
        ClusterAbstract refCluster = Standard.virialCluster(nRings, fRef, nPoints>3, eRef, true);
        double refIntegral = HSB[nRings];
        
        P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2HeSimplified p2Approx = new P2HeSimplified(space);
        Potential p3 = new P3CPSNonAdditiveHe(space);
        final P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);
        
        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2, new ApiIntergroupExchange(nBeads));
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }
        

        PotentialGroupPI pTargetApproxGroup = new PotentialGroupPI(beadFac);
        pTargetApproxGroup.addPotential(p2Approx, new ApiIntergroupExchange(nBeads));
        PotentialGroupPISkip[] pTargetApproxSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetApproxSkip[i] = pTargetApproxGroup.new PotentialGroupPISkip(i);
        }
        
        // p3 masquerading as a potential between 2 molecules
        PotentialGroupPI p3TargetGroup = new PotentialGroupPI(beadFac);
        p3TargetGroup.addPotential(p3, new ANIntergroupExchange(nBeads));
        p3TargetGroup.addPotential(p2, new ApiIntergroupExchange(nBeads));
        PotentialGroupPISkip[] p3TargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3TargetGroup.new PotentialGroupPISkip(i);
        }
        
        PotentialGroupPI p3TargetApproxGroup = new PotentialGroupPI(beadFac);
        p3TargetApproxGroup.addPotential(p3Approx, new ANIntergroupExchange(nBeads));
        p3TargetApproxGroup.addPotential(p2Approx, new ApiIntergroupExchange(nBeads));
        PotentialGroupPISkip[] p3TargetApproxSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetApproxSkip[i] = p3TargetApproxGroup.new PotentialGroupPISkip(i);
        }

        MayerGeneral fTarget = new MayerGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerGeneral fTargetApprox = new MayerGeneral(pTargetApproxGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        MayerGeneral[] fTargetSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            fTargetSkip[i] = new MayerGeneral(pTargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        
        MayerGeneral f3Target = new MayerGeneral(p3TargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };
        
        MayerGeneral f3TargetApprox = new MayerGeneral(p3TargetApproxGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return super.f(pair, r2, beta/nBeads);
            }
        };

        MayerGeneral[] f3TargetSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            f3TargetSkip[i] = new MayerGeneral(p3TargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }

        MayerFunctionThreeBody e3Target = new MayerFunctionMolecularThreeBody(p3TargetGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return super.f(molecules, r2, beta/nBeads) + 1;
            }
            public double calcF(double x) {
            	return Math.exp(x);
            }
        };
        
        int nDiagrams = 1;
        if (nRings == 3 && pairOnly) {
        	nDiagrams = 4;
        }
    	ClusterBonds[] clusters = new ClusterBonds[nDiagrams];
    	double[] weights = new double[nDiagrams];
    	
    	if (nRings == 2) {
    		// works for B3Odd and B4Odd diagrams with two clusters
    		int[][][] bondList = {{{0,1}}, {}};
            clusters[0] = new ClusterBonds(2, bondList, false);
            weights[0] = 1.0;
    	} else if (nRings == 3) {
    		//attempt at B4Odd diagram with 1 dimer and 2 monomers
    		if (pairOnly) {
    			int[][][] bondList0 = {{{0,1},{0,2},{1,2}}, {}};
    			int[][][] bondList1 = {{{0,1},{0,2}}, {}};
    			int[][][] bondList2 = {{{0,1},{1,2}}, {}};
    			int[][][] bondList3 = {{{0,2},{1,2}}, {}};
                clusters[0] = new ClusterBonds(3, bondList0, false);
                clusters[1] = new ClusterBonds(3, bondList1, false);
                clusters[2] = new ClusterBonds(3, bondList2, false);
                clusters[3] = new ClusterBonds(3, bondList3, false);
                weights[0] = 1.0;
                weights[1] = 1.0;
                weights[2] = 1.0;
                weights[3] = 1.0;
    		} else {
    			clusters[0] = new ClusterBondsNonAdditive(3, new int[1][0][0], new int[][]{{},{},{},{0}});
    			weights[0] = 1.0;
    		}
            
    	} else {
        	throw new RuntimeException("Error in rings array.");
        }
    	

    	double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*heMass*temperature);

    	ClusterAbstract targetCluster; 
    	
    	if (pairOnly) {
    		if (calcApprox) {
    			if (subtractHalf){
    				throw new RuntimeException("TBD");
    			}
    			targetCluster = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetApprox});
    		} else if (subtractHalf) {
				ClusterAbstract fullCluster = new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
				ClusterAbstract[] halfClusters = new ClusterAbstract[2];
    			halfClusters[0] = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetSkip[0]});
    			halfClusters[1] = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetSkip[1]});
    			targetCluster = new ClusterDifference(fullCluster, halfClusters);
			} else if (subtractApprox){
				ClusterAbstract[] approxCluster = new ClusterAbstract[1];
				approxCluster[0] =  new ClusterSum(clusters,weights,new MayerFunction[]{fTargetApprox});
    			ClusterAbstract totalCluster =  new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
    			targetCluster = new ClusterDifference(totalCluster, approxCluster);
			} else {
				targetCluster = new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
			}	
    	} else {  //nonadditive contribution
    		if (nRings == 2) {
    			ClusterAbstract[] addCluster = new ClusterAbstract[1];
        		if (calcApprox) {
        			if (subtractHalf){
        				throw new RuntimeException("TBD");
        			}
        			addCluster[0] =  new ClusterSum(clusters,weights,new MayerFunction[]{fTargetApprox});
        			ClusterAbstract totalCluster =  new ClusterSum(clusters,weights,new MayerFunction[]{f3TargetApprox});
        			targetCluster = new ClusterDifference(totalCluster, addCluster);
        			
        		} else if (subtractHalf) {
        			
        			ClusterAbstract addClusterSH =  new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
        			ClusterAbstract[] halfClustersAdd = new ClusterAbstract[2];
        			halfClustersAdd[0] = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetSkip[0]});
        			halfClustersAdd[1] = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetSkip[1]});
        			addCluster[0] = new ClusterDifference(addClusterSH, halfClustersAdd);
        			
        			ClusterAbstract totalClusterSH =  new ClusterSum(clusters,weights,new MayerFunction[]{f3Target});
                	ClusterAbstract[] halfClustersTotal = new ClusterAbstract[2];
        			halfClustersTotal[0] = new ClusterSum(clusters,weights,new MayerFunction[]{f3TargetSkip[0]});
        			halfClustersTotal[1] = new ClusterSum(clusters,weights,new MayerFunction[]{f3TargetSkip[1]});
        			ClusterAbstract totalCluster = new ClusterDifference(totalClusterSH, halfClustersTotal);
        			
        			targetCluster = new ClusterDifference(totalCluster, addCluster);
        			
        		} else if (subtractApprox) {
        			
        			ClusterAbstract totalClusterApprox =  new ClusterSum(clusters,weights,new MayerFunction[]{f3TargetApprox});
        			ClusterAbstract[] addClusterApprox = new ClusterAbstract[2];
        			addClusterApprox[0] = new ClusterSum(clusters,weights,new MayerFunction[]{fTargetApprox});
        			ClusterAbstract[] approxCluster = new ClusterAbstract[1];
        			approxCluster[0] = new ClusterDifference(totalClusterApprox, addClusterApprox);
        			
        			ClusterAbstract totalCluster1 =  new ClusterSum(clusters,weights,new MayerFunction[]{f3Target});
                	ClusterAbstract[] addCluster1 = new ClusterAbstract[1];
                	addCluster1[0] = new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
        			ClusterAbstract totalCluster = new ClusterDifference(totalCluster1, addCluster1);
        			
        			targetCluster = new ClusterDifference(totalCluster, approxCluster);
        			
        		} else {
        		
                	addCluster[0] =  new ClusterSum(clusters,weights,new MayerFunction[]{fTarget});
                	ClusterAbstract totalCluster =  new ClusterSum(clusters,weights,new MayerFunction[]{f3Target});
                	targetCluster = new ClusterDifference(totalCluster, addCluster);
            		
        		}
    		} else if (nRings == 3) {
    			if (calcApprox) {
    				throw new RuntimeException("TBD");
        		} else if (subtractHalf) {
        			throw new RuntimeException("TBD");
        		} else if (subtractApprox) {
        			throw new RuntimeException("TBD");
        		} else {
        		
        			ClusterAbstract totalCluster = new ClusterSumMultibody(clusters, new double[]{1.0}, new MayerFunction[]{fTarget}, new MayerFunctionNonAdditive[]{e3Target});
        			ClusterAbstract[] addCluster = new ClusterAbstract[1];
        	    	addCluster[0] =  new ClusterSum(new ClusterBonds[]{new ClusterBonds(3, new int[][][]{{{0,1},{0,2},{2,1}}})}, new double[]{1.0}, new MayerFunction[]{eTarget});
        	    	targetCluster = new ClusterDifference(totalCluster, addCluster);
            		
        		}
    			
    		} else {
    			throw new RuntimeException("Error in rings array.");
    		}

    	}
        
        // the cluster's temperature determines the factor multiplied in the exponential (f=e-1)
        // we want 1/(P*kT)
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        
        System.out.println("sigmaHSRef: "+sigmaHSRef);
        // overerr expects this string, BnHS
        System.out.println("B"+nRings+"HS: "+refIntegral);
        long steps = stepsPerBlock*blocks;
        if (steps%1000 != 0) {
            throw new RuntimeException("steps should be a multiple of 1000");
        }
        System.out.println(steps+" steps ("+blocks+" blocks of "+stepsPerBlock+" steps)");
        System.out.println(1000+" steps per overlap-sampling block");
        
        SpeciesGeneral[] species = new SpeciesGeneral[nRings];
        int[] simMolecules = new int[nRings];
        for (int i=0;i<nRings;i++) {
            species[i] = new SpeciesBuilder(space)
                    .withConformation(new ConformationLinear(space, 0))
                    .addCount(AtomType.simple("He" + rings[i] + i, heMass), rings[i] * nBeads)
                    .build();
            simMolecules[i] = 1;
        }

        // new int[]{nPoints+(doFlex?1:0)} needed at fourth (but not third) order
        final SimulationVirialOverlap2 sim = new SimulationVirialOverlap2(space, species, simMolecules, temperature, refCluster, targetCluster);
        sim.setDoFasterer(true);
        sim.init();

        //Garberoglio and Harvey's U12
        Potential2SoftSpherical p2B = new P2HePCKLJS(space) {
        	public double energy(IAtomList molecules) {
           	 	return (super.energy(molecules))/nBeads;
           	}
        };
        PotentialGroup uXC = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
        uXC.addPotential(p2B, new ANIntragroupExchange(2, nBeads));
        // integrators share a common potentialMaster.  so just add to one
        for (int i=0;i<nRings;i++) {
        	if (rings[i]>1) {
        		sim.integrators[1].getPotentialMaster().addPotential(uXC,new ISpecies[]{species[i]});
        	}
        }
       
        
        // we'll use substeps=1000 initially (to allow for better initialization)
        // and then later switch to 1000 overlap steps
        sim.integratorOS.setNumSubSteps(1000);
        steps /= 1000;

        
        // rotation is a bit pointless when we can regrow the chain completely
        sim.integrators[0].getMoveManager().removeMCMove(sim.mcMoveRotate[0]);
        sim.integrators[1].getMoveManager().removeMCMove(sim.mcMoveRotate[1]);
        
        System.out.println("regrow full ring");
        
                
        double energyFac = nBeads*Math.PI/(lambda*lambda);
        
        MCMoveClusterRingPartialRegrow[] ringMove = new MCMoveClusterRingPartialRegrow[2];
        MCMoveClusterRingScale[] ringScale = new MCMoveClusterRingScale[2];
        int numCBMCTrials = nBeads/20+20; 
        int numRegrowBeads = nBeads/3;
        System.out.println("regrow "+numRegrowBeads+" beads ("+numCBMCTrials+" trials)");
        for (int i=0;i<2;i++) {
        	ringMove[i] = new MCMoveClusterRingPartialRegrow(sim.integrators[1].getPotentialMaster(),sim.getRandom(), space);
            ringMove[i].setEnergyFactor(energyFac);
            ringMove[i].setNumTrial(numCBMCTrials);
            ringMove[i].setNumBeads(numRegrowBeads);
            sim.integrators[i].getMoveManager().addMCMove(ringMove[i]);
            ringScale[i] = new MCMoveClusterRingScale(sim.integrators[1].getPotentialMaster(),sim.getRandom(), space);
            ringScale[i].setEnergyFactor(energyFac);
            sim.integrators[i].getMoveManager().addMCMove(ringScale[i]);
            sim.integrators[i].getMoveManager().setFrequency(ringScale[i], 0.01);
        }
        
        if (refFreq >= 0) {
            sim.integratorOS.setAdjustStepFraction(false);
            sim.integratorOS.setRefStepFraction(refFreq);
        }
        
        //Ensure that beads are not initialized in overlapped configuration.
        //Exchange interaction would yield sampling weight of zero.
        double r = 2;
        for (int j=0;j<nRings;j++) {
        	IAtomList leafList = sim.box[0].getMoleculeList(species[j]).get(0).getChildList();
            for (int i = 0; i<leafList.size(); i++) {
                Vector p = leafList.get(i).getPosition();
                p.setX(0, r*Math.cos((2*Math.PI*i)/leafList.size()));
                p.setX(1, r*Math.sin((2*Math.PI*i)/leafList.size()));
            }
            IAtomList leafList1 = sim.box[1].getMoleculeList(species[j]).get(0).getChildList();
            for (int i = 0; i<leafList1.size(); i++) {
                Vector p = leafList1.get(i).getPosition();
                p.setX(0, j*10+r*Math.cos((2*Math.PI*i)/leafList1.size()));
                p.setX(1, r*Math.sin((2*Math.PI*i)/leafList1.size()));
            }
        }
        

       
        if (false) {
            // unnecessary because our MC move regrows the chain using the
            // probability distribution appropriate for the harmonic bonds
            
            // create the intramolecular potential here, add to it and add it to
            // the potential master if needed
            PotentialGroup pIntra = sim.integrators[1].getPotentialMaster().makePotentialGroup(1);
            // we want exp[-(pi*P/lambda^2) * sum(x^2)]
            // we set the integrator temperature=1 above, so when it does
            //   exp[-beta * U] = exp[-U]
            // so just make the spring constant whatever we need to get the above expression
            P2Harmonic p2Bond = new P2Harmonic(space, 2*Math.PI*nBeads/(lambda*lambda)*temperature);
            int[][] pairs = new int[nBeads][2];
            for (int i=0; i<nBeads-1; i++) {
                pairs[i][0] = i;
                pairs[i][1] = i+1;
            }
            pairs[nBeads-1][0] = nBeads-1;
            pairs[nBeads-1][1] = 0;
            pIntra.addPotential(p2Bond, new ApiIndexList(pairs));
            // integrators share a common potentialMaster.  so just add to one
            sim.integrators[1].getPotentialMaster().addPotential(pIntra,new ISpecies[]{sim.getSpecies(0)});
        }
/*
        if (false) {
            double vSize = 5;
            sim.box[0].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
            sim.box[1].getBoundary().setBoxSize(space.makeVector(new double[]{vSize,vSize,vSize}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, space, sim.getController());
            DisplayBox displayBox0 = simGraphic.getDisplayBox(sim.box[0]); 
            DisplayBox displayBox1 = simGraphic.getDisplayBox(sim.box[1]);
            displayBox0.setPixelUnit(new Pixel(300.0/vSize));
            displayBox1.setPixelUnit(new Pixel(300.0/vSize));
            displayBox0.setShowBoundary(false);
            displayBox1.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox0.canvas).setBackgroundColor(Color.WHITE);
            ((DisplayBoxCanvasG3DSys)displayBox1.canvas).setBackgroundColor(Color.WHITE);
            
            
            IAtomType type = species.getLeafType();
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox0.getDiameterHash();
            diameterManager.setDiameter(type, 0.02+1.0/nBeads);
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
//            sim.getController().removeAction(sim.ai);
//            sim.getController().addAction(new IAction() {
//                public void actionPerformed() {
//                    sim.initRefPref(null, 10, false);
//                    sim.equilibrate(null, 20);
//                    sim.ai.setMaxSteps(Long.MAX_VALUE);
//                }
//            });
//            sim.getController().addAction(sim.ai);
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
                public void actionPerformed() {
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    data.x = ratio;
                    averageBox.putData(data);
                    data.x = error;
                    errorBox.putData(data);
                }
                
                DataDouble data = new DataDouble();
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
        */
        
        // if running interactively, don't use the file
        String refFileName = null;
        if (isCommandline) {
            String tempString = ""+temperatureK;
            if (temperatureK == (int)temperatureK) {
                // temperature is an integer, use "200" instead of "200.0"
                tempString = ""+(int)temperatureK;
            }
            refFileName = "refpref"+nPoints;
            refFileName += pairOnly ? "_2b" : "_3b";
            refFileName += "_"+tempString+"_"+nBeads;
            if (subtractHalf) {
                refFileName += "_sh";
            }
           
        }
        long t1 = System.currentTimeMillis();
        // this will either read the refpref in from a file or run a short simulation to find it
        sim.initRefPref(refFileName, steps/40);
        /*
        // this can't be done after equilibration.  ClusterSumShell needs at least
        // one accepted move before it can collect real data.  we'll reset below
        MeterVirial meterDiagrams = new MeterVirial(targetDiagrams);
        meterDiagrams.setBox(sim.box[1]);
        AccumulatorAverageCovariance accumulatorDiagrams = null;
        if (targetDiagrams.length > 0) {
            // if we have 1 diagram, we don't need covariance, but it won't actually cause problems
            accumulatorDiagrams = new AccumulatorAverageCovariance(steps);
            DataPumpListener pumpDiagrams = new DataPumpListener(meterDiagrams, accumulatorDiagrams);
            sim.integrators[1].getEventManager().addListener(pumpDiagrams);
        }*/

        // run another short simulation to find MC move step sizes and maybe narrow in more on the best ref pref
        // if it does continue looking for a pref, it will write the value to the file
        sim.equilibrate(refFileName, blocksEq*stepsPerBlock/1000);
ActivityIntegrate ai = new ActivityIntegrate(sim.integratorOS, steps);
sim.setAccumulatorBlockSize(stepsPerBlock);
        //sim.integratorOS.setNumSubSteps((int)steps);
        //sim.integratorOS.setAgressiveAdjustStepFraction(true);

        System.out.println("equilibration finished");
        System.out.println("MC Move step sizes (ref)    "+sim.mcMoveTranslate[0].getStepSize());
        System.out.println("MC Move step sizes (target) "+sim.mcMoveTranslate[1].getStepSize());


        final HistogramNotSoSimple targHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple targPiHist = new HistogramNotSoSimple(70, new DoubleRange(-1, 8));
        final HistogramNotSoSimple hist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final HistogramNotSoSimple piHist = new HistogramNotSoSimple(100, new DoubleRange(0, sigmaHSRef));
        final ClusterAbstract finalTargetCluster = targetCluster.makeCopy();
        IntegratorListener histListenerRef = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                CoordinatePairSet cPairs = sim.box[0].getCPairSet();
                for (int i=0; i<2; i++) {
                    for (int j=i+1; j<2; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }
                double v = finalTargetCluster.value(sim.box[0]);
                hist.addValue(Math.sqrt(r2Max), v);
                piHist.addValue(Math.sqrt(r2Max), Math.abs(v));
            }

            public void integratorInitialized(IntegratorEvent e) {
            }
        };
        IntegratorListener histListenerTarget = new IntegratorListener() {
            public void integratorStepStarted(IntegratorEvent e) {}

            public void integratorStepFinished(IntegratorEvent e) {
                double r2Max = 0;
                double r2Min = Double.POSITIVE_INFINITY;
                CoordinatePairSet cPairs = sim.box[1].getCPairSet();
                for (int i=0; i<nRings; i++) {
                    for (int j=i+1; j<nRings; j++) {
                        double r2ij = cPairs.getr2(i, j);
                        if (r2ij < r2Min) r2Min = r2ij;
                        if (r2ij > r2Max) r2Max = r2ij;
                    }
                }

                double v = finalTargetCluster.value(sim.box[1]);
                double r = Math.sqrt(r2Max);
                if (r > 1) {
                    r = Math.log(r);
                }
                else {
                    r -= 1;
                }
                targHist.addValue(r, v);
                targPiHist.addValue(r, Math.abs(v));
            }

            public void integratorInitialized(IntegratorEvent e) {}
        };
        if (!isCommandline) {
            // if interactive, print intermediate results
            final double refIntegralF = refIntegral;
            IntegratorListener progressReport = new IntegratorListener() {
                public void integratorInitialized(IntegratorEvent e) {}
                public void integratorStepStarted(IntegratorEvent e) {}
                public void integratorStepFinished(IntegratorEvent e) {
                    if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                    System.out.print(sim.integratorOS.getStepCount()+" steps: ");
                    double[] ratioAndError = sim.dvo.getAverageAndError();
                    double ratio = ratioAndError[0];
                    double error = ratioAndError[1];
                    System.out.println("abs average: "+ratio*refIntegralF+" error: "+error*refIntegralF);
                    if (ratio == 0 || Double.isNaN(ratio)) {
                        throw new RuntimeException("oops");
                    }
                }
            };
            sim.integratorOS.getEventManager().addListener(progressReport);
            if (params.doHist) {
                IntegratorListener histReport = new IntegratorListener() {
                    public void integratorInitialized(IntegratorEvent e) {}
                    public void integratorStepStarted(IntegratorEvent e) {}
                    public void integratorStepFinished(IntegratorEvent e) {
                        if ((sim.integratorOS.getStepCount()*10) % ai.getMaxSteps() != 0) return;
                        System.out.println("**** reference ****");
                        double[] xValues = hist.xValues();
                        double[] h = hist.getHistogram();
                        double[] piH = piHist.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (!Double.isNaN(h[i])) {
                                System.out.println(xValues[i]+" "+h[i]+" "+piH[i]);
                            }
                        }
                        System.out.println("**** target ****");
                        xValues = targHist.xValues();
                        h = targHist.getHistogram();
                        piH = targPiHist.getHistogram();
                        for (int i=0; i<xValues.length; i++) {
                            if (!Double.isNaN(h[i])) {
                                double r = xValues[i];
                                if (r < 0) r += 1;
                                else r = Math.exp(r);
                                System.out.println(r+" "+h[i]+" "+piH[i]);
                            }
                        }
                    }
                };
                sim.integratorOS.getEventManager().addListener(histReport);
            }

        }
        if (params.doHist) {
            System.out.println("collecting histograms");
            // only collect the histogram if we're forcing it to run the reference system
            sim.integrators[0].getEventManager().addListener(histListenerRef);
            sim.integrators[1].getEventManager().addListener(histListenerTarget);
        }
sim.getController().runActivityBlocking(ai);
        /*
        if (accumulatorDiagrams != null) {
            accumulatorDiagrams.reset();
        }
*/
        // make the accumulator block size equal to the # of steps performed for each overlap step.
        // make the integratorOS aggressive so that it runs either reference or target
        // then, we'll have some number of complete blocks in the accumulator
        long t2 = System.currentTimeMillis();
        
        if (params.doHist) {
            double[] xValues = hist.xValues();
            double[] h = hist.getHistogram();
            for (int i=0; i<xValues.length; i++) {
                if (!Double.isNaN(h[i])) {
                    System.out.println(xValues[i]+" "+(-2*h[i]+1));
                }
            }
        }

        System.out.println("final reference step fraction "+sim.integratorOS.getIdealRefStepFraction());
        System.out.println("actual reference step fraction "+sim.integratorOS.getRefStepFraction());
        

        sim.printResults(refIntegral);

        DataGroup allData = (DataGroup)sim.accumulators[1].getData();
        IData dataAvg = allData.getData(sim.accumulators[1].AVERAGE.index);
        IData dataErr = allData.getData(sim.accumulators[1].ERROR.index);
        IData dataCov = allData.getData(sim.accumulators[1].BLOCK_COVARIANCE.index);
        // we'll ignore block correlation -- whatever effects are here should be in the full target results
        
        
        System.out.println();

        DataGroup allData0 = (DataGroup)sim.accumulators[0].getData();
        IData dataAuto0 = allData0.getData(sim.accumulators[0].BLOCK_CORRELATION.index);
        System.out.println("reference autocorrelation function: "+dataAuto0.getValue(0));
        System.out.println("reference overlap autocorrelation function: "+dataAuto0.getValue(1));
        
        System.out.println();

        IData dataAuto = allData.getData(sim.accumulators[1].BLOCK_CORRELATION.index);
        System.out.println("target autocorrelation function: "+dataAuto.getValue(0));
        System.out.println("target overlap autocorrelation function: "+dataAuto.getValue(1));
        
        System.out.println();
        
        System.out.println("time: "+(t2-t1)/1000.0);
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
    public static class VirialHePIParam extends ParameterBase {
        public int[] rings = {2,1,1}; //length of array indicates number of rings (and hard spheres); elements indicate number of molecules in each ring
        public int nBeads = 16;
        public double temperature = 7;   // Kelvin
        public long blocks = 1000;  //NOT overlap blocks
        public int stepsPerBlock = 1000;
        public long blocksEq=10; //NOT overlap steps
        public double refFrac = 0.5; //not adjustment of step freqency if positive
        public boolean doHist = false;
        public double sigmaHSRef = -1; // -1 means use equation for sigmaHSRef
        public int beadFac = 1;  //2 if subtract half, 1 if not
        public boolean subtractHalf = false;
        public boolean calcApprox = true;
        public boolean subtractApprox = false;
        public boolean pairOnly = true;
        public boolean doTotal = false;

    }
}
