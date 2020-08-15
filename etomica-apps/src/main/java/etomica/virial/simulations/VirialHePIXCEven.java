/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.simulations;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.iterator.ANIntergroupCoupled;
import etomica.atom.iterator.ApiIntergroupCoupled;
import etomica.chem.elements.ElementChemical;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorRatioAverageCovariance;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataGroup;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMoleculeList;
import etomica.potential.*;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheres;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.dimensions.*;
import etomica.units.dimensions.Dimension;
import etomica.util.Constants;
import etomica.util.Constants.CompassDirection;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;
import etomica.virial.*;
import etomica.virial.PotentialGroup3PI.PotentialGroup3PISkip;
import etomica.virial.PotentialGroupPI.PotentialGroupPISkip;

import javax.swing.*;
import java.awt.*;

/**
 * MSMC calculation for the "even" exchange component of helium-4 virial coefficients.
 * These pieces are single rings containing paths of two or more helium-4 atoms.
 * 
 * Direct sampling with sampling weight determined by reference diagram.
 * 
 * For additive components, reference value is the same diagram for an ideal gas 
 * (or system with some measure of IG character.)
 * 
 * For nonadditive components, reference value is the additive component.
 * 
 * Adapted from VirialHePIXC by Shaul and Schultz
 */
public class VirialHePIXCEven {

    public static void main(String[] args) {
        VirialHePIParam params = new VirialHePIParam();
        boolean isCommandline = args.length > 0;
        if (isCommandline) {
            ParseArgs parseArgs = new ParseArgs(params);
            parseArgs.parseArgs(args, true);
        }
        final int nPoints = params.nPoints;
        double temperatureK = params.temperature;
        final long steps = params.numSteps;
        final double aRef = params.aRef;
        final double aTarget = params.aTarget;
        if (aTarget > aRef) {
            throw new RuntimeException("aTarget should be less than aRef");
        }
        final int nBeads = (params.nBeads > -1) ? params.nBeads : ((int)(1200/temperatureK) + 7);
        Space space = Space3D.getInstance();
        final boolean pairOnly = params.nPoints == 2 || params.pairOnly;
        final boolean doTotal = params.doTotal;
        boolean subtractHalf = params.subtractHalf;
        final int beadFac = subtractHalf ? 2 : 1;
        final boolean calcApprox = params.calcApprox;
        final boolean subtractApprox = !calcApprox && !subtractHalf && params.subtractApprox;
        

        if (subtractHalf) {
            System.out.println("He Path Integral ("+nBeads+"-mer chains) B"+nPoints+"Even at "+temperatureK+" K");
            System.out.println("Calculating difference between "+nBeads/beadFac+" and "+nBeads+" beads");
            if (subtractApprox) {
            	throw new RuntimeException("calcApprox or subtractApprox only");
            }
        }
        else {
            System.out.println("He Path Integral ("+nBeads+"-mer chains) B"+nPoints+"Even at "+temperatureK+" K");
            if (subtractApprox) {  
                 System.out.println("computing difference from approximate He");
            }
        }
        System.out.println("perturbing from a = "+aRef+" to a = "+aTarget);
        if (calcApprox) System.out.println("Calculating coefficients for approximate potential");
        if (pairOnly) {
            System.out.println("computing pairwise contribution");
        }
        else {
            System.out.println("computing non-additive contribution");
        }
        double temperature = Kelvin.UNIT.toSim(temperatureK);
        
        P2HePCKLJS p2 = new P2HePCKLJS(space);
        P2HeSimplified p2Approx = new P2HeSimplified(space);
        Potential p3 = new P3CPSNonAdditiveHe(space);
        final P3CPSNonAdditiveHeSimplified p3Approx = new P3CPSNonAdditiveHeSimplified(space);
        p3Approx.setParameters(temperatureK);
        
        PotentialGroupPI pTargetGroup = new PotentialGroupPI(beadFac);
        pTargetGroup.addPotential(p2, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetSkip[i] = pTargetGroup.new PotentialGroupPISkip(i);
        }

        PotentialGroupPI pTargetApproxGroup = new PotentialGroupPI(beadFac);
        pTargetApproxGroup.addPotential(p2Approx, new ApiIntergroupCoupled());
        PotentialGroupPISkip[] pTargetApproxSkip = new PotentialGroupPISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            pTargetApproxSkip[i] = pTargetApproxGroup.new PotentialGroupPISkip(i);
        }
        
        PotentialGroup3PI p3TargetGroup = new PotentialGroup3PI(beadFac);
        p3TargetGroup.addPotential(p3, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }
        
        PotentialGroup3PI p3TargetApproxGroup = new PotentialGroup3PI(beadFac);
        p3TargetApproxGroup.addPotential(p3Approx, new ANIntergroupCoupled(3));
        PotentialGroup3PISkip[] p3TargetApproxSkip = new PotentialGroup3PISkip[beadFac];
        for (int i=0; i<beadFac; i++) {
            p3TargetApproxSkip[i] = p3TargetGroup.new PotentialGroup3PISkip(i);
        }
        
        MayerEGeneral eRef = new MayerEGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return aRef + (1-aRef)*super.f(pair, r2, beta/nBeads);
            }
        };
        
        MayerEGeneral eTarget = new MayerEGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return aTarget + (1-aTarget)*super.f(pair, r2, beta/nBeads);
            }
        };
        
        MayerGeneral fTarget = new MayerGeneral(pTargetGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return aTarget + (1-aTarget)*super.f(pair, r2, beta/nBeads);
            }
        };
        
        MayerGeneral fTargetApprox = new MayerGeneral(pTargetApproxGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return aTarget + (1-aTarget)*super.f(pair, r2, beta/nBeads);
            }
        };
        
        MayerEGeneral eTargetApprox = new MayerEGeneral(pTargetApproxGroup) {
            public double f(IMoleculeList pair, double r2, double beta) {
                return aTarget + (1-aTarget)*super.f(pair, r2, beta/nBeads);
            }
        };
       
        MayerEGeneral[] eTargetSkip = new MayerEGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            eTargetSkip[i] = new MayerEGeneral(pTargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        
        MayerEGeneral[] eTargetApproxSkip = new MayerEGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            eTargetApproxSkip[i] = new MayerEGeneral(pTargetApproxSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        
        MayerGeneral[] fTargetSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            fTargetSkip[i] = new MayerGeneral(pTargetSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        
        MayerGeneral[] fTargetApproxSkip = new MayerGeneral[beadFac];
        for (int i=0; i<beadFac; i++) {
            fTargetApproxSkip[i] = new MayerGeneral(pTargetApproxSkip[i]) {
                public double f(IMoleculeList pair, double r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(pair, r2, beta/(nBeads/beadFac));
                }
            };
        }
        
        MayerFunctionThreeBody e3Target = new MayerFunctionMolecularThreeBody(p3TargetGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return aTarget + (1-aTarget)*super.f(molecules, r2, beta/nBeads) + 1;
            }
            public double calcF(double x) {
            	return Math.exp(x);
            }
        };
        
        MayerFunctionThreeBody[] e3TargetSkip = new MayerFunctionMolecularThreeBody[beadFac];
        for (int i=0; i<beadFac; i++) {
            e3TargetSkip[i] = new MayerFunctionMolecularThreeBody(p3TargetSkip[i]) {
                public double f(IMoleculeList molecules, double[] r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(molecules, r2, beta/(nBeads/beadFac)) + 1;
                }
                public double calcF(double x) {
                	return Math.exp(x);
                }
            };
        }
        MayerFunctionThreeBody e3TargetApprox = new MayerFunctionMolecularThreeBody(p3TargetApproxGroup) {
            public double f(IMoleculeList molecules, double[] r2, double beta) {
                return aTarget + (1-aTarget)*super.f(molecules, r2, beta/nBeads) + 1;
            }
            public double calcF(double x) {
            	return Math.exp(x);
            }
        };
        
        MayerFunctionThreeBody[] e3TargetApproxSkip = new MayerFunctionMolecularThreeBody[beadFac];
        for (int i=0; i<beadFac; i++) {
            e3TargetSkip[i] = new MayerFunctionMolecularThreeBody(p3TargetApproxSkip[i]) {
                public double f(IMoleculeList molecules, double[] r2, double beta) {
                    return aTarget + (1-aTarget)*super.f(molecules, r2, beta/(nBeads/beadFac)) + 1;
                }
                public double calcF(double x) {
                	return Math.exp(x);
                }
            };
        }

        ClusterSum refCluster; ClusterAbstract targetCluster; ClusterWeight samplingCluster; 
        ClusterBonds[] clusters = new ClusterBonds[1];
        ClusterBonds[] clustersAdd = new ClusterBonds[1];
    	double[] weights = new double[1];
    	double[] weightsAdd = new double[1];
    	
    	if (nPoints == 2) {
    		int[][][] bondList = {{{0,1}}, {}};
            clusters[0] = new ClusterBonds(2, bondList, false);
            weights[0] = 1.0;
    	} else if (nPoints == 3) {
    		clusters[0] = new ClusterBonds(3, new int[][][]{{{0,1},{0,2},{2,1}}});
    		weights[0] = 1.0;
    		if (!pairOnly) {
    			clustersAdd[0] = clusters[0];
    			weightsAdd[0] = weights[0];
    			clusters[0] = new ClusterBondsNonAdditive(3, new int[1][0][0], new int[][]{{},{},{},{0}});
    		}
    	} else if (nPoints == 4) {
    		clusters[0] = new ClusterBonds(4, new int[][][]{{{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}});
    		weights[0] = 1.0;
    		if (!pairOnly) {
    			clustersAdd[0] = clusters[0];
    			weightsAdd[0] = weights[0];
    			clusters[0] = new ClusterBondsNonAdditive(4, new int[1][0][0], new int[][]{{},{},{},{0}});
    		}
    	} else {
        	throw new RuntimeException("Cannot consider that order yet.");
        }
    	
    	if (pairOnly || nPoints == 2) {
			if (calcApprox) {
				if (subtractHalf) {
					ClusterAbstract fullCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApprox});
	    			ClusterAbstract[] halfClusters = new ClusterAbstract[2];
	    			halfClusters[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApproxSkip[0]});
	    			halfClusters[1] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApproxSkip[1]});
	    			targetCluster = new ClusterDifference(fullCluster, halfClusters);
				} else {
					targetCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApprox});
				}
			} else {
				if (subtractHalf) {
					ClusterAbstract fullCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eTarget});
	    			ClusterAbstract[] halfClusters = new ClusterAbstract[2];
	    			halfClusters[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetSkip[0]});
	    			halfClusters[1] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetSkip[1]});
	    			targetCluster = new ClusterDifference(fullCluster, halfClusters);
				} else if (subtractApprox){
					ClusterAbstract fullCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eTarget});
	    			ClusterAbstract[] approxCluster = new ClusterAbstract[1];
	    			approxCluster[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApprox});
	    			targetCluster = new ClusterDifference(fullCluster, approxCluster);
				} else {
					targetCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eTarget});
				}
			}
        	refCluster = new ClusterSum(clusters, weights, new MayerFunction[]{eRef});
    	} else {
    		if (doTotal) {
    			throw new RuntimeException("Haven't set up total calculation yet.");
    		} else {
    			ClusterAbstract[] addCluster = new ClusterAbstract[1];
    			ClusterAbstract totalCluster;
    			if (calcApprox) {
    				totalCluster = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetApprox}, new MayerFunctionNonAdditive[]{e3TargetApprox});
        	    	addCluster[0] =  new ClusterSum(clustersAdd, weightsAdd, new MayerFunction[]{eTargetApprox}); 
    				if (subtractHalf) {
    					ClusterAbstract[] halfTotalClusters = new ClusterAbstract[2];
    	    			halfTotalClusters[0] = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetApproxSkip[0]}, new MayerFunctionNonAdditive[]{e3TargetApproxSkip[0]});
    	    			halfTotalClusters[1] = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetApproxSkip[1]}, new MayerFunctionNonAdditive[]{e3TargetApproxSkip[1]});
    	    			
    	    			ClusterAbstract[] halfAddClusters = new ClusterAbstract[2];
    	    			halfAddClusters[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApproxSkip[0]});
    	    			halfAddClusters[1] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApproxSkip[1]});
    	    			
    	    			totalCluster = new ClusterDifference(totalCluster, halfTotalClusters);
    	    			addCluster[0] = new ClusterDifference(addCluster[0], halfAddClusters);
    				} 
        		} else if (subtractApprox){
        			totalCluster = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTarget}, new MayerFunctionNonAdditive[]{e3Target});
        	    	addCluster[0] =  new ClusterSum(clustersAdd, weightsAdd, new MayerFunction[]{eTarget});  	
        			
        			ClusterAbstract[] approxTotalClusters = new ClusterAbstract[1];
        			approxTotalClusters[0] = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetApprox}, new MayerFunctionNonAdditive[]{e3TargetApprox});
	    			
	    			ClusterAbstract[] approxAddClusters = new ClusterAbstract[1];
	    			approxAddClusters[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetApprox});
	    			
	    			totalCluster = new ClusterDifference(totalCluster, approxTotalClusters);
	    			addCluster[0] = new ClusterDifference(addCluster[0], approxAddClusters);
        		} else {
        			totalCluster = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTarget}, new MayerFunctionNonAdditive[]{e3Target});
        	    	addCluster[0] =  new ClusterSum(clustersAdd, weightsAdd, new MayerFunction[]{eTarget});  	
        			if (subtractHalf) {
        				ClusterAbstract[] halfTotalClusters = new ClusterAbstract[2];
    	    			halfTotalClusters[0] = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetSkip[0]}, new MayerFunctionNonAdditive[]{e3TargetSkip[0]});
    	    			halfTotalClusters[1] = new ClusterSumMultibody(clusters, weights, new MayerFunction[]{fTargetSkip[1]}, new MayerFunctionNonAdditive[]{e3TargetSkip[1]});
    	    			
    	    			ClusterAbstract[] halfAddClusters = new ClusterAbstract[2];
    	    			halfAddClusters[0] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetSkip[0]});
    	    			halfAddClusters[1] = new ClusterSum(clusters, weights, new MayerFunction[]{eTargetSkip[1]});
    	    			
    	    			totalCluster = new ClusterDifference(totalCluster, halfTotalClusters);
    	    			addCluster[0] = new ClusterDifference(addCluster[0], halfAddClusters);
        			} 
        		}
    			
    			targetCluster = new ClusterDifference(totalCluster, addCluster);
        		
        		//additive contribution as reference
        		refCluster = new ClusterSum(clustersAdd, weightsAdd, new MayerFunction[]{eTarget});
    		}
    	}

        samplingCluster = ClusterWeightAbs.makeWeightCluster(refCluster);

        // the cluster's temperature determines the factor multiplied in the exponential (f=e-1)
        // we want 1/(P*kT)
        targetCluster.setTemperature(temperature);
        refCluster.setTemperature(temperature);
        
        System.out.println(steps+" steps");
        double heMass = 4.002602;
        double lambda = Constants.PLANCK_H/Math.sqrt(2*Math.PI*heMass*temperature);
        SpeciesSpheres species = new SpeciesSpheres(space, nBeads, new AtomType(new ElementChemical("He", heMass, 2)), new ConformationLinear(space, 0));
        // the temperature here goes to the integrator, which uses it for the purpose of intramolecular interactions
        // we handle that manually below, so just set T=1 here
        
        final SimulationVirial sim = new SimulationVirial(space, species, temperature, samplingCluster, refCluster, new ClusterAbstract[]{targetCluster});

        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveTranslate);
        sim.integrator.getMoveManager().removeMCMove(sim.mcMoveRotate);

        int[][] atoms;
        if (nPoints == 2){
        	atoms = new int[][]{{0,1}};
        } else if (nPoints == 3) {
        	atoms = new int[][]{{0,1,2}};
        } else if (nPoints == 4) {
        	atoms = new int[][]{{0,1,2,3}};
        } else {
        	throw new RuntimeException("Cannot consider that order yet.");
        }
        MCMoveBox ring;
        double energyFac = nBeads*Math.PI/(lambda*lambda);
        if (aRef < 1) {
            MCMoveClusterRingPartialRegrow ringMove;
            ringMove = new MCMoveClusterRingPartialRegrow(sim.integrator.getPotentialMaster(), sim.getRandom(), space, atoms);
            ringMove.setEnergyFactor(energyFac);
            int numRegrowBeads = nBeads/3;
            System.out.println("regrow "+numRegrowBeads+" beads");
            ringMove.setNumBeads(numRegrowBeads);
            ring = ringMove;
        }
        else {
            MCMoveClusterRingRegrow ringMove; 
            ringMove = new MCMoveClusterRingRegrow(sim.getRandom(), space, atoms);
            ringMove.setEnergyFactor(energyFac);
            int numCBMCTrials = nBeads/40+10;
            System.out.println("regrow full ring");
            ringMove.setNumTrial(numCBMCTrials);
            ring = ringMove;
        }
        sim.integrator.getMoveManager().addMCMove(ring);
        
        MCMoveClusterRingScale ringScale;
        ringScale = new MCMoveClusterRingScale(sim.integrator.getPotentialMaster(),sim.getRandom(), space, atoms);
        ringScale.setEnergyFactor(energyFac);
        sim.integrator.getMoveManager().addMCMove(ringScale);
        sim.integrator.getMoveManager().setFrequency(ringScale, 0.01);
//        AtomActionTranslateBy translator = new AtomActionTranslateBy(space);
//        IVectorRandom groupTranslationVector = (IVectorRandom)translator.getTranslationVector();
//        MoleculeChildAtomAction moveMoleculeAction = new MoleculeChildAtomAction(translator);
//        IMoleculeList molecules = sim.box.getMoleculeList();
//        for (int i=1; i<nPoints; i++) {
//            groupTranslationVector.setX(0, i*3);
//            moveMoleculeAction.actionPerformed(molecules.getMolecule(i));
//            molecules.getMolecule(i).getChildList().getAtom(i).getPosition().setX(1, 1);
//        }
//        sim.box.trialNotify();
//        double pi = sim.box.getSampleCluster().value(sim.box);
//        if (pi == 0) throw new RuntimeException("initialization failed");
//        sim.box.acceptNotify();


        AtomType type = species.getLeafType();
        
        
        
        double r = 2;
        IAtomList leafList = sim.box.getLeafList();
        for (int i = 0; i<leafList.size(); i++) {
            Vector p = leafList.get(i).getPosition();
            p.setX(0, r*Math.cos((2*Math.PI*i)/leafList.size()));
            p.setX(1, r*Math.sin((2*Math.PI*i)/leafList.size()));
        }
		
        if (false) {
            double vSize = 10;
            sim.box.getBoundary().setBoxSize(Vector.of(new double[]{vSize, vSize, vSize}));
            SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
            DisplayBox displayBox = simGraphic.getDisplayBox(sim.box); 
            displayBox.setPixelUnit(new Pixel(300.0/vSize));
            displayBox.setShowBoundary(false);
            ((DisplayBoxCanvasG3DSys)displayBox.canvas).setBackgroundColor(Color.WHITE);
            
//            IAtomList leafList = sim.box.getLeafList();
            AtomPair pair = new AtomPair();
            for (int i = 0; i<leafList.size()-1; i++) {
                pair.atom0 = leafList.get(i);
                pair.atom1 = leafList.get(i+1);
                ((DisplayBoxCanvasG3DSys)displayBox.canvas).makeBond(pair, null);
            }
//            pair.atom0 = leafList.getAtom(leafList.getAtomCount()-1);
//            pair.atom1 = leafList.getAtom(0);
//            ((DisplayBoxCanvasG3DSys)displayBox.canvas).makeBond(pair, null);
            
            DiameterHashByType diameterManager = (DiameterHashByType)displayBox.getDiameterHash();
            diameterManager.setDiameter(type, 1.0/nBeads);
            ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
            displayBox.setColorScheme(colorScheme);
            simGraphic.makeAndDisplayFrame();

            
            // if running interactively, set filename to null so that it doens't read
            // (or write) to a refpref file
            
            final DisplayTextBox averageBox = new DisplayTextBox();
            averageBox.setLabel("Average");
            final DisplayTextBox errorBox = new DisplayTextBox();
            errorBox.setLabel("Error");
            JLabel jLabelPanelParentGroup = new JLabel("ratio");
            final JPanel panelParentGroup = new JPanel(new BorderLayout());
            panelParentGroup.add(jLabelPanelParentGroup,CompassDirection.NORTH.toString());
            panelParentGroup.add(averageBox.graphic(), BorderLayout.WEST);
            panelParentGroup.add(errorBox.graphic(), BorderLayout.EAST);
            simGraphic.getPanel().controlPanel.add(panelParentGroup, SimulationPanel.getVertGBC());


            IAction pushAnswer = new IAction() {
                DataDouble data = new DataDouble();
                
                public void actionPerformed() {
                    DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
                    data.x = ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO.index)).getData()[1];
                    averageBox.putData(data);
                    data.x = ((DataDoubleArray) allYourBase.getData(AccumulatorRatioAverageCovariance.RATIO_ERROR.index)).getData()[1];
                    errorBox.putData(data);
                }
            };
            IDataInfo dataInfo = new DataDouble.DataInfoDouble("B"+nPoints, new CompoundDimension(new Dimension[]{new DimensionRatio(Volume.DIMENSION, Quantity.DIMENSION)}, new double[]{nPoints-1}));
            averageBox.putDataInfo(dataInfo);
            averageBox.setLabel("average");
            errorBox.putDataInfo(dataInfo);
            errorBox.setLabel("error");
            errorBox.setPrecision(2);
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pushAnswer));

            sim.addEquilibration(steps / 100);
            sim.getController().addActivity(new ActivityIntegrate(sim.integrator));

            return;
        }
        
        long t1 = System.currentTimeMillis();
        sim.equilibrate(steps/100);

        sim.setAccumulatorBlockSize(steps > 1000 ? steps/1000 : 1);
        
        System.out.println("equilibration finished");

//        if (false) {
//            final double refIntegralF = refIntegral;
//            IntegratorListener progressReport = new IntegratorListener() {
//                public void integratorInitialized(IIntegratorEvent e) {}
//                public void integratorStepStarted(IIntegratorEvent e) {}
//                public void integratorStepFinished(IIntegratorEvent e) {
//                    if ((sim.integrator.getStepCount()*10) % sim.getController2().getMaxSteps() != 0) return;
//                    System.out.print(sim.integrator.getStepCount()+" steps: ");
//                    double[] ratioAndError = sim.dsvo.getOverlapAverageAndError();
//                    double ratio = ratioAndError[0];
//                    double error = ratioAndError[1];
//                    System.out.println("abs average: "+ratio*refIntegralF+", error: "+error*refIntegralF);
//                }
//            };
//            sim.integratorOS.getEventManager().addListener(progressReport);
//        }

        sim.integrator.getMoveManager().setEquilibrating(false);
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), steps);
        long t2 = System.currentTimeMillis();

        System.out.println("Ring acceptance "+ring.getTracker().acceptanceRatio());
        if (aRef == 1 && nPoints==2) {
            double refIntegral = Math.pow(lambda, 3.0) * Math.pow(2.0, -2.5);
            
            System.out.println("reference integral "+refIntegral);
        }

        DataGroup allYourBase = (DataGroup)sim.accumulator.getData();
        
        System.out.println();
        System.out.println("reference average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[0]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[0]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.ERROR.index)).getData()[0]);

        double ratio = ((DataDoubleArray) allYourBase.getData(sim.accumulator.RATIO.index)).getData()[1];
        double error = ((DataDoubleArray) allYourBase.getData(sim.accumulator.RATIO_ERROR.index)).getData()[1];

        System.out.println("target average: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.AVERAGE.index)).getData()[1]
                + " stdev: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.STANDARD_DEVIATION.index)).getData()[1]
                + " error: " + ((DataDoubleArray) allYourBase.getData(sim.accumulator.ERROR.index)).getData()[1]);

        System.out.println();
        System.out.println("ratio average: "+ratio+", error: "+error);
        
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
        public int nPoints = 4;
        public int nBeads = 16; //-1;
        public double temperature = 6;   // Kelvin
        public long numSteps = 1000000;
        public double aRef = 1;       // a=1 => ideal gas,  a=0 => ebond
        public double aTarget = 0;
        public boolean subtractHalf = false;
        public boolean calcApprox = false;
        public boolean subtractApprox = true;
        public boolean pairOnly = false;
        public boolean doTotal = false;
        
    }
}
