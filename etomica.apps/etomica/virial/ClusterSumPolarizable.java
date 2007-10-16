package etomica.virial;

import etomica.atom.AtomArrayList;
import etomica.potential.Potential;
import etomica.potential.PotentialPolarizable;


public class ClusterSumPolarizable implements ClusterAbstract, java.io.Serializable {

    /**
     * Constructor for ClusterSum.
     */
    public ClusterSumPolarizable(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray) {
        if (subClusterWeights.length != subClusters.length) throw new IllegalArgumentException("number of clusters and weights must be the same");
        clusters = new ClusterBonds[subClusters.length];
        clusterWeights = subClusterWeights;
        int pointCount = subClusters[0].pointCount();
        for(int i=0; i<clusters.length; i++) {
            clusters[i] = subClusters[i];
            if(clusters[i].pointCount() != pointCount) throw new IllegalArgumentException("Attempt to construct ClusterSum with clusters having differing numbers of points");
        }
        f = fArray;
        fValues = new double[pointCount][pointCount][fArray.length];
        uijPol = new double[pointCount][pointCount];

        scfAtoms = new AtomArrayList(); // USE THIS LIST FOR ALL ATOMS, WHETHER 3 OR 4; KMB, 8/16/06
    }

    // equal point count enforced in constructor 
    public int pointCount() {
        return clusters[0].pointCount();
    }
    
    public ClusterAbstract makeCopy() {
        ClusterSumPolarizable copy = new ClusterSumPolarizable(clusters,clusterWeights,f);
        copy.setTemperature(1/beta);
        copy.setDeltaDCut(Math.sqrt(deltaDCut2));
        return copy;
    }

    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        int thisCPairID = cPairs.getID();
//        System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) return value;
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
            return value;
        }
        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;
        
        int nPoints = pointCount();
        
        // recalculate all f values for all pairs
        PotentialPolarizable scfPotential = null;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    fValues[i][j][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    fValues[j][i][k] = fValues[i][j][k];
                    scfPotential = (PotentialPolarizable) f[k].getPotential();
                    if (k==0) uijPol[i][j] = scfPotential.getLastPolarizationEnergy();
                }
            }
        }

        value = 0.0;

        for(int i=0; i<clusters.length; i++) {
            //System.out.println("clusters.length = " + clusters.length);
            // clusters.length = 1 for B3
            double v = clusters[i].value(fValues);
            value += clusterWeights[i] * v;
            //System.out.println("value = " + value);
            //System.out.println("clusterWeights["+i+"] = " + clusterWeights[i]);
        }
        
        if (nPoints == 3) {
            // check that no pair of molecules is overlapped (overlap => gij=0)
            // if any pair is overlapped, then deltaC=0
            double f12 = fValues[0][1][0];
            double f13 = fValues[0][2][0];
            double f23 = fValues[1][2][0];
            double g12 = f12+1; //Math.exp(-beta*u12);
            double g13 = f13+1; //Math.exp(-beta*u13);
            double g23 = f23+1; //Math.exp(-beta*u23);
    		if (g12*g13*g23 != 0) {
    		    // Get a handle on the list of atoms from the AtomPairSet
    	        scfAtoms.clear();
                scfAtoms.add(box.molecule(0));
                scfAtoms.add(box.molecule(1));
    	        scfAtoms.add(box.molecule(2));
    			double u123Pol = scfPotential.getPolarizationEnergy(scfAtoms);

    			//deltaC = Math.exp(-beta*u123) - Math.exp(-beta*(u12 + u13 + u23));
    			double deltau123 = u123Pol-(uijPol[0][1] + uijPol[0][2] + uijPol[1][2]);
    	        double betaU123 = beta*deltau123;
    	        double expBetaU123;
    	        if (Math.abs(betaU123) < 1.e-8) {
    	            // for small x, exp(-x)-1 ~= -x
    	            // for x < 1E-8, the approximation is value within machine precision
    	            // for x < 1E-15, exp(-x) is 1, so the approximation is more accurate
    	            //   than simply doing the math.
    	            expBetaU123 = -betaU123;
    	        }
    	        else {
    	            expBetaU123 = Math.exp(-beta*deltau123)-1;
    	        }

                double deltaC = expBetaU123*g12*g13*g23;
    			
    			// deltaC has to be multiplied by clusterWeights, just like v was multiplied by
    			// clusterWeights above to get value
    			deltaC = deltaC*clusterWeights[0];

    			if (Double.isInfinite(deltaC)) {
    				System.out.println("deltaC = " + deltaC);
    			}
                value += deltaC;
                //System.out.println("u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23 + ", u123 = " + u123 + ", deltaC = " + deltaC + ", value(before) = " + value);
    		}
		}

        else if (nPoints == 4) {
            // deltaD runs into precision problems for long distances
            boolean truncateLong = false;
            for (int i=0; i<nPoints-1; i++) {
                for (int j=i+1; j<nPoints; j++) {
                    if (cPairs.getr2(i,j) > deltaDCut2) {
                        truncateLong = true;
                        break;
                    }
                }
            }

            if (!truncateLong) {
                double f12 = fValues[0][1][0];
                double f13 = fValues[0][2][0];
                double f14 = fValues[0][3][0];
                double f23 = fValues[1][2][0];
                double f24 = fValues[1][3][0];
                double f34 = fValues[2][3][0];
                double g12 = f12+1; //Math.exp(-beta*u12);
                double g13 = f13+1; //Math.exp(-beta*u13);
                double g14 = f14+1; //Math.exp(-beta*u14);
                double g23 = f23+1; //Math.exp(-beta*u23);
                double g24 = f24+1; //Math.exp(-beta*u24);
                double g34 = f34+1; //Math.exp(-beta*u34);

                scfAtoms.clear();
                // we need to properly construct these lists even if we don't use them
                // (due to overlaps) because the next list is obtained by removing/adding
                // atoms from this one.
                scfAtoms.add(box.molecule(0));
                scfAtoms.add(box.molecule(1));
                scfAtoms.add(box.molecule(2));

                double deltaD = 0;
                // if 12 13 or 23 is overlapped, then we can't calculate u123Pol and
                // couldn't calculate the uijPol.  Fortunately, gij is 0, so the 123
                // term is 0.
                if (g12*g13*g23 != 0) {
    				double u123Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    double deltaU123 = u123Pol - (uijPol[0][1] + uijPol[0][2] + uijPol[1][2]);
                    double beta123 = beta*deltaU123;
                    // for small x, exp(-x)-1 ~= -x
                    // for x < 1E-8, the approximation is value within machine precision
                    // for x < 1E-15, exp(-x) is 1, so the approximation is more accurate
                    //   than simply doing the math.
                    double exp123 = -beta123;
                    if (Math.abs(beta123) > 1E-8) {
                        exp123 = Math.exp(-beta123) - 1;
                    }
                    // original formula has g14+g24+g34-2 instead of f14+f24+f34+1.
                    // Using f should have better precision since
                    // the g's will all be close to 1 (or even equal to 1) for 
                    // systems with molecules having large separations.
                    deltaD += -exp123*g12*g13*g23*((f14+f24+f34)+1);
                }

    			scfAtoms.remove(2);
    			scfAtoms.add(box.molecule(3));
    			if (g12*g14*g24 != 0) {
    				double u124Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    double deltaU124 = u124Pol-(uijPol[0][1]+uijPol[0][3]+uijPol[1][3]);
                    double beta124 = beta*deltaU124;
                    double exp124 = -beta124;
                    if (Math.abs(beta124) > 1E-8) {
                        exp124 = Math.exp(-beta124) - 1;
                    }
                    deltaD +=  -exp124*g12*g14*g24*((f13+f23+f34)+1);
    			}

                scfAtoms.remove(1);
                scfAtoms.add(box.molecule(2));
                if (g13*g14*g34 != 0) {
                    double u134Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    double deltaU134 = u134Pol-(uijPol[0][2]+uijPol[0][3]+uijPol[2][3]);
                    double beta134 = beta*deltaU134;
                    double exp134 = -beta134;
                    if (Math.abs(beta134) > 1E-8) {
                        exp134 = Math.exp(-beta134) - 1;
                    }
                    deltaD += -exp134*g13*g14*g34*((f12+f23+f24)+1);
                }

                scfAtoms.remove(0);
                scfAtoms.add(box.molecule(1));
                if (g23*g24*g34 != 0) {
                    double u234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    double deltaU234 = u234Pol-(uijPol[1][2]+uijPol[1][3]+uijPol[2][3]);
                    double beta234 = beta*deltaU234;
                    double exp234 = -beta234;
                    if (Math.abs(beta234) > 1E-8) {
                        exp234 = Math.exp(-beta234) - 1;
                    }
                    deltaD += -exp234*g23*g24*g34*((f12+f13+f14)+1);
                }

                scfAtoms.add(box.molecule(0));
                if (g12*g13*g14*g23*g24*g34 != 0) {
                    double u1234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    // deltaU1234 would have deltaUabc subtracted off, but we'd also add it back
                    // in for expU1234, so just don't subtract in the first place 
    				double deltaU1234 = u1234Pol-(uijPol[0][1]+uijPol[0][2]+uijPol[0][3]+uijPol[1][2]+uijPol[1][3]+uijPol[2][3]); //-deltaU123-deltaU124-deltaU134-deltaU234;
                    double beta1234 = beta*deltaU1234; //deltaU123+deltaU124+deltaU134+deltaU234+deltaU1234);
                    double exp1234 = -beta1234;
                    if (Math.abs(beta1234) > 1E-8) {
                        exp1234 = Math.exp(-beta1234) - 1;
                    }
    				deltaD += exp1234*g12*g13*g14*g23*g24*g34;
                }

                // Mason and Spurling book deltaD; 5/11/07

                // deltaU1234 = u1234-u12-u13-u14-u23-u24-u34;

//          deltaD = 2*Math.exp(-beta*(u12+u13+u23))*(Math.exp(-beta*deltaU123)-1) + 2*Math.exp(-beta*(u12+u14+u24))*(Math.exp(-beta*deltaU124)-1)
//					+ 2*Math.exp(-beta*(u13+u14+u34))*(Math.exp(-beta*deltaU134)-1) + 2*Math.exp(-beta*(u23+u24+u34))*(Math.exp(-beta*deltaU234)-1)
//					+ Math.exp(-beta*(u12+u13+u23-u14))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u13+u23-u24))*(1-Math.exp(-beta*deltaU123))
//					+ Math.exp(-beta*(u12+u13+u23-u34))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u14+u24-u13))*(1-Math.exp(-beta*deltaU124))
//					+ Math.exp(-beta*(u12+u14+u24-u23))*(1-Math.exp(-beta*deltaU124)) + Math.exp(-beta*(u12+u14+u24-u34))*(1-Math.exp(-beta*deltaU124))
//					+ Math.exp(-beta*(u13+u14+u34-u12))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u13+u14+u34-u23))*(1-Math.exp(-beta*deltaU134))
//					+ Math.exp(-beta*(u13+u14+u34-u24))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u23+u24+u34-u12))*(1-Math.exp(-beta*deltaU234))
//					+ Math.exp(-beta*(u23+u24+u34-u13))*(1-Math.exp(-beta*deltaU234)) + Math.exp(-beta*(u23+u24+u34-u14))*(1-Math.exp(-beta*deltaU234))
//					+ Math.exp(-beta*(u12+u13+u14+u23+u24+u34))*(Math.exp(-beta*deltaU1234)-1);


//        				 kmb added this line 8/16/06
                // deltaD has to be multiplied by weightPrefactor from Standard class, just like deltaC was multiplied by
                // clusterWeights above to get value; note, for B3 clusterWeights = weightPrefactor

                // -(1/8) is the B4 prefactor multiplying all diagrams.
                deltaD = -0.125*deltaD;
                value += deltaD;
            }
        }
        else if (nPoints == 5) {
            double f12 = fValues[0][1][0];
            double f13 = fValues[0][2][0];
            double f14 = fValues[0][3][0];
            double f15 = fValues[0][4][0];
            double f23 = fValues[1][2][0];
            double f24 = fValues[1][3][0];
            double f25 = fValues[1][4][0];
            double f34 = fValues[2][3][0];
            double f35 = fValues[2][4][0];
            double f45 = fValues[3][4][0];
            double g12 = f12+1; //Math.exp(-beta*u12);
            double g13 = f13+1; //Math.exp(-beta*u13);
            double g14 = f14+1; //Math.exp(-beta*u14);
            double g15 = f15+1; //Math.exp(-beta*u14);
            double g23 = f23+1; //Math.exp(-beta*u23);
            double g24 = f24+1; //Math.exp(-beta*u24);
            double g25 = f25+1; //Math.exp(-beta*u14);
            double g34 = f34+1; //Math.exp(-beta*u34);
            double g35 = f35+1; //Math.exp(-beta*u14);
            double g45 = f45+1; //Math.exp(-beta*u14);
                
            // can't do this!  B5 terms need to be separated out and each used
            // if it corresponds to an set of atoms with no overlaps
                if (g12*g13*g14*g15*g23*g24*g25*g34*g35*g45 != 0) {
                    for(int i=0; i<nPoints-1; i++) {
                        for(int j=i+1; j<nPoints; j++) {
                            uijPol[i][j] = ((Potential)scfPotential).energy(aPairs.getAPair(i,j));
                        }
                    }

                    
                    //System.out.println("Sum of pair energies is infinity: uijPol[0][1] = " + uijPol[0][1] + ", uijPol[0][2] = " + uijPol[0][2] + ", uijPol[1][2] = " + uijPol[1][2]);
                    scfAtoms.clear();
                    // we need to properly construct these lists even if we don't use them
                    // (due to overlaps) because the next list is obtained by removing/adding
                    // atoms from this one.
                    scfAtoms.add(box.molecule(0));
                    scfAtoms.add(box.molecule(1));
                    scfAtoms.add(box.molecule(2));  // 123
                    //System.out.println("Hey, in the heart of the delta B5 code");
                    double u123 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(2);
                    scfAtoms.add(box.molecule(3));  // 124
                    double u124 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(2);
                    scfAtoms.add(box.molecule(4));  // 125
                    double u125 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(2));  // 153
                    double u135 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(3));  // 134
                    double u134 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(4));  // 145
                    double u145 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(0);
                    scfAtoms.add(box.molecule(1));  // 452
                    double u245 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(2));  // 423
                    double u234 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(0);
                    scfAtoms.add(box.molecule(4));  // 235
                    double u235 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(0);
                    scfAtoms.add(box.molecule(3));  // 354
                    double u345 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.add(box.molecule(0));  // 3541
                    double u1345 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(1));  // 3412
                    double u1234 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(1);
                    scfAtoms.add(box.molecule(4));  // 3125
                    double u1235 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(0);
                    scfAtoms.add(box.molecule(3));  // 1254
                    double u1245 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.remove(0);
                    scfAtoms.add(box.molecule(2));  // 2543
                    double u2345 = ((Potential)scfPotential).energy(scfAtoms);
                    scfAtoms.add(box.molecule(0));  // 25431
                    double u12345 = ((Potential)scfPotential).energy(scfAtoms);
                    
                    double deltaU123 = u123-uijPol[0][1]-uijPol[0][2]-uijPol[1][2];
                    double deltaU124 = u124-uijPol[0][1]-uijPol[0][3]-uijPol[1][3];
                    double deltaU125 = u125-uijPol[0][1]-uijPol[0][4]-uijPol[1][4];
                    double deltaU134 = u134-uijPol[0][2]-uijPol[0][3]-uijPol[2][3];
                    double deltaU135 = u135-uijPol[0][2]-uijPol[0][4]-uijPol[2][4];
                    double deltaU145 = u145-uijPol[0][4]-uijPol[0][3]-uijPol[3][4];
                    double deltaU234 = u234-uijPol[1][2]-uijPol[1][3]-uijPol[2][3];
                    double deltaU235 = u235-uijPol[1][2]-uijPol[1][4]-uijPol[2][4];
                    double deltaU245 = u245-uijPol[1][4]-uijPol[1][3]-uijPol[3][4];
                    double deltaU345 = u345-uijPol[2][3]-uijPol[2][4]-uijPol[3][4];
                    
                    double deltaU1234 = u1234-uijPol[0][1]-uijPol[0][2]-uijPol[0][3]-uijPol[1][2]-uijPol[1][3]-uijPol[2][3]-deltaU123-deltaU124-deltaU134-deltaU234;
                    double deltaU1235 = u1235-uijPol[0][1]-uijPol[0][2]-uijPol[0][4]-uijPol[1][2]-uijPol[1][4]-uijPol[2][4]-deltaU123-deltaU125-deltaU135-deltaU235;
                    double deltaU1245 = u1245-uijPol[0][1]-uijPol[0][3]-uijPol[0][4]-uijPol[1][3]-uijPol[1][4]-uijPol[3][4]-deltaU124-deltaU125-deltaU145-deltaU245;
                    double deltaU1345 = u1345-uijPol[0][2]-uijPol[0][3]-uijPol[0][4]-uijPol[2][3]-uijPol[2][4]-uijPol[3][4]-deltaU134-deltaU135-deltaU145-deltaU345;
                    double deltaU2345 = u2345-uijPol[1][2]-uijPol[1][3]-uijPol[1][4]-uijPol[2][3]-uijPol[2][4]-uijPol[3][4]-deltaU234-deltaU235-deltaU245-deltaU345;
                    
                    double deltaU12345 = u12345-uijPol[0][1]-uijPol[0][2]-uijPol[0][3]-uijPol[0][4]-uijPol[1][2]-uijPol[1][3]-uijPol[1][4]-uijPol[2][3]-uijPol[2][4]-uijPol[3][4]-
                    deltaU123-deltaU124-deltaU125-deltaU134-deltaU135-deltaU145-deltaU234-deltaU235-deltaU245-deltaU345-deltaU1234-deltaU1235-deltaU1245-deltaU1345-deltaU2345;
                    
//                  System.out.println(uijPol[0][1] + " " + uijPol[0][2] + " " + uijPol[0][3] + " " + uijPol[0][4] + " " + uijPol[1][2] + " " + uijPol[1][3] + " " + uijPol[1][4] + " " + uijPol[2][3] + " " + uijPol[2][4] + " " + uijPol[3][4] + " " + u123 + " " + u124 + " " + u125 + " " + u134 + " " + u135 + " " + u145 + " " + u234 + " " + u235 + " " + u245 + " " + u345 + " " + u1234 + " " + u1235 + " " + u1245 + " " + u1345 + " " + u2345 + " " + u12345);
//                  System.out.println(deltaU123 + " " + deltaU124 + " "  + deltaU125 + " " + deltaU134 + " "  + deltaU135 + " " + deltaU145 + " " + deltaU234 + " " + deltaU235 + " " + deltaU245 + " " + deltaU345 + " " + deltaU1234 + " " + deltaU1235 + " " + deltaU1245 + " " + deltaU1345 + " " + deltaU2345 + " " + deltaU12345);
                    
/*                      double g12 = Math.exp(-beta*uijPol[0][1]);
                    double g13 = Math.exp(-beta*uijPol[0][2]);
                    double g14 = Math.exp(-beta*uijPol[0][3]);
                    double g23 = Math.exp(-beta*uijPol[1][2]);
                    double g24 = Math.exp(-beta*uijPol[1][3]);
                    double g34 = Math.exp(-beta*uijPol[2][3]);
*/                      

                    // DON'T FORGET THE BETAS! KMB 5/7/07
                    double betaU12 = uijPol[0][1]*beta;
                    double betaU13 = uijPol[0][2]*beta;
                    double betaU14 = uijPol[0][3]*beta;
                    double betaU15 = uijPol[0][4]*beta;
                    double betaU23 = uijPol[1][2]*beta;
                    double betaU24 = uijPol[1][3]*beta;
                    double betaU25 = uijPol[1][4]*beta;
                    double betaU34 = uijPol[2][3]*beta;
                    double betaU35 = uijPol[2][4]*beta;
                    double betaU45 = uijPol[3][4]*beta;

                    betaU12 = -Math.log(g12);
                    betaU13 = -Math.log(g13);
                    betaU14 = -Math.log(g14);
                    betaU15 = -Math.log(g15);
                    betaU23 = -Math.log(g23);
                    betaU24 = -Math.log(g24);
                    betaU25 = -Math.log(g25);
                    betaU34 = -Math.log(g34);
                    betaU35 = -Math.log(g35);
                    betaU45 = -Math.log(g45);

                    double du012 = deltaU123*beta;
                    double du013 = deltaU124*beta;
                    double du014 = deltaU125*beta;
                    double du023 = deltaU134*beta;
                    double du024 = deltaU135*beta;
                    double du034 = deltaU145*beta;
                    double du123 = deltaU234*beta;
                    double du124 = deltaU235*beta;
                    double du134 = deltaU245*beta;
                    double du234 = deltaU345*beta;
                    
                    double du0123 = deltaU1234*beta;
                    double du0124 = deltaU1235*beta;
                    double du0134 = deltaU1245*beta;
                    double du0234 = deltaU1345*beta;
                    double du1234 = deltaU2345*beta;
                    
                    double du01234 = deltaU12345*beta;
                    
                    double deltaE = ((-6)*Math.exp((-betaU12) - betaU13 - betaU23) + 6 
                            *Math.exp((-du012) - betaU12 - betaU13 - betaU23) + 
                            3*Math.exp((-betaU12) - betaU13 - betaU14 - betaU23) - 3*Math.exp
                        ((-du012) - betaU12 - betaU13 - betaU14 - 
                            betaU23) + 3*Math.exp((-betaU12) - betaU13 - betaU15 - betaU23) - 3 
                        *Math.exp((-du012) - betaU12 - betaU13 -
                             betaU15 - betaU23) - 2*Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU23) + 
                        2*Math.exp((-du012) - betaU12 - betaU13 - betaU14 -
                             betaU15 - betaU23) - 6*Math.exp((-betaU12) - betaU14 - betaU24) + 6 
                        *Math.exp((-du013) - betaU12 -
                             betaU14 - betaU24) + 3*Math.exp((-betaU12) - betaU13 - betaU14 - betaU24) - 3 
                        *Math.exp((-du013) - betaU12 - betaU13 - betaU14 - 
                            betaU24) + 3*Math.exp((-betaU12) - betaU14 - betaU15 - betaU24) - 3 
                        *Math.exp((-du013) - betaU12 - betaU14 - betaU15 - betaU24) - 
                            2*Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU24) + 2 
                        *Math.exp((-du013) - betaU12 - betaU13 - betaU14 - betaU15 - betaU24) + 3 
                        *Math.exp((-betaU12) - betaU13 - betaU23 - betaU24) - 3 
                        *Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU24) + 3 
                        *Math.exp((-betaU12) - betaU14 - betaU23 - betaU24) - 3 
                        *Math.exp((-du013) - betaU12 - betaU14 - betaU23 - betaU24) - Math.exp(
                        (-betaU12) - betaU13 - betaU15 - betaU23 - betaU24) + Math.exp((-du012) - betaU12 - 
                        betaU13 - betaU15 - betaU23 - betaU24) - Math.exp((-betaU12) - betaU14 - betaU15 - betaU23 - betaU24
                        ) + Math.exp((-du013) - betaU12 - betaU14 - betaU15 - betaU23 - 
                            betaU24) - 6*Math.exp((-betaU12) - betaU15 - betaU25) + 6*Math.exp
                        ((-du014) - betaU12 - betaU15 - betaU25) + 3*Math.exp((-betaU12) - betaU13 - 
                        betaU15 - betaU25) - 3*Math.exp((-du014) - betaU12 - betaU13 - betaU15 - betaU25) + 3 
                        *Math.exp((-betaU12) - betaU14 - betaU15 - betaU25) - 3 
                        *Math.exp((-du014) - betaU12 - betaU14 - betaU15 - betaU25) - 2 
                        *Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU25) + 2 
                        *Math.exp((-du014) - betaU12 - betaU13 - betaU14 - betaU15 - betaU25) + 3 
                        *Math.exp((-betaU12) - betaU13 - betaU23 - 
                            betaU25) - 3*Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU25) - 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU23 - betaU25) + 
                        Math.exp((-du012) - betaU12 - betaU13 - betaU14 - betaU23 - betaU25) + 3 
                        *Math.exp((-betaU12) - betaU15 - betaU23 - betaU25) - 3 
                        *Math.exp((-du014) - betaU12 - betaU15 - betaU23 - betaU25) - Math.exp(
                        (-betaU12) - betaU14 - betaU15 - betaU23 - betaU25) + Math.exp((-du014) - betaU12 - 
                        betaU14 - betaU15 - betaU23 - betaU25) + 3*Math.exp((-
                            betaU12) - betaU14 - betaU24 - betaU25) - 3*Math.exp((-du013) - betaU12 - betaU14 
                        - betaU24 - betaU25) - Math.exp((-
                            betaU12) - betaU13 - betaU14 - betaU24 - betaU25) + Math.exp((-du013) - betaU12 - 
                        betaU13 - betaU14 - betaU24 - betaU25) + 3*Math.exp((-betaU12) - betaU15 - betaU24 - betaU25) 
                        - 3*Math.exp((-du014) - betaU12 - betaU15 - betaU24 - betaU25) - 
                        Math.exp((-
                            betaU12) - betaU13 - betaU15 - betaU24 - betaU25) + Math.exp((-du014) - betaU12 - 
                        betaU13 - betaU15 - betaU24 - betaU25) - 2*Math.exp((-betaU12) - betaU13 - betaU23 - betaU24 - 
                        betaU25) + 2*Math.exp((-du012) - betaU12 - 
                            betaU13 - betaU23 - betaU24 - betaU25) - 2*Math.exp((-betaU12) - betaU14 - betaU23 - 
                        betaU24 - betaU25) + 2*Math.exp((-du013) -
                             betaU12 - betaU14 - betaU23 - betaU24 - betaU25) - 2*Math.exp((-betaU12) - betaU15 - 
                        betaU23 - betaU24 - betaU25) + 2*Math.exp((-du014) - betaU12 - betaU15 - betaU23 - betaU24 
                        - betaU25) - 6*Math.exp((-betaU13) - betaU14 - betaU34) + 6*Math.exp
                        ((-du023) - betaU13 -
                             betaU14 - betaU34) + 3*Math.exp((-betaU12) - betaU13 - betaU14 - betaU34) - 3 
                        *Math.exp((-du023) - betaU12 - betaU13 - betaU14 - betaU34) + 
                        3*Math.exp((-betaU13) - betaU14 - betaU15 - betaU34) - 3 
                        *Math.exp((-du023) - betaU13 - betaU14 - betaU15 - betaU34) - 2 
                        *Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU34) + 2 
                        *Math.exp((-du023) - betaU12 - betaU13 - betaU14 - betaU15 - betaU34) + 3 
                        *Math.exp((-betaU12) - betaU13 - betaU23 - betaU34) - 3 
                        *Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU34) + 3 
                        *Math.exp((-betaU13) - betaU14 - betaU23 - betaU34) - 3 
                        *Math.exp((-du023) - betaU13 - betaU14 - betaU23 - betaU34) - Math.exp(
                        (-betaU12) - betaU13 - betaU15 - betaU23 - betaU34) + Math.exp((-du012) - betaU12 - 
                        betaU13 - betaU15 - betaU23 - betaU34) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU23 - betaU34
                        ) + Math.exp((-du023) - betaU13 - betaU14 - betaU15 - betaU23 - 
                            betaU34) + 3*Math.exp((-betaU12) - betaU14 - betaU24 - betaU34) - 3 
                        *Math.exp((-du013) - betaU12 -
                             betaU14 - betaU24 - betaU34) + 3*Math.exp((-betaU13) - betaU14 - betaU24 - betaU34) - 
                        3*Math.exp((-du023) - betaU13 - betaU14 - 
                            betaU24 - betaU34) - Math.exp((-
                        betaU12) - betaU14 - betaU15 - betaU24 - betaU34) + Math.exp((-du013) - betaU12 -
                             betaU14 - betaU15 - betaU24 - betaU34) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU24 
                        - betaU34) + Math.exp((-du023) - betaU13 - 
                            betaU14 - betaU15 - betaU24 - betaU34) - 6*Math.exp((-betaU23) - betaU24 - betaU34) + 
                        6*Math.exp((-du123) - betaU23 - betaU24 - betaU34) + 3*Math.exp(
                        (-betaU12) - betaU23 - betaU24 - betaU34) - 3*Math.exp((-du123) - betaU12 - betaU23 - 
                        betaU24 - betaU34) + 
                            3*Math.exp((-betaU13) - betaU23 - betaU24 - betaU34) - 3*Math.exp
                        ((-du123) - betaU13 - betaU23 - betaU24 - betaU34) + 3*Math.exp((-betaU14) - betaU23 
                        - betaU24 - betaU34) - 3*Math.exp((-du123) - betaU14 - betaU23 - betaU24 - betaU34) - 
                        3*Math.exp((-betaU12) - betaU13 - betaU14 - betaU23 - betaU24 - betaU34) + 3 
                        *Math.exp((-du012) - 
                            du0123 - du013 - du023 - du123 - betaU12 - betaU13 - betaU14 - betaU23 - betaU24 - betaU34) + 2 
                        *Math.exp((-betaU15) - betaU23 - betaU24 - 
                            betaU34) - 2*Math.exp((-du123) - betaU15 - betaU23 - betaU24 - betaU34) - 
                        Math.exp((-betaU12) - betaU15 - betaU23 - betaU24 - betaU34) + 
                        Math.exp((-du123) - betaU12 - betaU15 - betaU23 - betaU24 - betaU34) - 
                        Math.exp((-betaU13) - betaU15 - 
                            betaU23 - betaU24 - betaU34) + Math.exp((-du123) - betaU13 - betaU15 - betaU23 - betaU24 
                        - betaU34) - Math.exp((-betaU14) - betaU15 - betaU23 - betaU24 - betaU34) + 
                        Math.exp((-du123) - betaU14 - betaU15 - betaU23 - betaU24 - betaU34) + 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU23 - betaU24 - betaU34) - 
                        Math.exp((-du012) - du0123 - 
                            du013 - du023 - du123 - betaU12 - betaU13 - betaU14 - 
                            betaU15 - betaU23 - betaU24 - betaU34) + 2*Math.exp((-betaU13) - betaU14 - betaU25 - 
                        betaU34) - 2*Math.exp((-du023) - betaU13 - betaU14 -
                             betaU25 - betaU34) - Math.exp((-betaU12) - betaU13 - betaU14 - betaU25 - betaU34) + 
                        Math.exp((-du023) - betaU12 - betaU13 - 
                            betaU14 - betaU25 - betaU34) + 2*Math.exp((-betaU12) - betaU15 - betaU25 - betaU34) - 
                        2*Math.exp((-du014) - betaU12 - betaU15 - betaU25 - betaU34) - 
                        Math.exp((-betaU12) - betaU13 - betaU15 - betaU25 - betaU34) + 
                        Math.exp((-du014) - betaU12 - betaU13 - betaU15 - betaU25 - betaU34) - 
                        Math.exp((-betaU12) - betaU14 - betaU15 - 
                        betaU25 - betaU34) + Math.exp((-du014) - betaU12 - betaU14 - betaU15 - betaU25 - betaU34) - 
                        Math.exp((-betaU13) - betaU14 - betaU15 -
                             betaU25 - betaU34) + Math.exp((-du023) - betaU13 - betaU14 - betaU15 - betaU25 - 
                        betaU34) + Math.exp((-betaU12) - betaU13 - betaU14 - 
                        betaU15 - betaU25 - betaU34) - Math.exp((-du014) - du023 - betaU12 - betaU13 - betaU14 - 
                        betaU15 - betaU25 - betaU34) - Math.exp((-betaU12) - betaU13 - betaU23 - betaU25 - betaU34) + 
                        Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU25 - 
                            betaU34) - Math.exp((-betaU13) - betaU14 - betaU23 - betaU25 - betaU34) + 
                        Math.exp((-du023) - betaU13 - betaU14 - betaU23 - betaU25 - betaU34) - 
                        Math.exp((-betaU12) - betaU15 - betaU23 - betaU25 - betaU34) + 
                        Math.exp((-du014) - betaU12 - betaU15 - betaU23 - betaU25 - betaU34) - 
                        Math.exp((-betaU12) - betaU14 - betaU24 - betaU25 - betaU34) + 
                        Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU25 - betaU34) - 
                        Math.exp((-betaU13) - betaU14 - betaU24 - betaU25 - betaU34) + 
                        Math.exp((-du023) - 
                            betaU13 - betaU14 - betaU24 - betaU25 - betaU34) - Math.exp((-betaU12) - betaU15 - betaU24 - 
                        betaU25 - betaU34) + Math.exp((-du014) - betaU12 -
                             betaU15 - betaU24 - betaU25 - betaU34) + 3*Math.exp((-betaU23) - betaU24 - betaU25 - 
                        betaU34) - 3*Math.exp((-du123) - 
                            betaU23 - betaU24 - betaU25 - betaU34) - 2*Math.exp((-betaU12) - betaU23 - betaU24 - 
                        betaU25 - betaU34) + 2*Math.exp((-du123) - betaU12 -
                             betaU23 - betaU24 - betaU25 - betaU34) - Math.exp((-betaU13) - betaU23 - betaU24 - betaU25 
                        - betaU34) + Math.exp((-
                            du123) - betaU13 - betaU23 - betaU24 - betaU25 - betaU34) - Math.exp((-betaU14) - 
                        betaU23 - betaU24 - betaU25 - betaU34) + Math.exp((-du123) - betaU14 - betaU23 - betaU24 - 
                        betaU25 - betaU34) + Math.exp((-betaU12) - betaU13 - betaU14 - betaU23 - betaU24 - betaU25 - betaU34
                        ) - Math.exp((-du012) - du0123 - du013 - du023 - du123 - betaU12 - 
                        betaU13 - betaU14 - betaU23 - betaU24 - betaU25 - betaU34) - Math.exp((-betaU15) - betaU23 - 
                            betaU24 - betaU25 - betaU34) + Math.exp((-du123) - betaU15 - betaU23 - betaU24 - betaU25 
                        - betaU34) + Math.exp((-betaU12) - betaU15 - betaU23 - betaU24 - betaU25 - betaU34) - 
                        Math.exp((-du014) - 
                            du123 - betaU12 - betaU15 - betaU23 - betaU24 - betaU25 - betaU34) - 6 
                        *Math.exp((-betaU13) - betaU15 - betaU35) + 6*Math.exp((-du024) 
                        - betaU13 - betaU15 - betaU35) + 3*Math.exp((-betaU12) - betaU13 - betaU15 - betaU35) - 3 
                        *Math.exp((-du024) - betaU12 - betaU13 - betaU15 - betaU35) + 3 
                        *Math.exp((-betaU13) - betaU14 - betaU15 - betaU35) - 3 
                        *Math.exp((-du024) - betaU13 - 
                            betaU14 - betaU15 - betaU35) - 2*Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - 
                        betaU35) + 2*Math.exp((-du024) - betaU12 - betaU13 - betaU14 - betaU15 - betaU35) + 3 
                        *Math.exp((-betaU12) - betaU13 - betaU23 - betaU35) - 3 
                        *Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU35) - Math.exp(
                        (-betaU12) - betaU13 - betaU14 - betaU23 - betaU35) + Math.exp((-
                            du012) - betaU12 - betaU13 - betaU14 - betaU23 - betaU35) + 3*Math.exp((-betaU13) 
                        - betaU15 - betaU23 - betaU35) - 3*Math.exp((-du024) - betaU13 - betaU15 - betaU23 - 
                        betaU35) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU23 - betaU35) + 
                        Math.exp((-du024) - betaU13 - betaU14 - betaU15 - betaU23 - betaU35) + 2 
                        *Math.exp((-betaU12) - betaU14 - betaU24 - 
                        betaU35) - 2*Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU35) - 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU24 - betaU35) + 
                        Math.exp((-du013) - betaU12 - betaU13 - betaU14 - betaU24 - betaU35) + 
                            2*Math.exp((-betaU13) - 
                        betaU15 - betaU24 - betaU35) - 2*Math.exp((-du024) - betaU13 - betaU15 - betaU24 - 
                        betaU35) - Math.exp((-betaU12) - betaU13 - betaU15 - betaU24 - betaU35) + 
                        Math.exp((-du024) - betaU12 - betaU13 - betaU15 - betaU24 - betaU35) - 
                        Math.exp((-betaU12) - betaU14 - betaU15 - betaU24 - betaU35) + 
                        Math.exp((-du013) - betaU12 - betaU14 - betaU15 - betaU24 - 
                            betaU35) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU24 - betaU35) + 
                        Math.exp((-du024) - betaU13 - betaU14 - betaU15 - betaU24 - betaU35) + 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU24 - betaU35) - Math.exp
                        ((-du013) - du024 - betaU12 - betaU13 - betaU14 - betaU15 - betaU24 - betaU35) - 
                        Math.exp((-betaU12) - betaU13 - betaU23 - betaU24 - betaU35) + 
                        Math.exp((-du012) - betaU12 - 
                            betaU13 - betaU23 - betaU24 - betaU35) - Math.exp((-betaU12) - betaU14 - betaU23 - betaU24 - 
                        betaU35) + Math.exp((-du013) - betaU12 - betaU14 - betaU23 - betaU24 - betaU35) - 
                        Math.exp((-betaU13) - betaU15 - betaU23 - betaU24 - betaU35) + 
                        Math.exp((-du024) - betaU13 - betaU15 - betaU23 - betaU24 - betaU35) + 3 
                        *Math.exp((-betaU12) - betaU15 - betaU25 - betaU35) - 3 
                        *Math.exp((-du014) - betaU12 - betaU15 - betaU25 - betaU35) + 3 
                        *Math.exp((-betaU13) - betaU15 - betaU25 - betaU35) - 3 
                        *Math.exp((-du024) - betaU13 - betaU15 - betaU25 - betaU35) - Math.exp(
                        (-betaU12) - betaU14 - betaU15 - betaU25 - betaU35) + Math.exp((-du014) - betaU12 - 
                        betaU14 - betaU15 - betaU25 - betaU35) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU25 - betaU35
                        ) + Math.exp((-du024) - betaU13 - betaU14 - betaU15 - betaU25 - betaU35) - 6 
                        *Math.exp((-betaU23) - betaU25 - betaU35) + 6*Math.exp((-du124) 
                        - betaU23 - betaU25 - betaU35) + 3*Math.exp((-betaU12) - betaU23 - betaU25 - betaU35) - 3 
                        *Math.exp((-du124) - betaU12 - betaU23 - betaU25 - betaU35) + 3 
                        *Math.exp((-betaU13) - betaU23 - betaU25 - betaU35) - 3 
                        *Math.exp((-du124) - betaU13 - betaU23 - betaU25 - betaU35) + 2 
                        *Math.exp((-betaU14) - betaU23 - betaU25 - betaU35) - 2 
                        *Math.exp((-du124) - betaU14 - betaU23 - betaU25 - betaU35) - Math.exp(
                        (-betaU12) - betaU14 - betaU23 - betaU25 - betaU35) + Math.exp((-
                            du124) - betaU12 - betaU14 - betaU23 - betaU25 -
                         betaU35) - Math.exp((-betaU13) - betaU14 - betaU23 - betaU25 - betaU35) + 
                        Math.exp((-du124) - betaU13 - betaU14 - betaU23 - betaU25 - betaU35) + 3 
                        *Math.exp((-betaU15) - betaU23 - betaU25 - betaU35) - 3*Math.exp((-
                            du124) - betaU15 - betaU23 - betaU25 - betaU35) - 3*Math.exp((-betaU12) - betaU13 
                        - betaU15 - betaU23 - betaU25 - betaU35) + 3*Math.exp((-du012) - du0124 - 
                        du014 - du024 - du124 - betaU12 - betaU13 - betaU15 - betaU23 - betaU25 - betaU35) - Math.exp
                        ((-betaU14) - betaU15 - betaU23 - betaU25 - betaU35) + Math.exp((-du124) - betaU14 
                        - betaU15 - betaU23 - betaU25 - betaU35) + Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - 
                        betaU23 - betaU25 - betaU35) - Math.exp((-du012) -
                             du0124 - du014 - du024 - du124 - betaU12 - betaU13 - betaU14 - betaU15 - betaU23 - betaU25 - betaU35
                        ) - Math.exp((-betaU12) - betaU14 - betaU24 - betaU25 - betaU35) + 
                        Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU25 - betaU35) - 
                        Math.exp((-betaU12) -
                         betaU15 - betaU24 - betaU25 - betaU35) + Math.exp((-du014) - betaU12 - betaU15 - betaU24 - 
                        betaU25 - betaU35) - Math.exp((-betaU13) - betaU15 - betaU24 - betaU25 - betaU35) + 
                        Math.exp((-du024) - betaU13 - betaU15 - betaU24 - betaU25 - betaU35) + 
                            3*Math.exp((-betaU23) - betaU24 - betaU25 - betaU35) - 3*Math.exp
                        ((-du124) - betaU23 - betaU24 - betaU25 - betaU35) - 2*Math.exp((-betaU12) - betaU23 
                        - betaU24 - betaU25 - betaU35) + 2*Math.exp((-du124) - betaU12 - betaU23 - betaU24 - 
                        betaU25 - betaU35) - Math.exp((-betaU13) - betaU23 -
                         betaU24 - betaU25 - betaU35) + Math.exp((-du124) - betaU13 - betaU23 - betaU24 - betaU25 - 
                        betaU35) - Math.exp((-betaU14) - betaU23 - betaU24 - betaU25 - betaU35) + 
                        Math.exp((-du124) - betaU14 - betaU23 - betaU24 - betaU25 - betaU35) + 
                        Math.exp((-betaU12) - betaU14 - betaU23 - betaU24 - betaU25 - betaU35) - Math.exp
                        ((-du013) - du124 - betaU12 - betaU14 - betaU23 - 
                            betaU24 - betaU25 - betaU35) - Math.exp((-betaU15) - betaU23 - betaU24 - betaU25 - betaU35) 
                        + Math.exp((-du124) - betaU15 - betaU23 - 
                            betaU24 - betaU25 - betaU35) + Math.exp((-betaU12) - betaU13 - betaU15 - betaU23 - betaU24 - 
                        betaU25 - betaU35) - Math.exp((-du012) - du0124 - du014 - du024 - du124 
                        - betaU12 - betaU13 - betaU15 - betaU23 - betaU24 - betaU25 - betaU35) + 3*Math.exp((-betaU13) 
                        - betaU14 - betaU34 - betaU35) - 3*Math.exp((-du023) - betaU13 - betaU14 - betaU34 - 
                        betaU35) - Math.exp((-betaU12) - betaU13 - betaU14 - betaU34 - betaU35) + 
                        Math.exp((-du023) - betaU12 - betaU13 - betaU14 - betaU34 - betaU35) + 
                            3*Math.exp((-betaU13) - betaU15 - betaU34 - betaU35) - 3*Math.exp
                        ((-du024) - betaU13 - betaU15 - betaU34 - betaU35) - Math.exp((-betaU12) - betaU13 - 
                        betaU15 - betaU34 - betaU35) + Math.exp((-du024) - betaU12 - betaU13 - betaU15 - betaU34 - 
                        betaU35) - 2*Math.exp((-
                            betaU12) - betaU13 - betaU23 - betaU34 - betaU35) + 2*Math.exp((-du012) - betaU12 
                        - betaU13 - betaU23 - betaU34 - betaU35) - 2*Math.exp((-betaU13) - betaU14 - betaU23 - betaU34 
                        - betaU35) + 2*Math.exp((-du023) - betaU13 - betaU14 - betaU23 - 
                            betaU34 - betaU35) - 2*Math.exp((-betaU13) - betaU15 - betaU23 - betaU34 - betaU35) + 
                        2*Math.exp((-du024) - betaU13 - betaU15 -
                             betaU23 - betaU34 - betaU35) - Math.exp((-betaU12) - betaU14 - betaU24 - betaU34 - 
                        betaU35) + Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU34 - betaU35) - 
                        Math.exp((-betaU13) - betaU14 - betaU24 - betaU34 - betaU35) + 
                        Math.exp((-du023) - betaU13 - betaU14 - betaU24 - betaU34 - betaU35) - 
                        Math.exp((-betaU13) - betaU15 - betaU24 - betaU34 - betaU35) + 
                        Math.exp((-du024) - betaU13 - betaU15 - betaU24 - betaU34 - betaU35) + 3 
                        *Math.exp((-betaU23) - betaU24 - betaU34 - betaU35) - 3 
                        *Math.exp((-du123) - betaU23 - betaU24 - 
                            betaU34 - betaU35) - Math.exp((-betaU12) - betaU23 - betaU24 - betaU34 - betaU35) + 
                        Math.exp((-du123) - betaU12 - betaU23 - betaU24 - betaU34 - betaU35) - 2 
                        *Math.exp((-betaU13) - betaU23 - betaU24 - betaU34 - betaU35) + 2 
                        *Math.exp((-du123) - betaU13 - betaU23 - betaU24 - betaU34 - betaU35) - 
                        Math.exp((-betaU14) - betaU23 - betaU24 - betaU34 - betaU35) + 
                        Math.exp((-du123) - betaU14 - betaU23 - betaU24 - betaU34 - betaU35) + 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU23 - betaU24 - betaU34 -
                             betaU35) - Math.exp((-du012) - du0123 - du013 - du023 - du123 - 
                        betaU12 - betaU13 - betaU14 - betaU23 - betaU24 - betaU34 - betaU35) - Math.exp((-betaU15) - betaU23 
                        - betaU24 - betaU34 - betaU35) + Math.exp((-du123) - betaU15 - betaU23 - betaU24 - betaU34 - 
                        betaU35) + Math.exp((-betaU13) - betaU15 - betaU23 - betaU24 - betaU34 - betaU35) - 
                        Math.exp((-du024) - du123 - betaU13 - betaU15 - betaU23 - betaU24 - betaU34 - betaU35) - 
                        Math.exp((-betaU13) - betaU14 - betaU25 - betaU34 - betaU35) + 
                        Math.exp((-du023) - betaU13 - betaU14 - betaU25 - betaU34 - betaU35) - 
                        Math.exp((-betaU12) - betaU15 - betaU25 - betaU34 - betaU35) + 
                        Math.exp((-du014) - 
                            betaU12 - betaU15 - betaU25 - betaU34 - betaU35) - Math.exp((-betaU13) - betaU15 - betaU25 - 
                        betaU34 - betaU35) + Math.exp((-du024) - betaU13 - betaU15 - betaU25 - betaU34 - betaU35) + 
                        3*Math.exp((-betaU23) - betaU25 - betaU34 - betaU35) - 3 
                        *Math.exp((-du124) - 
                            betaU23 - betaU25 - betaU34 - betaU35) - Math.exp((-betaU12) - betaU23 - betaU25 - betaU34 - 
                        betaU35) + Math.exp((-du124) - betaU12 - betaU23 - betaU25 - betaU34 - betaU35) - 2 
                        *Math.exp((-betaU13) - betaU23 - betaU25 - betaU34 - betaU35) + 2 
                        *Math.exp((-du124) - 
                            betaU13 - betaU23 - betaU25 - betaU34 - betaU35) - Math.exp((-betaU14) - betaU23 - betaU25 - 
                        betaU34 - betaU35) + Math.exp((-du124) - betaU14 - betaU23 - betaU25 - betaU34 - betaU35) + 
                        Math.exp((-betaU13) - betaU14 - betaU23 - betaU25 - betaU34 -
                         betaU35) - Math.exp((-du023) - du124 - betaU13 - betaU14 - betaU23 - betaU25 - betaU34 
                        - betaU35) - Math.exp((-betaU15) - betaU23 - betaU25 - betaU34 - betaU35) + 
                        Math.exp((-du124) - betaU15 - betaU23 - betaU25 - betaU34 - betaU35) + 
                        Math.exp((-betaU12) - betaU13 - betaU15 - betaU23 - betaU25 - betaU34 - betaU35) - 
                        Math.exp((-du012) - du0124 - du014 - du024 - du124 - betaU12 - betaU13 - 
                        betaU15 - betaU23 -
                             betaU25 - betaU34 - betaU35) - 6*Math.exp((-betaU14) - betaU15 - betaU45) + 6 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU45) + 
                            3*Math.exp((-betaU12) - betaU14 - betaU15 - betaU45) - 3*Math.exp
                        ((-du034) - betaU12 - betaU14 - betaU15 - betaU45) + 3*Math.exp((-betaU13) - betaU14 
                        - betaU15 - betaU45) - 3*Math.exp((-du034) - betaU13 - betaU14 - betaU15 - betaU45) - 
                        2*Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU45) + 2*Math.exp
                        ((-du034) - betaU12 - betaU13 - betaU14 - betaU15 - betaU45) + 
                            2*Math.exp((-betaU12) - betaU13 - betaU23 - betaU45) - 2*Math.exp
                        ((-du012) - betaU12 - betaU13 - betaU23 - betaU45) - Math.exp((-betaU12) - betaU13 - 
                        betaU14 - betaU23 - betaU45) + Math.exp((-du012) - betaU12 - betaU13 - betaU14 - betaU23 - 
                        betaU45) - Math.exp((-betaU12) - betaU13 - betaU15 - betaU23 - betaU45) + 
                        Math.exp((-du012) - betaU12 - betaU13 - betaU15 - betaU23 - betaU45) + 2 
                        *Math.exp((-betaU14) - betaU15 - betaU23 - betaU45) - 2 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU23 - betaU45) - Math.exp(
                        (-betaU12) - betaU14 - betaU15 - betaU23 - betaU45) + Math.exp((-du034) - betaU12 - 
                        betaU14 - betaU15 - betaU23 - betaU45) - Math.exp((-betaU13) - betaU14 - betaU15 - betaU23 - 
                        betaU45) + Math.exp((-du034) - betaU13 - betaU14 - betaU15 - betaU23 - betaU45) + 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU23 - betaU45) - Math.exp
                        ((-du012) - du034 - betaU12 - betaU13 - betaU14 - 
                            betaU15 - betaU23 - betaU45) + 3*Math.exp((-betaU12) - betaU14 - betaU24 - betaU45) - 
                        3*Math.exp((-du013) - betaU12 -
                             betaU14 - betaU24 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU14 - betaU24 - 
                        betaU45) + Math.exp((-du013) - betaU12 - betaU13 - betaU14 - 
                            betaU24 - betaU45) + 3*Math.exp((-betaU14) - betaU15 - betaU24 - betaU45) - 
                            3*Math.exp((-du034) - betaU14 - betaU15 - betaU24 - betaU45) - 
                        Math.exp((-betaU13) - betaU14 - betaU15 - betaU24 - betaU45) + 
                        Math.exp((-du034) - betaU13 - 
                        betaU14 - betaU15 - betaU24 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU23 - betaU24 - betaU45
                        ) + Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU24 - betaU45) - 
                        Math.exp((-betaU12) - betaU14 - betaU23 - betaU24 - betaU45) + 
                        Math.exp((-du013) - betaU12 - 
                        betaU14 - betaU23 - betaU24 - betaU45) - Math.exp((-betaU14) - betaU15 - betaU23 - betaU24 - betaU45
                        ) + Math.exp((-du034) - betaU14 - betaU15 - betaU23 - betaU24 - betaU45) + 3 
                        *Math.exp((-betaU12) - betaU15 - betaU25 - betaU45) - 3 
                        *Math.exp((-du014) - betaU12 - 
                        betaU15 - betaU25 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU15 - betaU25 - betaU45) + 
                        Math.exp((-du014) - betaU12 - betaU13 - betaU15 - betaU25 - betaU45) + 3 
                        *Math.exp((-betaU14) - betaU15 - betaU25 - betaU45) - 3 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU25 - betaU45) - Math.exp(
                        (-betaU13) - betaU14 - betaU15 - betaU25 - betaU45) + Math.exp((-du034) - betaU13 - 
                        betaU14 - betaU15 - betaU25 -
                             betaU45) - Math.exp((-betaU12) - betaU13 - betaU23 - betaU25 - betaU45) + 
                        Math.exp((-du012) - betaU12 - betaU13 - 
                            betaU23 - betaU25 - betaU45) - Math.exp((-betaU12) - betaU15 - betaU23 - betaU25 - betaU45) 
                        + Math.exp((-du014) - betaU12 - betaU15 - betaU23 -
                             betaU25 - betaU45) - Math.exp((-betaU14) - betaU15 - betaU23 - betaU25 - betaU45) + 
                        Math.exp((-du034) - betaU14 - betaU15 - 
                            betaU23 - betaU25 - betaU45) - 6*Math.exp((-betaU24) - betaU25 - betaU45) + 6 
                        *Math.exp((-du134) - betaU24 - betaU25 - betaU45) + 3 
                        *Math.exp((-betaU12) - betaU24 - betaU25 - betaU45) - 3 
                        *Math.exp((-du134) - betaU12 - betaU24 - betaU25 - betaU45) + 
                            2*Math.exp((-betaU13) - betaU24 - betaU25 - betaU45) - 2*Math.exp
                        ((-du134) - betaU13 - betaU24 - betaU25 - betaU45) - Math.exp((-betaU12) - betaU13 - 
                        betaU24 - betaU25 - betaU45) + Math.exp((-du134) - betaU12 - betaU13 - betaU24 - betaU25 - 
                        betaU45) + 3*Math.exp((-
                            betaU14) - betaU24 - betaU25 - betaU45) - 3*Math.exp((-du134) - betaU14 - betaU24 
                        - betaU25 - betaU45) - Math.exp((-betaU13) - betaU14 - betaU24 - betaU25 -
                             betaU45) + Math.exp((-du134) - betaU13 - betaU14 - betaU24 - betaU25 - betaU45) + 
                        3*Math.exp((-betaU15) - betaU24 - betaU25 -
                         betaU45) - 3*Math.exp((-du134) - betaU15 - betaU24 - betaU25 - betaU45) - 
                        Math.exp((-betaU13) - betaU15 - 
                        betaU24 - betaU25 - betaU45) + Math.exp((-du134) - betaU13 - betaU15 - betaU24 - betaU25 - 
                        betaU45) - 3*Math.exp((-betaU12) -
                         betaU14 - betaU15 - betaU24 - betaU25 - betaU45) + 3*Math.exp((-du013) - du0134 - 
                        du014 - du034 - du134 - betaU12 - betaU14 -
                         betaU15 - betaU24 - betaU25 - betaU45) + Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - 
                        betaU24 - betaU25 - betaU45) - Math.exp((-du013) - du0134 - du014 - du034 - 
                        du134 - betaU12 - betaU13 - betaU14 - betaU15 - betaU24 - betaU25 - betaU45) + 3*Math.exp((-
                            betaU23) - betaU24 - betaU25 - betaU45) - 3*Math.exp((-du134) - betaU23 - betaU24 
                        - betaU25 - betaU45) - 2*Math.exp((-betaU12) - betaU23 - betaU24 - 
                            betaU25 - betaU45) + 2*Math.exp((-du134) - betaU12 - betaU23 - betaU24 - betaU25 - 
                        betaU45) - Math.exp((-betaU13) - betaU23 - betaU24 - betaU25 - betaU45) + 
                        Math.exp((-du134) - betaU13 - betaU23 - betaU24 - betaU25 - betaU45) + 
                        Math.exp((-betaU12) - betaU13 - betaU23 - betaU24 - betaU25 - betaU45) - Math.exp
                        ((-du012) - du134 - betaU12 - betaU13 - betaU23 - betaU24 - betaU25 - betaU45) - 
                        Math.exp((-betaU14) - betaU23 - betaU24 - betaU25 - betaU45) + 
                        Math.exp((-du134) - betaU14 - betaU23 - betaU24 - betaU25 - betaU45) - 
                        Math.exp((-betaU15) - betaU23 - betaU24 - betaU25 - betaU45) + 
                        Math.exp((-du134) - 
                        betaU15 - betaU23 - betaU24 - betaU25 - betaU45) + Math.exp((-betaU12) - betaU14 - betaU15 - betaU23 
                        - betaU24 - betaU25 - betaU45) - Math.exp((-du013) - du0134 - du014 - du034 
                        - du134 - betaU12 - betaU14 - betaU15 - betaU23 - betaU24 - betaU25 - betaU45) + 3 
                        *Math.exp((-betaU13) - betaU14 - betaU34 - betaU45) - 3 
                        *Math.exp((-du023) - 
                            betaU13 - betaU14 - betaU34 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU14 - betaU34 - 
                        betaU45) + Math.exp((-du023) - betaU12 - betaU13 - betaU14 - betaU34 - betaU45) + 3 
                        *Math.exp((-betaU14) - betaU15 - betaU34 - betaU45) - 3 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU34 - 
                            betaU45) - Math.exp((-betaU12) - betaU14 - betaU15 - betaU34 - betaU45) + 
                        Math.exp((-du034) - betaU12 - betaU14 - betaU15 - betaU34 - betaU45) - 
                        Math.exp((-betaU12) - betaU13 - betaU23 - betaU34 - betaU45) + 
                        Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU34 - betaU45) - 
                        Math.exp((-betaU13) - betaU14 - betaU23 - betaU34 - 
                            betaU45) + Math.exp((-du023) - betaU13 - betaU14 - betaU23 - betaU34 - betaU45) - 
                        Math.exp((-betaU14) - betaU15 - 
                            betaU23 - betaU34 - betaU45) + Math.exp((-du034) - betaU14 - betaU15 - betaU23 - betaU34 
                        - betaU45) - 2*Math.exp((-betaU12) - betaU14 - betaU24 -
                             betaU34 - betaU45) + 2*Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU34 - 
                        betaU45) - 2*Math.exp((-betaU13) - betaU14 - betaU24 - betaU34 - betaU45) + 2 
                        *Math.exp((-du023) - betaU13 - betaU14 - betaU24 - betaU34 - betaU45) - 2 
                        *Math.exp((-betaU14) - betaU15 - betaU24 - betaU34 - betaU45) + 2 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU24 - betaU34 - betaU45) + 3 
                        *Math.exp((-betaU23) - betaU24 - betaU34 - betaU45) - 3 
                        *Math.exp((-du123) - betaU23 - betaU24 - betaU34 - betaU45) - Math.exp(
                        (-betaU12) - betaU23 - betaU24 - betaU34 - 
                            betaU45) + Math.exp((-du123) - betaU12 - betaU23 - betaU24 - betaU34 - betaU45) - 
                        Math.exp((-betaU13) - betaU23 - betaU24 - 
                            betaU34 - betaU45) + Math.exp((-du123) - betaU13 - betaU23 - betaU24 - betaU34 - betaU45
                        ) - 2*Math.exp((-betaU14) - betaU23 - betaU24 - betaU34 - betaU45) + 2 
                        *Math.exp((-du123) - betaU14 - betaU23 - betaU24 - betaU34 - betaU45) + 
                        Math.exp((-betaU12) - betaU13 - betaU14 - betaU23 - betaU24 - betaU34 - betaU45) - 
                        Math.exp((-du012) - du0123 - du013 - 
                            du023 - du123 - betaU12 - betaU13 - betaU14 - 
                        betaU23 - betaU24 - betaU34 - betaU45) - Math.exp((-betaU15) - betaU23 - betaU24 - betaU34 - betaU45
                        ) + Math.exp((-du123) - betaU15 - betaU23 - betaU24 - betaU34 - betaU45) + 
                        Math.exp((-betaU14) - betaU15 - betaU23 - betaU24 - betaU34 - betaU45) - Math.exp
                        ((-du034) - du123 - betaU14 - betaU15 - betaU23 - betaU24 - betaU34 - betaU45) - 
                        Math.exp((-betaU13) - betaU14 -
                             betaU25 - betaU34 - betaU45) + Math.exp((-du023) - betaU13 - betaU14 - betaU25 - 
                        betaU34 - betaU45) - Math.exp((-betaU12) - 
                            betaU15 - betaU25 - betaU34 - betaU45) + Math.exp((-du014) - betaU12 - betaU15 - betaU25 
                        - betaU34 - betaU45) - Math.exp((-betaU14) - betaU15 -
                             betaU25 - betaU34 - betaU45) + Math.exp((-du034) - betaU14 - betaU15 - betaU25 - 
                        betaU34 - betaU45) + 3*Math.exp((-betaU24) - betaU25 - 
                            betaU34 - betaU45) - 3*Math.exp((-du134) - betaU24 - betaU25 - betaU34 - betaU45) 
                        - Math.exp((-betaU12) - betaU24 - betaU25 - betaU34 - betaU45) + Math.exp(
                        (-du134) - betaU12 - betaU24 - betaU25 - betaU34 - betaU45) - Math.exp((-betaU13) - 
                        betaU24 - betaU25 - betaU34 - betaU45) + Math.exp((-du134) - betaU13 - betaU24 - betaU25 - 
                        betaU34 - betaU45) - 2*Math.exp((-betaU14) - betaU24 -
                         betaU25 - betaU34 - betaU45) + 2*Math.exp((-du134) - betaU14 - betaU24 - betaU25 - betaU34 
                        - betaU45) + Math.exp((-betaU13) - betaU14 - betaU24 - betaU25 - betaU34 - betaU45) - 
                        Math.exp((-du023) - du134 - betaU13 - betaU14 - betaU24 -
                             betaU25 - betaU34 - betaU45) - Math.exp((-betaU15) - betaU24 - betaU25 - betaU34 - 
                        betaU45) + Math.exp((-du134) - betaU15 - betaU24 - 
                            betaU25 - betaU34 - betaU45) + Math.exp((-betaU12) - betaU14 - betaU15 - betaU24 - betaU25 - 
                        betaU34 - betaU45) - Math.exp((-du013) - du0134 - du014 - du034 - du134 
                        - betaU12 - betaU14 - betaU15 - betaU24 - betaU25 - betaU34 - betaU45) + 3*Math.exp((-betaU13) 
                        - betaU15 - betaU35 - betaU45) - 3*Math.exp((-du024) - betaU13 - betaU15 - betaU35 - 
                        betaU45) - Math.exp((-betaU12) - betaU13 - betaU15 - betaU35 - betaU45) + 
                        Math.exp((-du024) - betaU12 - betaU13 - betaU15 - betaU35 - betaU45) + 3 
                        *Math.exp((-betaU14) - betaU15 - betaU35 - betaU45) - 3 
                        *Math.exp((-du034) - betaU14 - betaU15 - betaU35 - betaU45) - Math.exp(
                        (-betaU12) - betaU14 - betaU15 - betaU35 - betaU45) + Math.exp((-du034) - betaU12 - 
                        betaU14 - betaU15 - betaU35 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU23 - betaU35 - betaU45
                        ) + Math.exp((-du012) - betaU12 - betaU13 - betaU23 - betaU35 - betaU45) - 
                        Math.exp((-betaU13) - betaU15 - betaU23 - betaU35 - betaU45) + 
                        Math.exp((-du024) - betaU13 - betaU15 - betaU23 - betaU35 - betaU45) - 
                        Math.exp((-betaU14) - betaU15 - betaU23 - betaU35 - betaU45) + 
                        Math.exp((-du034) - betaU14 - betaU15 - betaU23 - betaU35 - betaU45) - 
                        Math.exp((-betaU12) - betaU14 - betaU24 - betaU35 - betaU45) + 
                        Math.exp((-du013) - betaU12 - betaU14 - betaU24 - betaU35 - betaU45) - 
                        Math.exp((-betaU13) - betaU15 - betaU24 - betaU35 - betaU45) + 
                        Math.exp((-du024) - betaU13 - betaU15 - betaU24 - betaU35 - betaU45) - 
                        Math.exp((-betaU14) - betaU15 - betaU24 - betaU35 - betaU45) + 
                        Math.exp((-du034) - betaU14 - betaU15 - betaU24 - betaU35 - betaU45) - 2 
                        *Math.exp((-betaU12) - betaU15 - betaU25 - betaU35 - betaU45) + 2 
                        *Math.exp((-du014) - betaU12 - betaU15 -
                         betaU25 - betaU35 - betaU45) - 2*Math.exp((-betaU13) - betaU15 - betaU25 - betaU35 - betaU45) 
                        + 2*Math.exp((-du024) - betaU13 - betaU15 - betaU25 - betaU35 - betaU45) - 2 
                        *Math.exp((-betaU14) - betaU15 - betaU25 - betaU35 - betaU45) + 2 
                        *Math.exp((-du034) - betaU14 - betaU15 -
                         betaU25 - betaU35 - betaU45) + 3*Math.exp((-betaU23) - betaU25 - betaU35 - betaU45) - 3 
                        *Math.exp((-
                            du124) - betaU23 - betaU25 - betaU35 - betaU45) - Math.exp((-betaU12) - betaU23 - 
                        betaU25 - betaU35 - betaU45) + Math.exp((-
                            du124) - betaU12 - betaU23 - betaU25 - betaU35 - betaU45) - Math.exp((-betaU13) - 
                        betaU23 - betaU25 - betaU35 - betaU45) + Math.exp((-du124) - betaU13 - betaU23 - betaU25 - 
                        betaU35 - betaU45) - Math.exp((-betaU14) - betaU23 - betaU25 - betaU35 - betaU45) + 
                        Math.exp((-
                            du124) - betaU14 - betaU23 - betaU25 - betaU35 - betaU45) - 2*Math.exp((-betaU15) 
                        - betaU23 - betaU25 - betaU35 - betaU45) + 2*Math.exp((-
                            du124) - betaU15 - betaU23 - 
                            betaU25 - betaU35 - betaU45) + Math.exp((-betaU12) - betaU13 - betaU15 - betaU23 - betaU25 - 
                        betaU35 - betaU45) - Math.exp((-du012) - du0124 - du014 - du024 - du124 
                        - betaU12 - betaU13 - betaU15 - betaU23 - betaU25 - betaU35 - betaU45) + Math.exp((-betaU14) - 
                        betaU15 - betaU23 - betaU25 - betaU35 - betaU45) - Math.exp((-du034) - du124 - betaU14 - 
                        betaU15 - betaU23 - betaU25 - betaU35 - betaU45) + 3*Math.exp((-betaU24) - betaU25 - betaU35 - 
                        betaU45) - 3*Math.exp((-du134) - betaU24 - betaU25 - betaU35 - betaU45) - 
                        Math.exp((-betaU12) - betaU24 - betaU25 - betaU35 - betaU45) + 
                        Math.exp((-du134) - betaU12 - betaU24 - betaU25 -
                             betaU35 - betaU45) - Math.exp((-betaU13) - betaU24 - betaU25 - betaU35 - betaU45) + 
                        Math.exp((-du134) - betaU13 - 
                            betaU24 - betaU25 - betaU35 - betaU45) - Math.exp((-betaU14) - betaU24 - betaU25 - betaU35 - 
                        betaU45) + Math.exp((-du134) - betaU14 -
                             betaU24 - betaU25 - betaU35 - betaU45) - 2*Math.exp((-betaU15) - betaU24 - betaU25 - 
                        betaU35 - betaU45) + 2*Math.exp((-du134) - betaU15 - betaU24 - betaU25 - betaU35 - 
                        betaU45) + Math.exp((-betaU13) - betaU15 - betaU24 - betaU25 - betaU35 - betaU45) - 
                        Math.exp((-du024) - du134 - betaU13 - betaU15 - betaU24 - betaU25 - betaU35 - betaU45) + 
                        Math.exp((-betaU12) - betaU14 - betaU15 -
                             betaU24 - betaU25 - betaU35 - betaU45) - Math.exp((-du013) - du0134 - du014 
                        - du034 - du134 - betaU12 - betaU14 - betaU15 - betaU24 - betaU25 - betaU35 - betaU45) - 6 
                        *Math.exp((-betaU34) - betaU35 - betaU45) + 6*Math.exp((-du234) 
                        - betaU34 - betaU35 - betaU45) + 2*Math.exp((-betaU12) - betaU34 - betaU35 - betaU45) - 2 
                        *Math.exp((-du234) - betaU12 - betaU34 - betaU35 - betaU45) + 3 
                        *Math.exp((-betaU13) - betaU34 - betaU35 - betaU45) - 3 
                        *Math.exp((-du234) - betaU13 - betaU34 - 
                            betaU35 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU12 - betaU13 - betaU34 - betaU35 - betaU45) + 3 
                        *Math.exp((-betaU14) - betaU34 - betaU35 - betaU45) - 3 
                        *Math.exp((-du234) - betaU14 - 
                            betaU34 - betaU35 - betaU45) - Math.exp((-betaU12) - betaU14 - betaU34 - betaU35 - betaU45) 
                        + Math.exp((-du234) - betaU12 - betaU14 - betaU34 - betaU35 - betaU45) + 3 
                        *Math.exp((-betaU15) - betaU34 - betaU35 - betaU45) - 3 
                        *Math.exp((-du234) - betaU15 - betaU34 - betaU35 - betaU45) - Math.exp(
                        (-betaU12) - 
                            betaU15 - betaU34 - betaU35 - betaU45) + Math.exp((-du234) - betaU12 - betaU15 - betaU34 
                        - betaU35 - betaU45) - 3*Math.exp((-betaU13) - betaU14 - betaU15 - 
                        betaU34 - betaU35 - betaU45) + 3*Math.exp((-du023) - du0234 - du024 - du034 
                        - du234 - betaU13 - 
                            betaU14 - betaU15 - betaU34 - betaU35 - betaU45) + Math.exp((-betaU12) - betaU13 - betaU14 - 
                        betaU15 - betaU34 - betaU35 - betaU45) - Math.exp((-du023) - du0234 - du024 - 
                        du034 - du234 - betaU12 - betaU13 - betaU14 - betaU15 - betaU34 - betaU35 - betaU45) + 
                            3*Math.exp((-betaU23) - 
                        betaU34 - betaU35 - betaU45) - 3*Math.exp((-du234) - betaU23 - betaU34 - betaU35 - 
                        betaU45) - Math.exp((-
                            betaU12) - betaU23 - betaU34 - betaU35 - betaU45) + Math.exp((-du234) - betaU12 - 
                        betaU23 - betaU34 - betaU35 - betaU45) - 2*Math.exp((-betaU13) - betaU23 - betaU34 - betaU35 - 
                        betaU45) + 2*Math.exp((-du234) - betaU13 - betaU23 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-betaU12) - betaU13 - betaU23 - betaU34 - betaU35 - betaU45) - Math.exp
                        ((-du012) - du234 - betaU12 - betaU13 - betaU23 - 
                            betaU34 - betaU35 - betaU45) - Math.exp((-betaU14) - betaU23 - betaU34 - betaU35 - betaU45) 
                        + Math.exp((-du234) - betaU14 - betaU23 - betaU34 - betaU35 - 
                            betaU45) - Math.exp((-betaU15) - betaU23 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU15 - betaU23 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-betaU13) - betaU14 - betaU15 - betaU23 - betaU34 - betaU35 - betaU45) - 
                        Math.exp((-du023) - du0234 - du024 - du034 - du234 - betaU13 - betaU14 - 
                        betaU15 - betaU23 - betaU34 - betaU35 - betaU45) + 3*Math.exp((-betaU24) - 
                            betaU34 - betaU35 - betaU45) - 3*Math.exp((-du234) - betaU24 - betaU34 - betaU35 - 
                        betaU45) - Math.exp((-betaU12) - betaU24 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU12 - betaU24 - betaU34 - betaU35 - betaU45) - 
                        Math.exp((-betaU13) - betaU24 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU13 - betaU24 - betaU34 - betaU35 - betaU45) - 2 
                        *Math.exp((-betaU14) - betaU24 - betaU34 - betaU35 - betaU45) + 2 
                        *Math.exp((-du234) - betaU14 - betaU24 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-betaU12) - betaU14 - betaU24 - betaU34 - betaU35 - betaU45) - Math.exp
                        ((-du013) - du234 - betaU12 - betaU14 - betaU24 - betaU34 - betaU35 - 
                            betaU45) - Math.exp((-betaU15) - betaU24 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU15 - betaU24 -
                             betaU34 - betaU35 - betaU45) + Math.exp((-betaU13) - betaU14 - betaU15 - betaU24 - betaU34 
                        - betaU35 - betaU45) - Math.exp((-du023) -
                             du0234 - du024 - du034 - du234 -
                         betaU13 - betaU14 - betaU15 - betaU24 - betaU34 - betaU35 - betaU45) + 3*Math.exp((-betaU25) - 
                        betaU34 - betaU35 - betaU45) - 3*Math.exp((-du234) - betaU25 - betaU34 - betaU35 - 
                        betaU45) - Math.exp((-betaU12) - 
                            betaU25 - betaU34 - betaU35 - betaU45) + Math.exp((-du234) - betaU12 - betaU25 - betaU34 
                        - betaU35 - betaU45) - Math.exp((-betaU13) - betaU25 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU13 - betaU25 - betaU34 - betaU35 - betaU45) - 
                        Math.exp((-betaU14) - betaU25 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-du234) - betaU14 - betaU25 - betaU34 - betaU35 - betaU45) - 2 
                        *Math.exp((-betaU15) - betaU25 - betaU34 - betaU35 - betaU45) + 2 
                        *Math.exp((-du234) - betaU15 - betaU25 - betaU34 - betaU35 - betaU45) + 
                        Math.exp((-betaU12) - betaU15 - betaU25 -
                         betaU34 - betaU35 - betaU45) - Math.exp((-du014) - du234 - betaU12 - betaU15 - betaU25 -
                             betaU34 - betaU35 - betaU45) + Math.exp((-betaU13) - betaU14 - betaU15 - betaU25 - betaU34 
                        - betaU35 - betaU45) - Math.exp((-du023) - du0234 -
                             du024 - du034 - du234 - betaU13 - betaU14 - betaU15 - betaU25 - betaU34 - betaU35 - betaU45) - 3 
                        *Math.exp((-betaU23) - 
                            betaU24 - betaU25 - betaU34 - betaU35 - betaU45) + 3*Math.exp((-du123) - du1234 
                        - du124 - du134 - du234 - betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - 
                            betaU45) + Math.exp((-betaU12) - betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - betaU45) 
                        - Math.exp((-du123) - du1234 - du124 - du134 - du234 - betaU12 - betaU23 
                        - betaU24 - betaU25 - betaU34 - betaU35 - betaU45) + Math.exp((-betaU13) - 
                            betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - 
                        betaU45) - Math.exp((-du123) - du1234 - du124 - du134 - du234 - betaU13 
                        - betaU23 - betaU24 - betaU25 - betaU34 - 
                            betaU35 - betaU45) + Math.exp((-betaU14) - betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - 
                        betaU45) - Math.exp((-du123) - du1234 - du124 - du134 - du234 - betaU14 
                        - betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - betaU45) + Math.exp((-betaU15) - 
                            betaU23 - betaU24 - betaU25 - betaU34 - betaU35 - betaU45) - Math.exp((-du123) - 
                        du1234 - du124 - du134 - du234 - betaU15 - betaU23 - betaU24 - betaU25 - betaU34 - 
                            betaU35 - betaU45) - Math.exp((-betaU12) - betaU13 - betaU14 - betaU15 - betaU23 - betaU24 - 
                        betaU25 - betaU34 - betaU35 - betaU45) + Math.exp((-du012) - du0123 - du01234 - 
                        du0124 - du013 - du0134 - du014 - du023 - du0234 - du024 - du034 - du123 - 
                        du1234 - du124 - du134 - du234 - betaU12 - betaU13 - betaU14 - betaU15 - betaU23 - betaU24 - betaU25 - 
                        betaU34 - betaU35 - betaU45));

                    //System.out.println("deltaE = " + deltaE);
                    // coefficient is -(1-5)/5! = -1/30
                        value += -deltaE/30.0;
                }
        }

        return value;
    }
    
    public ClusterBonds[] getClusters() {return clusters;}
    /**
     * @return Returns the temperature.
     */
    public double getTemperature() {
        return 1/beta;
    }
    /**
     * @param temperature The temperature to set.
     */
    public void setTemperature(double temperature) {
        beta = 1/temperature;
    }
    
    public void setDeltaDCut(double newDeltaDCut) {
        deltaDCut2 = newDeltaDCut*newDeltaDCut;
    }
    
    public double getDeltaDCut() {
        return Math.sqrt(deltaDCut2);
    }
    
    private static final long serialVersionUID = 1L;
    private final ClusterBonds[] clusters;
    private final double[] clusterWeights;
    private final MayerFunction[] f;
    private final double[][][] fValues;
    private int cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
    private double beta;
    protected final AtomArrayList scfAtoms;
    protected double deltaDCut2 = Double.POSITIVE_INFINITY;
    protected final double[][] uijPol;
}
