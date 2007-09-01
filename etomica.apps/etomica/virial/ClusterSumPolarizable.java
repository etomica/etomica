package etomica.virial;

import etomica.atom.AtomArrayList;
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
        
        // kmb added 7/10/06
        double u12Pol = Double.NaN;
        double u13Pol = Double.NaN;
        double u23Pol = Double.NaN;
        double u14Pol = Double.NaN;
        double u24Pol = Double.NaN;
        double u34Pol = Double.NaN;
        // recalculate all f values for all pairs
        PotentialPolarizable scfPotential = null;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    fValues[i][j][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    scfPotential = (PotentialPolarizable) f[k].getPotential();
                    // NEED k==0 IN BOOLEAN FOR B4 CALC, TO PREVENT OVERWRITING THE ATOM LISTS FOR THE
                    // TRIMER AND TETRAMER SCF ENERGY CALCS; KMB, 8/16/06
                    if (nPoints > 2 && k == 0) {
                        if (i==0) {
                            if (j==1) {
                                u12Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else if (j==2) {
                                u13Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else if (j==3) {
                                u14Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else {
                                throw new RuntimeException("Need to add B5 handling");
                            }
                        }
                        else if (i==1) {
                            if (j==2) {
                                u23Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else if (j==3) {
                                u24Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else {
                                throw new RuntimeException("Need to add B5 handling");
                            }
                        }
                        else if (i==2) {
                            if (j==3) {
                                u34Pol = scfPotential.getLastPolarizationEnergy();
                            }
                            else {
                                throw new RuntimeException("Need to add B5 handling");
                            }
                        }
                        else {
                            throw new RuntimeException("Need to add B5 handling");
                        }
                    }

                    fValues[j][i][k] = fValues[i][j][k];
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
    			double deltau123 = u123Pol-(u12Pol+u13Pol+u23Pol);
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

        if (nPoints == 4) {
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
                    double deltaU123 = u123Pol-u12Pol-u13Pol-u23Pol;
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
                    double deltaU124 = u124Pol-u12Pol-u14Pol-u24Pol;
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
                    double deltaU134 = u134Pol-u13Pol-u14Pol-u34Pol;
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
                    double deltaU234 = u234Pol-u23Pol-u24Pol-u34Pol;
                    double beta234 = beta*deltaU234;
                    double exp234 = -beta234;
                    if (Math.abs(beta234) > 1E-8) {
                        exp234 = Math.exp(-beta234) - 1;
                    }
                    deltaD += -exp234*g23*g24*g34*((f12+f13+f14)+1);
                }

                if (g12*g13*g14*g23*g24*g34 != 0) {
                    scfAtoms.add(box.molecule(0));
                    double u1234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                    // deltaU1234 would have deltaUabc subtracted off, but we'd also add it back
                    // in for expU1234, so just don't subtract in the first place 
    				double deltaU1234 = u1234Pol-u12Pol-u13Pol-u14Pol-u23Pol-u24Pol-u34Pol; //-deltaU123-deltaU124-deltaU134-deltaU234;
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
}
