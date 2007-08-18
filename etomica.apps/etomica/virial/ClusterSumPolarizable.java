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
        double deltaC = 0.0;
        double deltaD = 0.0;
        double u12 = 0.0;
        double u13 = 0.0;
        double u14 = 0.0;
        double u23 = 0.0;
        double u24 = 0.0;
        double u34 = 0.0;
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

        if (nPoints > 2) {
	        u12 = (-1)/beta*Math.log(fValues[1][0][0]+1);
	        u13 = (-1)/beta*Math.log(fValues[2][0][0]+1);
	        u23 = (-1)/beta*Math.log(fValues[2][1][0]+1);
        }
        
        if (nPoints > 3) {
	        u14 = (-1)/beta*Math.log(fValues[3][0][0]+1);
	        u24 = (-1)/beta*Math.log(fValues[3][1][0]+1);
	        u34 = (-1)/beta*Math.log(fValues[3][2][0]+1);
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
    		if (u12+u13+u23 == Double.POSITIVE_INFINITY) {
    			//System.out.println("Sum of pair energies is infinity: u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23);
    			deltaC = 0.0;
    		}
    		else {
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
    			deltaC = expBetaU123*Math.exp(-beta*(u12 + u13 + u23));       
    			
    			// deltaC has to be multiplied by clusterWeights, just like v was multiplied by
    			// clusterWeights above to get value
    			deltaC = deltaC*clusterWeights[0];
    			
    			if (Double.isInfinite(deltaC)) {
    				System.out.println("deltaC = " + deltaC);
    			}
    		}
	        		
    		//System.out.println("u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23 + ", u123 = " + u123 + ", deltaC = " + deltaC + ", value(before) = " + value);
    		value += deltaC;
		}

		if (nPoints == 4) {
			
			if (u12+u13+u14+u23+u24+u34 == Double.POSITIVE_INFINITY) {
    			//System.out.println("Sum of pair energies is infinity: u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23);
				deltaD = 0.0;
    		}
			else {
			    scfAtoms.clear();
		        scfAtoms.add(box.molecule(0));
                scfAtoms.add(box.molecule(1));
                scfAtoms.add(box.molecule(2));
                
				double u123Pol = scfPotential.getPolarizationEnergy(scfAtoms);
				scfAtoms.remove(2);
				scfAtoms.add(box.molecule(3));
				double u124Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                scfAtoms.remove(1);
                scfAtoms.add(box.molecule(2));
                double u134Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                scfAtoms.remove(0);
                scfAtoms.add(box.molecule(1));
                double u234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                scfAtoms.add(box.molecule(0));
                double u1234Pol = scfPotential.getPolarizationEnergy(scfAtoms);

				double deltaU123 = u123Pol-u12Pol-u13Pol-u23Pol;
				double deltaU124 = u124Pol-u12Pol-u14Pol-u24Pol;
				double deltaU134 = u134Pol-u13Pol-u14Pol-u34Pol;
				double deltaU234 = u234Pol-u23Pol-u24Pol-u34Pol;

				double deltaU1234 = u1234Pol-u12Pol-u13Pol-u14Pol-u23Pol-u24Pol-u34Pol-deltaU123-deltaU124-deltaU134-deltaU234;
				
				double g12 = Math.exp(-beta*u12);
				double g13 = Math.exp(-beta*u13);
				double g14 = Math.exp(-beta*u14);
				double g23 = Math.exp(-beta*u23);
				double g24 = Math.exp(-beta*u24);
				double g34 = Math.exp(-beta*u34);

				deltaD = (Math.exp(-beta*(deltaU123+deltaU124+deltaU134+deltaU234+deltaU1234))-1)*g12*g13*g14*g23*g24*g34
				        -(Math.exp(-beta*deltaU123)-1)*g12*g13*g23*(g14+g24+g34-2)
				        -(Math.exp(-beta*deltaU124)-1)*g12*g14*g24*(g13+g23+g34-2)
				        -(Math.exp(-beta*deltaU134)-1)*g13*g14*g34*(g12+g23+g24-2)
				        -(Math.exp(-beta*deltaU234)-1)*g23*g24*g34*(g12+g13+g14-2);

        				// Mason and Spurling book deltaD; 5/11/07

//					deltaU1234 = u1234-u12-u13-u14-u23-u24-u34;

/*					deltaD = 2*Math.exp(-beta*(u12+u13+u23))*(Math.exp(-beta*deltaU123)-1) + 2*Math.exp(-beta*(u12+u14+u24))*(Math.exp(-beta*deltaU124)-1)
						+ 2*Math.exp(-beta*(u13+u14+u34))*(Math.exp(-beta*deltaU134)-1) + 2*Math.exp(-beta*(u23+u24+u34))*(Math.exp(-beta*deltaU234)-1)
						+ Math.exp(-beta*(u12+u13+u23-u14))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u13+u23-u24))*(1-Math.exp(-beta*deltaU123))
						+ Math.exp(-beta*(u12+u13+u23-u34))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u14+u24-u13))*(1-Math.exp(-beta*deltaU124))
						+ Math.exp(-beta*(u12+u14+u24-u23))*(1-Math.exp(-beta*deltaU124)) + Math.exp(-beta*(u12+u14+u24-u34))*(1-Math.exp(-beta*deltaU124))
						+ Math.exp(-beta*(u13+u14+u34-u12))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u13+u14+u34-u23))*(1-Math.exp(-beta*deltaU134))
						+ Math.exp(-beta*(u13+u14+u34-u24))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u23+u24+u34-u12))*(1-Math.exp(-beta*deltaU234))
						+ Math.exp(-beta*(u23+u24+u34-u13))*(1-Math.exp(-beta*deltaU234)) + Math.exp(-beta*(u23+u24+u34-u14))*(1-Math.exp(-beta*deltaU234))
						+ Math.exp(-beta*(u12+u13+u14+u23+u24+u34))*(Math.exp(-beta*deltaU1234)-1);
*/

//        				 kmb added this line 8/16/06
	        			// deltaD has to be multiplied by weightPrefactor from Standard class, just like deltaC was multiplied by
	        			// clusterWeights above to get value; note, for B3 clusterWeights = weightPrefactor
        				
				deltaD = -0.125*deltaD;  //XXX - IS THIS DIFFERENT NOW?; KMB, 7/27/07  -32 is value of weightPrefactor in Standard for a B4 calc
			}
        			
			value += deltaD;

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
    
    private static final long serialVersionUID = 1L;
    private final ClusterBonds[] clusters;
    private final double[] clusterWeights;
    private final MayerFunction[] f;
    private final double[][][] fValues;
    private int cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
    private double beta;
    protected final AtomArrayList scfAtoms;
}
