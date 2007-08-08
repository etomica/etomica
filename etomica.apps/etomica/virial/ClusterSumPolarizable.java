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
        System.out.println("Hey, you initialized maxOOPiValues array!  Only once, I hope!");
        maxOOPiValues = new double[9];
        maxOOdistances = new double[9];
        binCounts = new double[9];
        for(int z=0; z<9; z++) {
        		maxOOPiValues[z] = Double.NEGATIVE_INFINITY;//0;
        		maxOOdistances[z] = 0;
        		binCounts[z] = 0;
        }

        scfAtoms = new AtomArrayList(); // USE THIS LIST FOR ALL ATOMS, WHETHER 3 OR 4; KMB, 8/16/06
        scfAtoms123 = new AtomArrayList();
        scfAtoms124 = new AtomArrayList();
        scfAtoms134 = new AtomArrayList();
        scfAtoms234 = new AtomArrayList();
    }

    // equal point count enforced in constructor 
    public int pointCount() {
        return clusters[0].pointCount();
    }
    
    public ClusterAbstract makeCopy() {
        return new ClusterSumPolarizable(clusters,clusterWeights,f);
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
        double u123 = 0.0;
        double u124 = 0.0;
        double u134 = 0.0;
        double u234 = 0.0;
        double u1234 = 0.0;
        double u12Pol = Double.NaN;
        double u13Pol = Double.NaN;
        double u23Pol = Double.NaN;
        // recalculate all f values for all pairs
        PotentialPolarizable scfPotential = null;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                for(int k=0; k<f.length; k++) {
                    fValues[i][j][k] = f[k].f(aPairs.getAPair(i,j),beta);
                    scfPotential = (PotentialPolarizable) f[k].getPotential();
                    // Here is how I get the individual atoms to pass to the energySCF method; kmb 7/10/06
                    // Check how this all works with a pairwise potential first!
                    // Do I grab the right atoms?
                    // Do they yield a value of 0 for energySCF and deltaC below?; kmb 7/10/06
                    
                    // NEED k==0 IN BOOLEAN FOR B4 CALC, TO PREVENT OVERWRITING THE ATOM LISTS FOR THE
                    // TRIMER AND TETRAMER SCF ENERGY CALCS; KMB, 8/16/06
                    if (nPoints > 2 && i==0 && j==1 && k==0) {
                    	scfAtoms.add(aPairs.getAPair(i,j).getAtom(0));
                    	scfAtoms.add(aPairs.getAPair(i,j).getAtom(1));
                        	
//                        	O1.E(((IAtomGroup)aPairs.getAPair(i,j).getAtom(0)).getChildList()node.l);
  //                              O2.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());

				// THESE LINES ADDED FOR B4 CASE; KMB, 8/15/06
                    	scfAtoms123.add(aPairs.getAPair(i,j).getAtom(0));
                    	scfAtoms123.add(aPairs.getAPair(i,j).getAtom(1));
                    	
                    	scfAtoms124.add(aPairs.getAPair(i,j).getAtom(0));
                    	scfAtoms124.add(aPairs.getAPair(i,j).getAtom(1));
                    	
                    	scfAtoms134.add(aPairs.getAPair(i,j).getAtom(0));
                    	
                    	scfAtoms234.add(aPairs.getAPair(i,j).getAtom(1));
                    	
                    	u12Pol = scfPotential.getPolarizationEnergy();

                    }
                    if (i==0 && j==2 && k==0) {
                    	//thirdAtom = aPairs.getAPair(i,j).getAtom(1);  // Does this ensure the third Atom is unique to the other two
                    												// obtained above? kmb, 7/11/06
                    	scfAtoms.add(aPairs.getAPair(i,j).getAtom(1));
/*                        	u13 = (-1)/beta*Math.log(fValues[0][2][0]+1);
                    		//u13Coulombic = f[k].energyCoul(aPairs.getAPair(i,j));
                    		O1.E(aPairs.getAPair(i,j).getAtom(0).node.lastLeafAtom().coord.position());
                    		O3.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());
                    		double rO1O3 = Math.sqrt(O1.Mv1Squared(O3));
                    		if (rO1O3 > 100) {
//                    			System.out.println("O1-O3 distance = " + rO1O3 + ", u13 = " + u13);
                    		}
*/

//                        	O3.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());

				// THESE LINES ADDED FOR B4 CASE; KMB, 8/15/06
                    	scfAtoms123.add(aPairs.getAPair(i,j).getAtom(1));
                    	scfAtoms134.add(aPairs.getAPair(i,j).getAtom(1));
                    	scfAtoms234.add(aPairs.getAPair(i,j).getAtom(1));
                    	u13Pol = scfPotential.getPolarizationEnergy();
                    }
                    if (i==1 && j==2 && k==0) {
//                        	Atom thirdAtom = aPairs.getAPair(i,j).getAtom(0);
/*                      	u23 = (-1)/beta*Math.log(fValues[1][2][0]+1);
                      	//u23Coulombic = f[k].energyCoul(aPairs.getAPair(i,j));
                      	O2.E(aPairs.getAPair(i,j).getAtom(0).node.lastLeafAtom().coord.position());
                      	O3.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());
                      	double rO2O3 = Math.sqrt(O2.Mv1Squared(O3));
                      	if (rO2O3 > 100) {
    //                  		System.out.println("O2-O3 distance = " + rO2O3 + ", u23 = " + u23);
                      	}
*/
                    	u23Pol = scfPotential.getPolarizationEnergy();
                    }
                    
                    if (i==0 && j==3 && k==0) {
//                        	Atom thirdAtom = aPairs.getAPair(i,j).getAtom(0);
/*                      	u23 = (-1)/beta*Math.log(fValues[1][2][0]+1);
                      	//u23Coulombic = f[k].energyCoul(aPairs.getAPair(i,j));
                      	O2.E(aPairs.getAPair(i,j).getAtom(0).node.lastLeafAtom().coord.position());
                      	O3.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());
                      	double rO2O3 = Math.sqrt(O2.Mv1Squared(O3));
                      	if (rO2O3 > 100) {
    //                  		System.out.println("O2-O3 distance = " + rO2O3 + ", u23 = " + u23);
                      	}
*/
                        	
//				O4.E(aPairs.getAPair(i,j).getAtom(1).node.lastLeafAtom().coord.position());

                        scfAtoms124.add(aPairs.getAPair(i,j).getAtom(1));
                    	scfAtoms134.add(aPairs.getAPair(i,j).getAtom(1));
                    	scfAtoms234.add(aPairs.getAPair(i,j).getAtom(1));
                    	scfAtoms.add(aPairs.getAPair(i,j).getAtom(1));
                    }

                    fValues[j][i][k] = fValues[i][j][k];
                }
            }
        }

// Need lines below for B3, but not for B2; kmb, 7/27/06
        if (nPoints > 2) {
	        u12 = (-1)/beta*Math.log(fValues[1][0][0]+1);
	        u13 = (-1)/beta*Math.log(fValues[2][0][0]+1);
	        u23 = (-1)/beta*Math.log(fValues[2][1][0]+1);
        }
        
        if (nPoints > 3) {
	        // FOR B4; KMB, 8/15/06
	        u14 = (-1)/beta*Math.log(fValues[3][0][0]+1);
	        u24 = (-1)/beta*Math.log(fValues[3][1][0]+1);
	        u34 = (-1)/beta*Math.log(fValues[3][2][0]+1);
        }
        
      	
        value = 0.0;

        for(int i=0; i<clusters.length; i++) {
            //System.out.println("clusters.length = " + clusters.length);
        		// clusters.length = 1 for B3
        		double v = clusters[i].value(fValues);
        		//System.out.println("v = " + v);
            value += clusterWeights[i] * v;
            //System.out.println("value = " + value);
            //System.out.println("clusterWeights[i] = " + clusterWeights[i]);
        }
        

        
        /*
         * HOW DO I WANT TO HANDLE THE DIFFERENT CASES - B3 VS. B4?
         * 
         * CAN I USE THE VARIABLE NPOINTS TO TELL ME WHETHER I HAVE 3 OR 4 MOLECULES?
         * 
         * LEFT OFF HERE; KMB, 8/15/06
         */
        
        
        
        // KMB 7/10/06; Presumably, I can call and add the SCF contribution here, after
        // calculation of the f values and multiplying them to get a number for value.
        		
/*	        	if (nPoints == 2) {
		        		
	    			arrayOODistanceIndex = findOORange2(O1,O2);
	    			
	    			if (arrayOODistanceIndex < 9) {
	            		if (Math.abs(value) > maxOOPiValues[arrayOODistanceIndex]) {
	            			maxOOPiValues[arrayOODistanceIndex] = value;
	            			maxOOdistances[arrayOODistanceIndex] = Math.sqrt(O1.Mv1Squared(O2));
	                 }
	    			}
	        }
*/        	
        	
        	
		if (nPoints == 3) {
    		if (u12+u13+u23 == Double.POSITIVE_INFINITY) {
    			//System.out.println("Sum of pair energies is infinity: u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23);
    			deltaC = 0.0;
    		}
    		else {
    		// Get a handle on the list of atoms from the AtomPairSet
    			
    			u123 = ((Potential)scfPotential).energy(scfAtoms);
    			double u123Pol = scfPotential.getPolarizationEnergy();
    			 //This conditional handles the case when the triplet energy does
    			 //not quite agree with the dimer energy, for the case where the third
    			 //water molecule is 1000 angstroms of farther away (see code in 
    			 //PotentialWaterPPC8forB3.java).
    			 //kmb, 8/3/06
    			 
    			//if (u123 == -123456789.0) {
    			//	u123 = u12 + u13 + u23;
    			//}

    			//u123 += u12 + u13 + u23;  // for testing TIP4P potential only; kmb, 7/26/06	
    		
    			//deltaC = Math.exp(-beta*u123) - Math.exp(-beta*(u12 + u13 + u23));
    			double deltau123 = u123Pol-(u12Pol+u13Pol+u23Pol);
    			deltaC = (Math.exp(-beta*deltau123)-1)*Math.exp(-beta*(u12 + u13 + u23));       
    			
    			// kmb added this line 7/18/06
    			// deltaC has to be multiplied by clusterWeights, just like v was multiplied by
    			// clusterWeights above to get value
    			deltaC = deltaC*clusterWeights[0];
    			
    			if (Double.isInfinite(deltaC)) {
    				System.out.println("deltaC = " + deltaC);
    			}
    		}
	        		
    		//System.out.println("u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23 + ", u123 = " + u123 + ", deltaC = " + deltaC + ", value(before) = " + value);
    		value += deltaC;
	        		
/*        			arrayOODistanceIndex = findOORange3(O1,O2,O3);
        			
//        			System.out.println("Updating pi vs. O-O array for LJ fluid, too!");
        			
        			if (arrayOODistanceIndex < 9) {
	            		if (value > maxOOPiValues[arrayOODistanceIndex]) {
	            			maxOOPiValues[arrayOODistanceIndex] = value;					
	                 }
        			}
*/
		}

		if (nPoints == 4) {
			
			if (u12+u13+u14+u23+u24+u34 == Double.POSITIVE_INFINITY) {
    			//System.out.println("Sum of pair energies is infinity: u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23);
				deltaD = 0.0;
    		}
			else {
				u123 = ((Potential)scfPotential).energy(scfAtoms123);
				u124 = ((Potential)scfPotential).energy(scfAtoms124);
				u134 = ((Potential)scfPotential).energy(scfAtoms134);
				u234 = ((Potential)scfPotential).energy(scfAtoms234);
				u1234 = ((Potential)scfPotential).energy(scfAtoms);

				double deltaU123 = u123-u12-u13-u23;
				double deltaU124 = u124-u12-u14-u24;
				double deltaU134 = u134-u13-u14-u34;
				double deltaU234 = u234-u23-u24-u34;

				double deltaU1234 = u1234-u12-u13-u14-u23-u24-u34-deltaU123-deltaU124-deltaU134-deltaU234;
				
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
					// Old assignment 6/2/07, KMB; REPLACED ABOVE
//					if (u1234 == -123456789.0) {
//                                                deltaD = 0.0;
//                                        }


//        				 kmb added this line 8/16/06
	        			// deltaD has to be multiplied by weightPrefactor from Standard class, just like deltaC was multiplied by
	        			// clusterWeights above to get value; note, for B3 clusterWeights = weightPrefactor
        				
				deltaD = deltaD*(-32.0);  //XXX - IS THIS DIFFERENT NOW?; KMB, 7/27/07  -32 is value of weightPrefactor in Standard for a B4 calc
			}
        			
			value += deltaD;

/*        			arrayOODistanceIndex = findOORange(O1,O2,O3,O4);
        			
        			if (arrayOODistanceIndex < 9) {
	            		if (value > maxOOPiValues[arrayOODistanceIndex]) {
	            			maxOOPiValues[arrayOODistanceIndex] = value;					
	                 }
        			}
*/            		

    	}
        
        if (u12 < Double.POSITIVE_INFINITY && u13 < Double.POSITIVE_INFINITY && u23 < Double.POSITIVE_INFINITY) {
        		//System.out.println("u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23);
    			//System.out.println("u12 = " + u12 + ", u13 = " + u13 + ", u23 = " + u23 + ", u123 = " + u123 + ", deltaC = " + deltaC + ", value(before) = " + value);
        }
        
//        System.out.println("PI (from Cluster Sum) = " + value);
//        	System.out.println(value);
        
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
    
/*    // KMB methods for plotting pi value vs. max O-O distance; 6/11/07
    public int findOORange(Vector3D O1, Vector3D O2, Vector3D O3, Vector3D O4) {
        
		// Adding code to find cutoff for B5 calc; KMB 5/25/07
		double rO1O2 = Math.sqrt(O1.Mv1Squared(O2));
		double rO1O3 = Math.sqrt(O1.Mv1Squared(O3));
         double rO1O4 = Math.sqrt(O1.Mv1Squared(O4));
		double rO2O3 = Math.sqrt(O2.Mv1Squared(O3));
		double rO2O4 = Math.sqrt(O2.Mv1Squared(O4));
		double rO3O4 = Math.sqrt(O3.Mv1Squared(O4));
		double maxDistance;
		
		if (rO1O2 > rO1O3) {
			maxDistance = rO1O2;					
         }
		else {
			maxDistance = rO1O3;
		}
		
		if (rO1O4 > maxDistance) {
			maxDistance = rO1O4;					
         }

		if (rO2O3 > maxDistance) {
			maxDistance = rO2O3;					
         }
		
		if (rO2O4 > maxDistance) {
			maxDistance = rO2O4;					
         }

		if (rO3O4 > maxDistance) {
			maxDistance = rO3O4;					
         }
		
//		System.out.println("Max O-O distance = " + maxDistance);
		
		double logMaxDistance = Math.log(maxDistance);
		
		if (logMaxDistance <= 1.2) {
			//bin1Count++;
			binCounts[0]++;
			return 0;					
         }
		
		if (logMaxDistance > 1.2 && logMaxDistance <= 1.4) {
			//bin2Count++;
			binCounts[1]++;
			return 1;					
         }
		
		if (logMaxDistance > 1.4 && logMaxDistance <= 1.6) {
			//bin3Count++;
			binCounts[2]++;
			return 2;					
         }
		
		if (logMaxDistance > 1.6 && logMaxDistance <= 1.8) {
			//bin4Count++;
			binCounts[3]++;
			return 3;					
         }
		
		if (logMaxDistance > 1.8 && logMaxDistance <= 2.0) {
			//bin5Count++;
			binCounts[4]++;
			return 4;					
         }
		
		if (logMaxDistance > 2.0 && logMaxDistance <= 2.2) {
			//bin6Count++;
			binCounts[5]++;
			return 5;					
         }

		if (logMaxDistance > 2.2 && logMaxDistance <= 2.4) {
			//bin7Count++;
			binCounts[6]++;
			return 6;					
         }
		
		if (logMaxDistance > 2.4 && logMaxDistance <= 2.6) {
			//bin8Count++;
			binCounts[7]++;
			return 7;					
         }
		
		if (logMaxDistance > 2.6 && logMaxDistance <= 2.8) {
			//bin9Count++;
			binCounts[8]++;
			return 8;					
         }
		else {
			return 9;
		}
		
		if (dist12 || dist13 || dist14 || dist23 || dist24 || dist34) {
			System.out.println(rO1O2 + " " + rO1O3 + " " + rO1O4 + " " + " " + rO2O3 + " " + rO2O4 + " " + " " + rO3O4 + " " + value);
//			System.out.println(dist12 + " " + dist13 + " " + dist14 + " " + dist15 + " " + dist23 + " " + dist24 + " " + dist25 + " " + dist34 + " " + dist35 + " " + dist45);
		}

    }
    
    public void updateMaxValue(int rangeIndex, int piValue, double[] maxValues) {
        
    		if (piValue > maxValues[rangeIndex]) {
			maxValues[rangeIndex] = piValue;					
         }
    	
    }
    
    public double[] getMaxValueArray() {
        
    		return maxOOPiValues;    	
    }

    public double[] getBinCountArray() {
        
    		return binCounts;    	
    }
    
    
//  KMB methods for plotting pi value vs. max O-O distance; 6/11/07
    public int findOORange3(Vector3D O1, Vector3D O2, Vector3D O3) {
        
		// Adding code to find cutoff for B5 calc; KMB 5/25/07
		double rO1O2 = Math.sqrt(O1.Mv1Squared(O2));
		double rO1O3 = Math.sqrt(O1.Mv1Squared(O3));
		double rO2O3 = Math.sqrt(O2.Mv1Squared(O3));
		double maxDistance;
		
		if (rO1O2 > rO1O3) {
			maxDistance = rO1O2;					
         }
		else {
			maxDistance = rO1O3;
		}
		
		if (rO2O3 > maxDistance) {
			maxDistance = rO2O3;					
         }
		
	
//		System.out.println("Max O-O distance = " + maxDistance);
		
		double logMaxDistance = Math.log(maxDistance);
		
		if (logMaxDistance <= 2.5) {
			//bin1Count++;
			binCounts[0]++;
			return 0;					
         }
		
		if (logMaxDistance > 2.5 && logMaxDistance <= 3.0) {
			//bin2Count++;
			binCounts[1]++;
			return 1;					
         }
		
		if (logMaxDistance > 3.0 && logMaxDistance <= 3.5) {
			//bin3Count++;
			binCounts[2]++;
			return 2;					
         }
		
		if (logMaxDistance > 3.5 && logMaxDistance <= 4.0) {
			//bin4Count++;
			binCounts[3]++;
			return 3;					
         }
		
		if (logMaxDistance > 4.0 && logMaxDistance <= 4.5) {
			//bin5Count++;
			binCounts[4]++;
			return 4;					
         }
		
		if (logMaxDistance > 4.5 && logMaxDistance <= 5.0) {
			//bin6Count++;
			binCounts[5]++;
			return 5;					
         }

		if (logMaxDistance > 5.0 && logMaxDistance <= 5.5) {
			//bin7Count++;
			binCounts[6]++;
			return 6;					
         }
		
		if (logMaxDistance > 5.5 && logMaxDistance <= 6.0) {
			//bin8Count++;
			binCounts[7]++;
			return 7;					
         }
		
		if (logMaxDistance > 6.0 && logMaxDistance <= 6.5) {
			//bin9Count++;
			binCounts[8]++;
			return 8;					
         }
		else {
			return 9;
		}
		
		if (dist12 || dist13 || dist14 || dist23 || dist24 || dist34) {
			System.out.println(rO1O2 + " " + rO1O3 + " " + rO1O4 + " " + " " + rO2O3 + " " + rO2O4 + " " + " " + rO3O4 + " " + value);
//			System.out.println(dist12 + " " + dist13 + " " + dist14 + " " + dist15 + " " + dist23 + " " + dist24 + " " + dist25 + " " + dist34 + " " + dist35 + " " + dist45);
		}

    }
    
//  KMB methods for plotting pi value vs. max O-O distance; 6/11/07
    public int findOORange2(Vector3D O1, Vector3D O2) {
        
		// Adding code to find cutoff for B5 calc; KMB 5/25/07
		double rO1O2 = Math.sqrt(O1.Mv1Squared(O2));
		double maxDistance;
		
			maxDistance = rO1O2;					
		
		
//		System.out.println("Max O-O distance = " + maxDistance);
		
		double logMaxDistance = Math.log(maxDistance);
		
		if (logMaxDistance <= 2.5) {
			//bin1Count++;
			binCounts[0]++;
			return 0;					
         }
		
		if (logMaxDistance > 2.5 && logMaxDistance <= 3.0) {
			//bin2Count++;
			binCounts[1]++;
			return 1;					
         }
		
		if (logMaxDistance > 3.0 && logMaxDistance <= 3.5) {
			//bin3Count++;
			binCounts[2]++;
			return 2;					
         }
		
		if (logMaxDistance > 3.5 && logMaxDistance <= 4.0) {
			//bin4Count++;
			binCounts[3]++;
			return 3;					
         }
		
		if (logMaxDistance > 4.0 && logMaxDistance <= 4.5) {
			//bin5Count++;
			binCounts[4]++;
			return 4;					
         }
		
		if (logMaxDistance > 4.5 && logMaxDistance <= 5.0) {
			//bin6Count++;
			binCounts[5]++;
			return 5;					
         }

		if (logMaxDistance > 5.0 && logMaxDistance <= 5.5) {
			//bin7Count++;
			binCounts[6]++;
			return 6;					
         }
		
		if (logMaxDistance > 5.5 && logMaxDistance <= 6.0) {
			//bin8Count++;
			binCounts[7]++;
			return 7;					
         }
		
		if (logMaxDistance > 6.0 && logMaxDistance <= 6.5) {
			//bin9Count++;
			binCounts[8]++;
			return 8;					
         }
		else {
			return 9;
		}
		
		if (dist12 || dist13 || dist14 || dist23 || dist24 || dist34) {
			System.out.println(rO1O2 + " " + rO1O3 + " " + rO1O4 + " " + " " + rO2O3 + " " + rO2O4 + " " + " " + rO3O4 + " " + value);
//			System.out.println(dist12 + " " + dist13 + " " + dist14 + " " + dist15 + " " + dist23 + " " + dist24 + " " + dist25 + " " + dist34 + " " + dist35 + " " + dist45);
		}

    }
*/    
    private static final long serialVersionUID = 1L;
    private final ClusterBonds[] clusters;
    private final double[] clusterWeights;
    private final MayerFunction[] f;
    private final double[][][] fValues;
    private int cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
    private double beta;
    private double[] maxOOPiValues, binCounts, maxOOdistances; // final or not?
    protected final AtomArrayList scfAtoms;
    protected final AtomArrayList scfAtoms123;
    protected final AtomArrayList scfAtoms124;
    protected final AtomArrayList scfAtoms134;
    protected final AtomArrayList scfAtoms234;
}
