/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
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

        scfAtoms = new MoleculeArrayList(5); // USE THIS LIST FOR ALL ATOMS, WHETHER 3 OR 4; KMB, 8/16/06
    }

    // equal point count enforced in constructor 
    public int pointCount() {
        return clusters[0].pointCount();
    }
    
    public ClusterAbstract makeCopy() {
        ClusterSumPolarizable copy = new ClusterSumPolarizable(clusters,clusterWeights,f);
        copy.setTemperature(1/beta);
        copy.setDeltaCut(Math.sqrt(deltaCut2));
        copy.setCaching(doCaching);
        return copy;
    }
    
    public void setCaching(boolean doCaching) {
        this.doCaching = doCaching;
    }

    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        AtomPairSet aPairs = box.getAPairSet();
        IMoleculeList atomSet = box.getMoleculeList();
        long thisCPairID = cPairs.getID();
//        System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (doCaching) {
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
        }
        
        for (int k=0; k<f.length; k++) {
            f[k].setBox(box);
        }

        int nPoints = pointCount();
        
        // recalculate all f values for all pairs
        PotentialPolarizable scfPotential = null;
        // deltaD and deltaE run into precision problems for long distances
        boolean truncateDelta = false;
        boolean tooClose = pushR2 > 0;
        for(int i=0; i<nPoints-1; i++) {
            for(int j=i+1; j<nPoints; j++) {
                double r2 = cPairs.getr2(i,j);
                if (r2 > pushR2) {
                    tooClose = false;
                }
                else if (r2 > deltaCut2) {
                    truncateDelta = true;
                    break;
                }
                for(int k=0; k<f.length; k++) {
                    fValues[i][j][k] = f[k].f(aPairs.getAPair(i,j),cPairs.getr2(i,j), beta);
                    fValues[j][i][k] = fValues[i][j][k];
                    scfPotential = (PotentialPolarizable) f[k].getPotential();
                    if (k==0) uijPol[i][j] = scfPotential.getLastPolarizationEnergy();
                }
            }
        }

        value = 0.0;
        if (tooClose) {
            return value;
        }

        for(int i=0; i<clusters.length; i++) {
            //System.out.println("clusters.length = " + clusters.length);
            // clusters.length = 1 for B3
            double v = clusters[i].value(fValues);
            value += clusterWeights[i] * v;
            //System.out.println("value = " + value);
            //System.out.println("clusterWeights["+i+"] = " + clusterWeights[i]);
        }
        
        if (truncateDelta) {
            return value;
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
                scfAtoms.add(atomSet.get(0));
                scfAtoms.add(atomSet.get(1));
                scfAtoms.add(atomSet.get(2));
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
            scfAtoms.add(atomSet.get(0));
            scfAtoms.add(atomSet.get(1));
            scfAtoms.add(atomSet.get(2));

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
            scfAtoms.add(atomSet.get(3));
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
            scfAtoms.add(atomSet.get(2));
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
            scfAtoms.add(atomSet.get(1));
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

            scfAtoms.add(atomSet.get(0));
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
//                  + 2*Math.exp(-beta*(u13+u14+u34))*(Math.exp(-beta*deltaU134)-1) + 2*Math.exp(-beta*(u23+u24+u34))*(Math.exp(-beta*deltaU234)-1)
//                  + Math.exp(-beta*(u12+u13+u23-u14))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u13+u23-u24))*(1-Math.exp(-beta*deltaU123))
//                  + Math.exp(-beta*(u12+u13+u23-u34))*(1-Math.exp(-beta*deltaU123)) + Math.exp(-beta*(u12+u14+u24-u13))*(1-Math.exp(-beta*deltaU124))
//                  + Math.exp(-beta*(u12+u14+u24-u23))*(1-Math.exp(-beta*deltaU124)) + Math.exp(-beta*(u12+u14+u24-u34))*(1-Math.exp(-beta*deltaU124))
//                  + Math.exp(-beta*(u13+u14+u34-u12))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u13+u14+u34-u23))*(1-Math.exp(-beta*deltaU134))
//                  + Math.exp(-beta*(u13+u14+u34-u24))*(1-Math.exp(-beta*deltaU134)) + Math.exp(-beta*(u23+u24+u34-u12))*(1-Math.exp(-beta*deltaU234))
//                  + Math.exp(-beta*(u23+u24+u34-u13))*(1-Math.exp(-beta*deltaU234)) + Math.exp(-beta*(u23+u24+u34-u14))*(1-Math.exp(-beta*deltaU234))
//                  + Math.exp(-beta*(u12+u13+u14+u23+u24+u34))*(Math.exp(-beta*deltaU1234)-1);


//                       kmb added this line 8/16/06
            // deltaD has to be multiplied by weightPrefactor from Standard class, just like deltaC was multiplied by
            // clusterWeights above to get value; note, for B3 clusterWeights = weightPrefactor

            // -(1/8) is the B4 prefactor multiplying all diagrams.
            // coefficient is -(1-4)/4! = -1/8
            deltaD = -0.125*deltaD;
            value += deltaD;
        }
        else if (nPoints == 5) {
            final double f12 = fValues[0][1][0];
            final double f13 = fValues[0][2][0];
            final double f14 = fValues[0][3][0];
            final double f15 = fValues[0][4][0];
            final double f23 = fValues[1][2][0];
            final double f24 = fValues[1][3][0];
            final double f25 = fValues[1][4][0];
            final double f34 = fValues[2][3][0];
            final double f35 = fValues[2][4][0];
            final double f45 = fValues[3][4][0];
            final double g12 = f12+1; //Math.exp(-beta*u12);
            final double g13 = f13+1; //Math.exp(-beta*u13);
            final double g14 = f14+1; //Math.exp(-beta*u14);
            final double g15 = f15+1; //Math.exp(-beta*u14);
            final double g23 = f23+1; //Math.exp(-beta*u23);
            final double g24 = f24+1; //Math.exp(-beta*u24);
            final double g25 = f25+1; //Math.exp(-beta*u14);
            final double g34 = f34+1; //Math.exp(-beta*u34);
            final double g35 = f35+1; //Math.exp(-beta*u14);
            final double g45 = f45+1; //Math.exp(-beta*u14);
            
            double beta123 = Double.NaN;
            double beta124 = Double.NaN;
            double beta125 = Double.NaN;
            double beta134 = Double.NaN;
            double beta135 = Double.NaN;
            double beta145 = Double.NaN;
            double beta234 = Double.NaN;
            double beta235 = Double.NaN;
            double beta245 = Double.NaN;
            double beta345 = Double.NaN;
            
            double deltaE = 0;

            scfAtoms.clear();
            // we need to properly construct these lists even if we don't use them
            // (due to overlaps) because the next list is obtained by removing/adding
            // atoms from this one.
            scfAtoms.add(atomSet.get(0));
            scfAtoms.add(atomSet.get(1));
            scfAtoms.add(atomSet.get(2));  // 123

            if (g12*g13*g23 != 0) {
                double u123Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU123 = u123Pol-uijPol[0][1]-uijPol[0][2]-uijPol[1][2];
                beta123 = deltaU123*beta;
                double exp123 = -beta123;
                if (Math.abs(beta123) > 1E-8) {
                    exp123 = Math.exp(-beta123) - 1;
                }

                deltaE += exp123*g12*g13*g23*(2*(f14 + f15 + f24 + f25 + f34 + f35 + 2*f45 + 2)
                        + f14*f25 + f14*f35 + f14*f45 + f15*f24 + f15*f34 + f15*f45
                        + f24*f35 + f24*f45 + f25*f34 + f25*f45 + f34*f45 + f35*f45
                        + 2*f14*f15 + 2*f24*f25 + 2*f34*f35);
                
                

//                deltaE += exp123*g12*g13*g23*(2*(g14*g15 + g24*g25 + g34*g35 - g45));
//                        + f45*(g24 + g14 + g15 + g25 + g34 + g35)
//                        + f14*f25 + f14*f35 + f15*f24 + f15*f34 + f24*f35 + f25*f34 
//                        + 2*f14*f15 + 2*f24*f25 + 2*f34*f35);
            }
            
            scfAtoms.remove(2);
            scfAtoms.add(atomSet.get(3));  // 124

            if (g12*g14*g24 != 0) {
                double u124Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU124 = u124Pol-uijPol[0][1]-uijPol[0][3]-uijPol[1][3];
                beta124 = deltaU124*beta;
                double exp124 = -beta124;
                if (Math.abs(beta124) > 1E-8) {
                    exp124 = Math.exp(-beta124) - 1;
                }

                deltaE += exp124*g12*g14*g24*(2*(f13 + f15 + f23 + f25 + f34 + f45 + 2*f35 + 2)
                        + f13*f25 + f13*f35 + f13*f45 + f15*f23 + f15*f34 + f15*f35
                        + f23*f35 + f23*f45 + f25*f34 + f25*f35 + f34*f35 + f35*f45
                        + 2*f13*f15 + 2*f23*f25 + 2*f34*f45);
            }

            scfAtoms.remove(2);
            scfAtoms.add(atomSet.get(4));  // 125

            if (g12*g15*g25 != 0) {
                double u125Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU125 = u125Pol-uijPol[0][1]-uijPol[0][4]-uijPol[1][4];
                beta125 = deltaU125*beta;
                double exp125 = -beta125;
                if (Math.abs(beta125) > 1E-8) {
                    exp125 = Math.exp(-beta125) - 1;
                }

                deltaE += exp125*g12*g15*g25*(2*(f13 + f14 + f23 + f24 + 2*f34 + f45 + f35 + 2)
                        + f13*f24 + f13*f34 + f13*f45 + f14*f23 + f14*f34 + f14*f35
                        + f23*f34 + f23*f45 + f24*f34 + f24*f35 + f35*f34 + f45*f34
                        + 2*f13*f14 + 2*f23*f24 + 2*f35*f45);
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(2));  // 153
            
            if (g13*g15*g35 != 0) {
                double u135Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU135 = u135Pol-uijPol[0][2]-uijPol[0][4]-uijPol[2][4];
                beta135 = deltaU135*beta;
                double exp135 = -beta135;
                if (Math.abs(beta135) > 1E-8) {
                    exp135 = Math.exp(-beta135) - 1;
                }

                deltaE += exp135*g13*g15*g35*(2*(f12 + f14 + f23 + f25 + f34 + f45 + 2*f24 + 2)
                        + f12*f24 + f12*f34 + f12*f45 + f14*f23 + f14*f25 + f14*f24
                        + f23*f24 + f23*f45 + f25*f24 + f25*f34 + f24*f34 + f24*f45
                        + 2*f12*f14 + 2*f23*f34 + 2*f25*f45);
                
                if (!Double.isNaN(beta124)) {
                    double exp124_135 = beta124+beta135;
                    if (Math.abs(exp124_135) > 1E-8) {
                        exp124_135 = 1 - Math.exp(-exp124_135);
                    }
                    deltaE += exp124_135*g12*g14*g24*g13*g15*g35;
                }
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(3));  // 134
            
            if (g13*g14*g34 != 0) {
                double u134Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU134 = u134Pol-uijPol[0][2]-uijPol[0][3]-uijPol[2][3];
                beta134 = deltaU134*beta;
                double exp134 = -beta134;
                if (Math.abs(beta134) > 1E-8) {
                    exp134 = Math.exp(-beta134) - 1;
                }

                deltaE += exp134*g13*g14*g34*(2*(f12 + f15 + f23 + f24 + f35 + f45 + 2*f25 + 2)
                        + f12*f25 + f12*f35 + f12*f45 + f15*f23 + f15*f24 + f15*f25
                        + f23*f25 + f23*f45 + f24*f25 + f24*f35 + f25*f35 + f25*f45
                        + 2*f12*f15 + 2*f23*f35 + 2*f24*f45);

                if (!Double.isNaN(beta125)) {
                    double exp125_134 = beta125+beta134;
                    if (Math.abs(exp125_134) > 1E-8) {
                        exp125_134 = 1 - Math.exp(-exp125_134);
                    }
                    deltaE += exp125_134*g12*g15*g25*g13*g14*g34;
                }
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(4));  // 145

            if (g14*g15*g45 != 0) {
                double u145Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU145 = u145Pol-uijPol[0][4]-uijPol[0][3]-uijPol[3][4];
                beta145 = deltaU145*beta;
                double exp145 = -beta145;
                if (Math.abs(beta145) > 1E-8) {
                    exp145 = Math.exp(-beta145) - 1;
                }

                deltaE += exp145*g14*g15*g45*(2*(f12 + f13 + f24 + f25 + f34 + f35 + 2*f23 + 2)
                        + f12*f23 + f12*f34 + f12*f35 + f13*f24 + f13*f25 + f13*f23
                        + f23*f24 + f24*f35 + f25*f23 + f25*f34 + f23*f34 + f23*f35
                        + 2*f12*f13 + 2*f24*f34 + 2*f25*f35);

                if (!Double.isNaN(beta123)) {
                    double exp123_145 = beta123+beta145;
                    if (Math.abs(exp123_145) > 1E-8) {
                        exp123_145 = 1 - Math.exp(-exp123_145);
                    }
                    deltaE += exp123_145*g12*g13*g23*g14*g15*g45;
                }
            }

            scfAtoms.remove(0);
            scfAtoms.add(atomSet.get(1));  // 452
            
            if (g24*g25*g45 != 0) {
                double u245Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU245 = u245Pol-uijPol[1][4]-uijPol[1][3]-uijPol[3][4];
                beta245 = deltaU245*beta;
                double exp245 = -beta245;
                if (Math.abs(beta245) > 1E-8) {
                    exp245 = Math.exp(-beta245) - 1;
                }
                deltaE += exp245*g24*g25*g45*(2*(f12 + f14 + f15 + f23 + f34 + f35 + 2*f13 + 2)
                        + f12*f13 + f12*f34 + f12*f35 + f14*f23 + f14*f13 + f14*f35
                        + f15*f13 + f15*f23 + f15*f34 + f13*f23 + f13*f34 + f13*f35
                        + 2*f12*f23 + 2*f14*f34 + 2*f15*f35);

                if (!Double.isNaN(beta123)) {
                    double exp123_245 = beta123+beta245;
                    if (Math.abs(exp123_245) > 1E-8) {
                        exp123_245 = 1 - Math.exp(-exp123_245);
                    }
                    deltaE += exp123_245*g12*g13*g23*g24*g25*g45;
                }
                if (!Double.isNaN(beta134)) {
                    double exp134_245 = beta134+beta245;
                    if (Math.abs(exp134_245) > 1E-8) {
                        exp134_245 = 1 - Math.exp(-exp134_245);
                    }
                    deltaE += exp134_245*g13*g14*g34*g24*g25*g45;
                }
                if (!Double.isNaN(beta135)) {
                    double exp135_245 = beta135+beta245;
                    if (Math.abs(exp135_245) > 1E-8) {
                        exp135_245 = 1 - Math.exp(-exp135_245);
                    }
                    deltaE += exp135_245*g13*g15*g35*g24*g25*g45;
                }
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(2));  // 423
            
            if (g23*g24*g34 != 0) {
                double u234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU234 = u234Pol-uijPol[1][2]-uijPol[1][3]-uijPol[2][3];
                beta234 = deltaU234*beta;
                double exp234 = -beta234;
                if (Math.abs(beta234) > 1E-8) {
                    exp234 = Math.exp(-beta234) - 1;
                }
                deltaE += exp234*g23*g24*g34*(2*(f12 + f13 + f14 + f25 + f35 + f45 + 2*f15 + 2)
                        + f12*f15 + f12*f35 + f12*f45 + f13*f25 + f13*f15 + f13*f45
                        + f14*f15 + f14*f25 + f14*f35 + f15*f25 + f15*f35 + f15*f45
                        + 2*f12*f25 + 2*f13*f35 + 2*f14*f45);

                if (!Double.isNaN(beta125)) {
                    double exp125_234 = beta125+beta234;
                    if (Math.abs(exp125_234) > 1E-8) {
                        exp125_234 = 1 - Math.exp(-exp125_234);
                    }
                    deltaE += exp125_234*g12*g15*g25*g23*g24*g34;
                }
                if (!Double.isNaN(beta135)) {
                    double exp135_234 = beta135+beta234;
                    if (Math.abs(exp135_234) > 1E-8) {
                        exp135_234 = 1 - Math.exp(-exp135_234);
                    }
                    deltaE += exp135_234*g13*g15*g35*g23*g24*g34;
                }
                if (!Double.isNaN(beta145)) {
                    double exp145_234 = beta145+beta234;
                    if (Math.abs(exp145_234) > 1E-8) {
                        exp145_234 = 1 - Math.exp(-exp145_234);
                    }
                    deltaE += exp145_234*g14*g15*g45*g23*g24*g34;
                }
            }

            scfAtoms.remove(0);
            scfAtoms.add(atomSet.get(4));  // 235
            
            if (g23*g25*g35 != 0) {
                double u235Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU235 = u235Pol-uijPol[1][2]-uijPol[1][4]-uijPol[2][4];
                beta235 = deltaU235*beta;
                double exp235 = -beta235;
                if (Math.abs(beta235) > 1E-8) {
                    exp235 = Math.exp(-beta235) - 1;
                }
                deltaE += exp235*g23*g25*g35*(2*(f12 + f13 + f15 + f24 + f34 + f45 + 2*f14 + 2)
                        + f12*f14 + f12*f34 + f12*f45 + f13*f24 + f13*f14 + f13*f45
                        + f15*f14 + f15*f24 + f15*f34 + f14*f24 + f14*f34 + f14*f45
                        + 2*f12*f24 + 2*f13*f34 + 2*f15*f45);

                if (!Double.isNaN(beta124)) {
                    double exp124_235 = beta124+beta235;
                    if (Math.abs(exp124_235) > 1E-8) {
                        exp124_235 = 1 - Math.exp(-exp124_235);
                    }
                    deltaE += exp124_235*g12*g14*g24*g23*g25*g35;
                }
                if (!Double.isNaN(beta134)) {
                    double exp134_235 = beta134+beta235;
                    if (Math.abs(exp134_235) > 1E-8) {
                        exp134_235 = 1 - Math.exp(-exp134_235);
                    }
                    deltaE += exp134_235*g13*g14*g34*g23*g25*g35;
                }
                if (!Double.isNaN(beta145)) {
                    double exp145_235 = beta145+beta235;
                    if (Math.abs(exp145_235) > 1E-8) {
                        exp145_235 = 1 - Math.exp(-exp145_235);
                    }
                    deltaE += exp145_235*g14*g15*g45*g23*g25*g35;
                }
            }

            scfAtoms.remove(0);
            scfAtoms.add(atomSet.get(3));  // 354
            
            if (g34*g35*g45 != 0) {
                double u345Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU345 = u345Pol-uijPol[2][3]-uijPol[2][4]-uijPol[3][4];
                beta345 = deltaU345*beta;
                double exp345 = -beta345;
                if (Math.abs(beta345) > 1E-8) {
                    exp345 = Math.exp(-beta345) - 1;
                }
                deltaE += exp345*g34*g35*g45*(2*(f13 + f14 + f15 + f23 + f24 + f25 + 2*f12 + 2)
                        + f13*f12 + f13*f24 + f13*f25 + f14*f23 + f14*f12 + f14*f25
                        + f15*f12 + f15*f23 + f15*f24 + f12*f23 + f12*f24 + f12*f25
                        + 2*f13*f23 + 2*f14*f24 + 2*f15*f25);
                
                if (!Double.isNaN(beta123)) {
                    double exp123_345 = beta123+beta345;
                    if (Math.abs(exp123_345) > 1E-8) {
                        exp123_345 = 1 - Math.exp(-exp123_345);
                    }
                    deltaE += exp123_345*g12*g13*g23*g34*g35*g45;
                }
                if (!Double.isNaN(beta124)) {
                    double exp124_345 = beta124+beta345;
                    if (Math.abs(exp124_345) > 1E-8) {
                        exp124_345 = 1 - Math.exp(-exp124_345);
                    }
                    deltaE += exp124_345*g12*g14*g24*g34*g35*g45;
                }
                if (!Double.isNaN(beta125)) {
                    double exp125_345 = beta125+beta345;
                    if (Math.abs(exp125_345) > 1E-8) {
                        exp125_345 = 1 - Math.exp(-exp125_345);
                    }
                    deltaE += exp125_345*g12*g15*g25*g34*g35*g45;
                }
            }

            scfAtoms.add(atomSet.get(0));  // 3541
            if (g13*g14*g15*g34*g35*g45 != 0) {
                double u1345Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU1345 = u1345Pol-uijPol[0][2]-uijPol[0][3]-uijPol[0][4]-uijPol[2][3]-uijPol[2][4]-uijPol[3][4];
                double beta1345 = deltaU1345*beta;
                double exp1345 = -beta1345;
                if (Math.abs(beta1345) > 1E-8) {
                    exp1345 = Math.exp(-beta1345) - 1;
                }
                deltaE += -exp1345*g13*g14*g15*g34*g35*g45*(f12 + f23 + f24 + f25 + 1);
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(1));  // 3412
            if (g12*g13*g14*g23*g24*g34 != 0) {
                double u1234Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU1234 = u1234Pol-uijPol[0][1]-uijPol[0][2]-uijPol[0][3]-uijPol[1][2]-uijPol[1][3]-uijPol[2][3];
                double beta1234 = deltaU1234*beta;
                double exp1234 = -beta1234;
                if (Math.abs(beta1234) > 1E-8) {
                    exp1234 = Math.exp(-beta1234) - 1;
                }
                deltaE += -exp1234*g12*g13*g14*g23*g24*g34*(f15 + f25 + f35 + f45 + 1);
            }

            scfAtoms.remove(1);
            scfAtoms.add(atomSet.get(4));  // 3125
            if (g12*g13*g15*g23*g25*g35 != 0) {
                double u1235Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU1235 = u1235Pol-uijPol[0][1]-uijPol[0][2]-uijPol[0][4]-uijPol[1][2]-uijPol[1][4]-uijPol[2][4];
                double beta1235 = deltaU1235*beta;
                double exp1235 = -beta1235;
                if (Math.abs(beta1235) > 1E-8) {
                    exp1235 = Math.exp(-beta1235) - 1;
                }
                deltaE += -exp1235*g12*g13*g15*g23*g25*g35*(f14 + f24 + f34 + f45 + 1);
            }

            scfAtoms.remove(0);
            scfAtoms.add(atomSet.get(3));  // 1254
            if (g12*g14*g15*g24*g25*g45 != 0) {
                double u1245Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU1245 = u1245Pol-uijPol[0][1]-uijPol[0][3]-uijPol[0][4]-uijPol[1][3]-uijPol[1][4]-uijPol[3][4];
                double beta1245 = deltaU1245*beta;
                double exp1245 = -beta1245;
                if (Math.abs(beta1245) > 1E-8) {
                    exp1245 = Math.exp(-beta1245) - 1;
                }
                deltaE += -exp1245*g12*g14*g15*g24*g25*g45*(f13 + f23 + f34 + f35 + 1);
            }

            scfAtoms.remove(0);
            scfAtoms.add(atomSet.get(2));  // 2543
            if (g23*g24*g25*g34*g35*g45 != 0) {
                double u2345Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU2345 = u2345Pol-uijPol[1][2]-uijPol[1][3]-uijPol[1][4]-uijPol[2][3]-uijPol[2][4]-uijPol[3][4];
                double beta2345 = deltaU2345*beta;
                double exp2345 = -beta2345;
                if (Math.abs(beta2345) > 1E-8) {
                    exp2345 = Math.exp(-beta2345) - 1;
                }
                deltaE += -exp2345*g23*g24*g25*g34*g35*g45*(f12 + f13 + f14 + f15 + 1);
            }

            scfAtoms.add(atomSet.get(0));  // 25431
            if (g12*g13*g14*g15*g23*g24*g25*g34*g35*g45 != 0) {
                double u12345Pol = scfPotential.getPolarizationEnergy(scfAtoms);
                double deltaU12345 = u12345Pol-uijPol[0][1]-uijPol[0][2]-uijPol[0][3]-uijPol[0][4]-uijPol[1][2]-uijPol[1][3]-uijPol[1][4]
                                              -uijPol[2][3]-uijPol[2][4]-uijPol[3][4];
                double beta12345 = deltaU12345*beta;
                double exp12345 = -beta12345;
                if (Math.abs(beta12345) > 1E-8) {
                    exp12345 = Math.exp(-beta12345) - 1;
                }
                deltaE += exp12345*g12*g13*g14*g15*g23*g24*g25*g34*g35*g45;
            }

            //System.out.println("deltaE = " + deltaE);
            // coefficient is -(1-5)/5! = -1/30
            value += -deltaE/30.0;
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
    
    public void setDeltaCut(double newDeltaDCut) {
        deltaCut2 = newDeltaDCut*newDeltaDCut;
    }
    
    public double getDeltaCut() {
        return Math.sqrt(deltaCut2);
    }
    
    private static final long serialVersionUID = 1L;
    private final ClusterBonds[] clusters;
    private final double[] clusterWeights;
    private final MayerFunction[] f;
    private final double[][][] fValues;
    private long cPairID = -1, lastCPairID = -1;
    private double value, lastValue;
    private double beta;
    protected final MoleculeArrayList scfAtoms;
    protected double deltaCut2 = Double.POSITIVE_INFINITY;
    protected final double[][] uijPol;
    public double pushR2 = 0;
    protected boolean doCaching;
}
