/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IMoleculeList;
import etomica.atom.MoleculeArrayList;
import etomica.util.Arrays;
import etomica.virial.cluster.VirialDiagrams;


/**
 * Cluster class capable of computing multibody (up to 6th order) contributions
 * to diagrams.
 * 
 * Currently, only one type of non-additive Mayer function may be used, although
 * some framework is in place to handle more than one.
 * 
 * @author Andrew Schultz
 */
public class ClusterSumMultibody extends ClusterSum {

    public ClusterSumMultibody(ClusterBonds[] subClusters, double[] subClusterWeights, MayerFunction[] fArray, MayerFunctionNonAdditive[] fNonAdditive) {
        super(subClusters, subClusterWeights, fArray);
        
        if (!clusters[0].isUsePermutations()) {
            int pointCount = subClusters[0].pointCount();
            // determine which fbonds are actually needed by the diagrams
            for (int i=0; i<pointCount-1; i++) {
                for (int j=i+1; j<pointCount; j++) {
                    int[] ff = fullBondIndexArray[i][j];
                    boolean has0 = false;
                    for (int k=0; k<ff.length; k++) {
                        if (ff[k] == 0) {
                            // we'll already calculate MayerFunction 0 for the i-j pair
                            has0 = true;
                            break;
                        }
                    }
                    if (!has0) {
                        // force f-bond to be calculated (we need it to get e)
                        fullBondIndexArray[i][j] = Arrays.resizeArray(ff, ff.length+1);
                        fullBondIndexArray[i][j][ff.length] = 0;
                    }
                }
            }
        }

        
        this.fNonAdditive = fNonAdditive;
        moleculeList = new MoleculeArrayList(3);
        int nPoints = subClusters[0].nPoints;
        fNonAdditiveValues = new double[nPoints+1][];
        fNonAdditiveValues[3] = new double[VirialDiagrams.tripletId(nPoints-3, nPoints-2, nPoints-1, nPoints)+1];
        if (nPoints>3) {
            fNonAdditiveValues[4] = new double[VirialDiagrams.quadId(nPoints-4, nPoints-3, nPoints-2, nPoints-1, nPoints)+1];
            if (nPoints>4) {
                fNonAdditiveValues[5] = new double[VirialDiagrams.quintId(nPoints-5, nPoints-4, nPoints-3, nPoints-2, nPoints-1, nPoints)+1];
                if (nPoints>5) {
                    fNonAdditiveValues[6] = new double[VirialDiagrams.sixId(nPoints-6, nPoints-5, nPoints-4, nPoints-3, nPoints-2, nPoints-1, nPoints)+1];
                }
            }
        }
        // brute force search, but will help if someone has diagrams where
        // some of the atom groups have no multibody interaction, or a
        // different multibody interaction.
        fNonAdditiveNeeded = new int[nPoints+1][0];
        numNonAdditiveNeeded = new int[nPoints+1];
        for (int i=3; i<nPoints+1; i++) {
            fNonAdditiveNeeded[i] = new int[fNonAdditiveValues[i].length];
            for (int k=0; k<fNonAdditiveValues[i].length; k++) {
                doSearch(i, k);
            }
        }
        r2 = new double[nPoints*(nPoints-1)/2];
    }

    protected void doSearch(int size, int id) {
        boolean found = false;
        for (int j=0; !found && j<clusters.length; j++) {
            if (clusters[j] instanceof ClusterBondsNonAdditive) {
                int[][] multiBonds = ((ClusterBondsNonAdditive)clusters[j]).getMultiBonds();
                if (multiBonds.length <= size) continue;
                for (int l=0; l<multiBonds[size].length; l++) {
                    if (multiBonds[size][l] == id) {
                        found = true;
                        break;
                    }
                }
            }
        }
        if (found) {
            fNonAdditiveNeeded[size][numNonAdditiveNeeded[size]] = id;
            numNonAdditiveNeeded[size]++;
        }
    }

    public ClusterAbstract makeCopy() {
        ClusterSumMultibody copy = new ClusterSumMultibody(clusters,clusterWeights,f,fNonAdditive);
        copy.setTemperature(1/beta);
        return copy;
    }
    
    public double value(BoxCluster box) {
        super.value(box);
        if (pushR2 > 0) {
            CoordinatePairSet cPairs = box.getCPairSet();
            double r2max = 0;
            for (int i=0; i<clusters[0].nPoints-1; i++) {
                for (int j=i+1; j<clusters[0].nPoints; j++) {
                    double r2ij = cPairs.getr2(i,j);
                    if (r2ij > r2max) {
                        r2max = r2ij;
                    }
                }
            }
            if (r2max < pushR2) {
                value *= Math.pow((r2max/pushR2),4);
            }
            
        }
        return value;
    }

    protected void calcValue() {
        value = 0.0;
        for(int i=0; i<clusters.length; i++) {
            double v = 0;
            if (clusters[i] instanceof ClusterBondsNonAdditive) {
                v = ((ClusterBondsNonAdditive)clusters[i]).value(fValues, fNonAdditiveValues);
            }
            else {
                v = clusters[i].value(fValues);
            }
            value += clusterWeights[i] * v;
        }
    }

    protected void updateF(BoxCluster box) {
        super.updateF(box);
        int nPoints = pointCount();
        CoordinatePairSet cPairs = box.getCPairSet();

        for (int k=0; k<fNonAdditive.length; k++) {
            fNonAdditive[k].setBox(box);
        }

        IMoleculeList molecules = box.getMoleculeList();
        // we do all 3-body, then all 4-body, etc
        // MayerFunctionThreeBody assumes that all 3-body calls will happen before the rest.
        for (int i=3; i<=nPoints; i++) {
            doUpdateFMulti(i, molecules, cPairs);
        }
    }

    public static String m2s(IMoleculeList list) {
        String s = ""+list.getMolecule(0).getIndex();
        for (int i=1; i<list.getMoleculeCount(); i++) {
            s += " "+list.getMolecule(i).getIndex();
        }
        return s;
    }

    /*
     * This looks nasty, but is probably still less expensive than computing
     * the multibody potentials.  A more generic while loop implementation
     * could be done, but wouldn't really help (getting up to 6 molecules with
     * for loops isn't crazy).
     */
    protected void doUpdateFMulti(int size, IMoleculeList molecules, CoordinatePairSet cPairs) {
        int nPoints = pointCount();
        int nextNeeded = fNonAdditiveNeeded[size].length > 0 ? fNonAdditiveNeeded[size][0] : -1;
        int groupID = 0, iGroup = 0;
        boolean debugme = false; //nfar >= 9 && nclose >= 3;
        if (debugme && size == 3) {
            for(int i=0; i<nPoints-1; i++) {
                for(int j=i+1; j<nPoints; j++) {
                    System.out.println(i+" "+j+" "+Math.sqrt(cPairs.getr2(i,j))+" "+fValues[i][j][0]);
                }
            }
        }
        for(int i=0; i<nPoints-2; i++) {
            if (nextNeeded == -1) return;
            moleculeList.clear();
            moleculeList.add(molecules.getMolecule(i));
            for(int j=i+1; j<nPoints-1; j++) {
                r2[0] = cPairs.getr2(i,j);
                moleculeList.add(molecules.getMolecule(j));
                for (int k=j+1; k<nPoints; k++) {
                    r2[1] = cPairs.getr2(i,k);
                    r2[2] = cPairs.getr2(j,k);
                    moleculeList.add(molecules.getMolecule(k));
                    double eProduct3 = (fValues[i][j][0]+1)*(fValues[i][k][0]+1)*(fValues[j][k][0]+1);
                    if (size==3) {
                        if (groupID == nextNeeded) {
                            if (eProduct3 != 0) {
                                fNonAdditiveValues[3][groupID] = eProduct3*fNonAdditive[0].f(moleculeList, r2, beta);
                            }
                            else {
                                fNonAdditiveValues[3][groupID] = 0;
                            }
                            if (debugme) {
                                System.out.println("3 "+groupID+" "+m2s(moleculeList)+" "+fNonAdditiveValues[3][groupID]);
                            }
                            
                            iGroup++;
                            nextNeeded = iGroup == fNonAdditiveNeeded[3].length ? -1 : fNonAdditiveNeeded[3][iGroup];
                        }
                        groupID++;
                        moleculeList.remove(2);
                        continue;
                    }
                    for (int l=k+1; l<nPoints; l++) {
                        r2[3] = cPairs.getr2(i,l);
                        r2[4] = cPairs.getr2(j,l);
                        r2[5] = cPairs.getr2(k,l);
                        moleculeList.add(molecules.getMolecule(l));
                        double eProduct4 = eProduct3 * (fValues[i][l][0]+1)*(fValues[j][l][0]+1)*(fValues[k][l][0]+1);
                        if (size==4) {
                            if (groupID == nextNeeded) {
                                if (eProduct4!= 0) {
                                    fNonAdditiveValues[4][groupID] = eProduct4*fNonAdditive[0].f(moleculeList, r2, beta);
                                }
                                else {
                                    fNonAdditiveValues[4][groupID] = 0;
                                }
                                if (debugme) {
                                    System.out.println("4 "+groupID+" "+m2s(moleculeList)+" "+fNonAdditiveValues[4][groupID]);
                                }
                                iGroup++;
                                nextNeeded = iGroup == fNonAdditiveNeeded[4].length ? -1 : fNonAdditiveNeeded[4][iGroup];
                            }
                            groupID++;
                            moleculeList.remove(3);
                            continue;
                        }
                        
                        for (int m=l+1; m<nPoints; m++) {
                            r2[6] = cPairs.getr2(i,m);
                            r2[7] = cPairs.getr2(j,m);
                            r2[8] = cPairs.getr2(k,m);
                            r2[9] = cPairs.getr2(l,m);
                            moleculeList.add(molecules.getMolecule(m));
                            double eProduct5 = eProduct4 * (fValues[i][m][0]+1)*(fValues[j][m][0]+1)*(fValues[k][m][0]+1)*(fValues[l][m][0]+1);
                            if (size==5) {
                                if (groupID == nextNeeded) {
                                    if (eProduct5!= 0) {
                                        fNonAdditiveValues[5][groupID] = eProduct5*fNonAdditive[0].f(moleculeList, r2, beta);
                                    }
                                    else {
                                        fNonAdditiveValues[5][groupID] = 0;
                                    }
                                    if (debugme) {
                                        System.out.println("5 "+groupID+" "+m2s(moleculeList)+" "+fNonAdditiveValues[5][groupID]);
                                    }
                                    iGroup++;
                                    nextNeeded = iGroup == fNonAdditiveNeeded[5].length ? -1 : fNonAdditiveNeeded[5][iGroup];
                                }
                                groupID++;
                                moleculeList.remove(4);
                                continue;
                            }

                            for (int mm=m+1; mm<nPoints; mm++) {
                                r2[10] = cPairs.getr2(i,mm);
                                r2[11] = cPairs.getr2(j,mm);
                                r2[12] = cPairs.getr2(k,mm);
                                r2[13] = cPairs.getr2(l,mm);
                                r2[14] = cPairs.getr2(m,mm);
                                moleculeList.add(molecules.getMolecule(mm));
                                double eProduct6 = eProduct5 * (fValues[i][mm][0]+1)*(fValues[j][mm][0]+1)*(fValues[k][mm][0]+1)*
                                                               (fValues[l][mm][0]+1)*(fValues[m][mm][0]+1);
                                if (size==6) {
                                    if (groupID == nextNeeded) {
                                        if (eProduct6!= 0) {
                                            fNonAdditiveValues[6][groupID] = eProduct6*fNonAdditive[0].f(moleculeList, r2, beta);
                                        }
                                        else {
                                            fNonAdditiveValues[6][groupID] = 0;
                                        }
                                        if (debugme) {
                                            System.out.println("6 "+groupID+" "+m2s(moleculeList)+" "+fNonAdditiveValues[6][groupID]);
                                        }
                                        iGroup++;
                                        nextNeeded = iGroup == fNonAdditiveNeeded[6].length ? -1 : fNonAdditiveNeeded[6][iGroup];
                                    }
                                    groupID++;
                                    moleculeList.remove(5);
                                    continue;
                                }
                            }
                            moleculeList.remove(4);
                        }
                        moleculeList.remove(3);
                    }
                    moleculeList.remove(2);
                }
                moleculeList.remove(1);
            }
            moleculeList.clear();
        }
    }

    public double[][] getFNonAdditiveValues() {
        return fNonAdditiveValues;
    }

    private static final long serialVersionUID = 1L;
    protected final MayerFunctionNonAdditive[] fNonAdditive;
    protected double[][] fNonAdditiveValues;
    protected final MoleculeArrayList moleculeList;
    protected final int[][] fNonAdditiveNeeded;
    protected final int[] numNonAdditiveNeeded;
    protected final double[] r2;
    protected final double pushR2 = 0;
}
