/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.chem.elements.Argon;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.IPotentialMolecular;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.virial.cluster.VirialDiagrams;

/**
 * Potential class which computes non-additive energies from potentials
 * that return the full energies for various-sized groups of molecules.  It is
 * assumed that a single potential can handle all sets of molecules of a given
 * size.
 * 
 * @author Andrew Schultz
 */
public class PotentialNonAdditive implements IPotentialMolecular {

    protected int nBody;
    protected final IPotentialMolecular[] potentials;
    protected final IPotentialMolecular[] p;
    protected final double[][] una;
    protected final MoleculeArrayList moleculeList;
    
    public PotentialNonAdditive(IPotentialMolecular[] potentials) {
        this.potentials = potentials;
        int nb = 0;
        for (int i=0; i<potentials.length; i++) {
            if (potentials[i].nBody() > nb) nb = potentials[i].nBody();
        }
        nBody = nb;
        p = new IPotentialMolecular[nBody+1];
        for (int i=0; i<potentials.length; i++) {
            nb = potentials[i].nBody();
            if (p[nb] != null) {
                throw new RuntimeException("we can handle zero or one "+nb+"-body potentials, but not more");
            }
            p[nb] = potentials[i];
        }
        una = new double[nBody+1][0];
        for (int i=2; i<una.length; i++) {
            int[] ids = new int[i];
            for (int j=0; j<i; j++) {
                ids[j] = nBody-i+j;
            }
            una[i] = new double[VirialDiagrams.getGroupID(ids, nBody)+1];
        }
        moleculeList = new MoleculeArrayList(nBody);
    }

    public PotentialNonAdditive(IPotentialMolecular[] potentials, int[] nbody) {
        this.potentials = potentials;
        int nb = 0;
        for (int i=0; i<potentials.length; i++) {
            if (nbody[i] > nb) nb = nbody[i];
        }
        nBody = nb;
        p = new IPotentialMolecular[nBody+1];
        for (int i=0; i<potentials.length; i++) {
            nb = nbody[i];
            if (p[nb] != null) {
                throw new RuntimeException("we can handle zero or one "+nb+"-body potentials, but not more");
            }
            p[nb] = potentials[i];
        }
        una = new double[nBody+1][0];
        for (int i=2; i<una.length; i++) {
            int[] ids = new int[i];
            for (int j=0; j<i; j++) {
                ids[j] = nBody-i+j;
            }
            una[i] = new double[VirialDiagrams.getGroupID(ids, nBody)+1];
        }
        moleculeList = new MoleculeArrayList(nBody);
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        for (int i=0; i<potentials.length; i++) {
            potentials[i].setBox(box);
        }
    }

    public int nBody() {
        return nBody;
    }

    public double energy(IMoleculeList molecules) {
        
        for (int i=1; i<=nBody; i++) {
            if (p[i] != null) {
                energy(i, molecules);
            }
        }
        return una[nBody][0];
    }
    
    public void energy(int size, IMoleculeList molecules) {
        int nPoints = molecules.getMoleculeCount();
        int groupID = 0;
        boolean debugme = false;
        for(int i=0; i<nPoints; i++) {
            moleculeList.clear();
            moleculeList.add(molecules.getMolecule(i));
            if (size==1) {
                una[1][i] = p[1].energy(moleculeList);
            }
            for(int j=i+1; j<nPoints; j++) {
                moleculeList.add(molecules.getMolecule(j));
                if (size==2) {
                    una[2][VirialDiagrams.pairId(i,j,nBody)] = p[2].energy(moleculeList);
                }
                for (int k=j+1; k<nPoints; k++) {
                    moleculeList.add(molecules.getMolecule(k));
                    double u23 = una[2][VirialDiagrams.pairId(i,j,nBody)] + una[2][VirialDiagrams.pairId(i,k,nBody)] + una[2][VirialDiagrams.pairId(j,k,nBody)];
                    double u33 = 0;
                    if (size==3) {
                        if (u23 < Double.POSITIVE_INFINITY) {
                            una[3][groupID] = p[3].energy(moleculeList) - u23;
                        }
                        else {
                            una[3][groupID] = 0;
                        }
                        u33 = una[3][groupID];
                        if (debugme) {
                            System.out.println("3 "+groupID+" "+ClusterSumMultibody.m2s(moleculeList)+" "+una[3][groupID]);
                        }
                        groupID++;
                        moleculeList.remove(2);
                        continue;
                    }
                    for (int l=k+1; l<nPoints; l++) {
                        moleculeList.add(molecules.getMolecule(l));
                        double u24 = u23 + una[2][VirialDiagrams.pairId(i,l,nBody)] + una[2][VirialDiagrams.pairId(j,l,nBody)] + una[2][VirialDiagrams.pairId(k,l,nBody)];
                        double u34 = u33 + una[3][VirialDiagrams.tripletId(i,j,l,nBody)] + una[3][VirialDiagrams.tripletId(i,k,l,nBody)] + una[3][VirialDiagrams.tripletId(j,k,l,nBody)];
                        double u44 = 0;
                        if (size==4) {
                            if (u24 < Double.POSITIVE_INFINITY) {
                                una[4][groupID] = p[4].energy(moleculeList) - u34 - u24;
                            }
                            else {
                                una[4][groupID] = 0;
                            }
                            u44 = una[4][groupID];
                            if (debugme) {
                                System.out.println("4 "+groupID+" "+ClusterSumMultibody.m2s(moleculeList)+" "+una[4][groupID]);
                            }
                            groupID++;
                            moleculeList.remove(3);
                            continue;
                        }
                        
                        for (int m=l+1; m<nPoints; m++) {
                            moleculeList.add(molecules.getMolecule(m));
                            double u25 = u24 + una[2][VirialDiagrams.pairId(i,m,nBody)] + una[2][VirialDiagrams.pairId(j,m,nBody)] + una[2][VirialDiagrams.pairId(k,m,nBody)]
                                    + una[2][VirialDiagrams.pairId(l,m,nBody)];
                            double u35 = u34 + una[3][VirialDiagrams.tripletId(i,j,m,nBody)] + una[3][VirialDiagrams.tripletId(i,k,m,nBody)] + una[3][VirialDiagrams.tripletId(i,l,m,nBody)]
                                             + una[3][VirialDiagrams.tripletId(j,k,m,nBody)] + una[3][VirialDiagrams.tripletId(j,l,m,nBody)] + una[3][VirialDiagrams.tripletId(k,l,m,nBody)];
                            double u45 = u44 + una[4][VirialDiagrams.quadId(i,j,k,m,nBody)] + una[4][VirialDiagrams.quadId(i,j,l,m,nBody)] + una[4][VirialDiagrams.quadId(i,k,l,m,nBody)]
                                    + una[4][VirialDiagrams.quadId(j,k,l,m,nBody)];
                            double u55 = 0;
                            if (size==5) {
                                if (u25 < Double.POSITIVE_INFINITY) {
                                    una[5][groupID] = p[5].energy(moleculeList) - u45 - u35 - u25;
                                }
                                else {
                                    una[5][groupID] = 0;
                                }
                                u55 = una[5][groupID];
                                if (debugme) {
                                    System.out.println("5 "+groupID+" "+ClusterSumMultibody.m2s(moleculeList)+" "+una[5][groupID]);
                                }
                                groupID++;
                                moleculeList.remove(4);
                                continue;
                            }

                            for (int mm=m+1; mm<nPoints; mm++) {
                                moleculeList.add(molecules.getMolecule(mm));
                                double u26 = u25 + una[2][VirialDiagrams.pairId(i,mm,nBody)] + una[2][VirialDiagrams.pairId(j,mm,nBody)] + una[2][VirialDiagrams.pairId(k,mm,nBody)]
                                        + una[2][VirialDiagrams.pairId(l,mm,nBody)] + una[2][VirialDiagrams.pairId(m,mm,nBody)];
                                double u36 = u35 + una[3][VirialDiagrams.tripletId(i,j,mm,nBody)] + una[3][VirialDiagrams.tripletId(i,k,mm,nBody)] + una[3][VirialDiagrams.tripletId(i,l,mm,nBody)]
                                                 + una[3][VirialDiagrams.tripletId(j,k,mm,nBody)] + una[3][VirialDiagrams.tripletId(j,l,mm,nBody)] + una[3][VirialDiagrams.tripletId(k,l,mm,nBody)]
                                                 + una[3][VirialDiagrams.tripletId(i,m,mm,nBody)] + una[3][VirialDiagrams.tripletId(j,m,mm,nBody)] + una[3][VirialDiagrams.tripletId(k,m,mm,nBody)]
                                                 + una[3][VirialDiagrams.tripletId(l,m,mm,nBody)];
                                double u46 = u45 + una[4][VirialDiagrams.quadId(i,j,k,mm,nBody)] + una[4][VirialDiagrams.quadId(i,j,l,mm,nBody)] + una[4][VirialDiagrams.quadId(i,k,l,mm,nBody)]
                                                 + una[4][VirialDiagrams.quadId(j,k,l,mm,nBody)] + una[4][VirialDiagrams.quadId(i,j,m,mm,nBody)] + una[4][VirialDiagrams.quadId(i,k,m,mm,nBody)]
                                                 + una[4][VirialDiagrams.quadId(i,l,m,mm,nBody)] + una[4][VirialDiagrams.quadId(j,k,l,mm,nBody)] + una[4][VirialDiagrams.quadId(j,k,m,mm,nBody)]
                                                 + una[4][VirialDiagrams.quadId(k,l,m,mm,nBody)];
                                double u56 = u55 + una[5][VirialDiagrams.quintId(i,j,k,l,mm,nBody)] + una[5][VirialDiagrams.quintId(i,j,l,m,mm,nBody)] + una[5][VirialDiagrams.quintId(i,k,l,m,mm,nBody)]
                                                 + una[5][VirialDiagrams.quintId(j,k,l,m,mm,nBody)];
                                if (size==6) {
                                    if (u26 < Double.POSITIVE_INFINITY) {
                                        una[6][groupID] = p[6].energy(moleculeList) - u56 - u46 - u36 - u26;
                                    }
                                    else {
                                        una[6][groupID] = 0;
                                    }
                                    if (debugme) {
                                        System.out.println("6 "+groupID+" "+ClusterSumMultibody.m2s(moleculeList)+" "+una[6][groupID]);
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

    public PotentialNonAdditiveNB makeNB(int nBodyPretend) {
        return new PotentialNonAdditiveNB(nBodyPretend);
    }

    public class PotentialNonAdditiveNB implements IPotentialMolecular {
        protected int nBodyPretend;
        public PotentialNonAdditiveNB(int nBody) {
            nBodyPretend = nBody;
        }
        public int nBody() {
            return nBodyPretend;
        }

        public double getRange() {
            return PotentialNonAdditive.this.getRange();
        }

        public void setBox(Box box) {
            PotentialNonAdditive.this.setBox(box);
        }

        public double energy(IMoleculeList molecules) {
            return PotentialNonAdditive.this.energy(molecules);
        }
    }

    public static void main(String[] args) {

        Space space = Space3D.getInstance();
        PotentialEmulCached p3 = new PotentialEmulCached(space, "3body_template.in", 3);
        PotentialEmulCached p2 = new PotentialEmulCached(space, "2body_template.in", 3);
        Simulation sim = new Simulation(space);
        BoxCluster box = new BoxCluster(null, space);
        sim.addBox(box);
        SpeciesSpheresMono species = new SpeciesSpheresMono(space, Argon.INSTANCE);
        sim.addSpecies(species);
        box.setNMolecules(species, 3);
        p2.setBox(box);
        p3.setBox(box);
        box.getLeafList().getAtom(0).getPosition().setX(0, -4);
        box.getLeafList().getAtom(1).getPosition().setX(0, -4);
        box.getLeafList().getAtom(0).getPosition().setX(1, -3);
        box.getLeafList().getAtom(1).getPosition().setX(1, 3);
        MoleculeArrayList molecules02 = new MoleculeArrayList(2);
        molecules02.add(box.getMoleculeList().getMolecule(0));
        molecules02.add(box.getMoleculeList().getMolecule(1));
        double u01 = p2.energy(molecules02);
        molecules02.clear();
        molecules02.add(box.getMoleculeList().getMolecule(0));
        molecules02.add(box.getMoleculeList().getMolecule(2));
        MoleculeArrayList molecules12 = new MoleculeArrayList(2);
        molecules12.add(box.getMoleculeList().getMolecule(1));
        molecules12.add(box.getMoleculeList().getMolecule(2));
        
        PotentialNonAdditive p3NA = new PotentialNonAdditive(new IPotentialMolecular[]{p2, p3});
        p3NA.setBox(box);
        
        for (int i=0; i<100; i++) {
            double r = -3 + i*0.1;
            box.getLeafList().getAtom(2).getPosition().setX(0, r);
            box.trialNotify();
            box.acceptNotify();
            p3NA.setBox(box);
            double uEmul = p3.energy(box.getMoleculeList(species));
            double u02 = p2.energy(molecules02);
            double u12 = p2.energy(molecules12);
            System.out.printf("%5.2f ", r);
            System.out.println(uEmul+" "+(uEmul - (u01+u02+u12))+" "+p3NA.energy(box.getMoleculeList(species)));
        }
    }
}
