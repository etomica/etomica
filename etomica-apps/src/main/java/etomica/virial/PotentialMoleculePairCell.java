/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.NeighborIteratorCell;
import etomica.nbr.cell.NeighborIteratorCellFaster;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMoleculePair;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;
import etomica.nbr.cell.NeighborCellManagerMulti;

import java.util.ArrayList;

public class PotentialMoleculePairCell implements IPotentialMolecular {

    protected final IPotential2[][] atomPotentials;
    protected BoxCluster box;
    protected double [][]energies,lastEnergies;
    protected long cPairID = -1, lastCPairID = -1;
    protected NeighborIterator iterator;
    public PotentialCompute pcheck;

    public PotentialMoleculePairCell( SpeciesManager sm) {
        this( makeAtomPotentials(sm));
    }
    public boolean doDebug = false;

    public PotentialMoleculePairCell( IPotential2[][] atomPotentials) {

        this.atomPotentials = atomPotentials;


    }

    public void setBox(BoxCluster box) {
        this.box = box;
        int n=box.getMoleculeList().size();
        energies=new double[n][n];
        lastEnergies=new double[n][n];
    }

    public void setNeighborManager(NeighborManager neighborManager){
        iterator = neighborManager.makeNeighborIterator();
    }

    private static IPotential2[][] makeAtomPotentials(SpeciesManager sm) {
        // we could try to store the potentials more compactly, but it doesn't really matter
        ISpecies species = sm.getSpecies(sm.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        return new IPotential2[lastTypeIndex + 1][lastTypeIndex + 1];
    }

    public void setAtomPotential(AtomType atomType1, AtomType atomType2, IPotential2 p2) {
        atomPotentials[atomType1.getIndex()][atomType2.getIndex()] = p2;
        atomPotentials[atomType2.getIndex()][atomType1.getIndex()] = p2;
    }

    public IPotential2[][] getAtomPotentials() {
        return atomPotentials;
    }

    @Override
    public double energy(IMoleculeList molecules) {
        return energy(molecules.get(0), molecules.get(1));
    }
public  long callCount,pairCount;
    public IntegratorBox integrator;
    public double energy(IMolecule molecule1, IMolecule molecule2) {
       // System.out.println("computing energy " +molecule1+ " " +molecule2);
        //new


        CoordinatePairSet cPairs = box.getCPairSet();
        long thisCPairID = cPairs.getID();
//            System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) {
//                System.out.println("clusterSum "+cPairID+" returning recent "+value);
            return energies[molecule1.getIndex()][molecule2.getIndex()];
        }
        if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            for (int i =0;i < energies.length;i ++){
                for (int j =0;j< energies.length; j ++){
                    energies[i][j]=lastEnergies[i][j];
                }
            }
//                System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
            return energies[molecule1.getIndex()][molecule2.getIndex()];
        }
callCount++;
        // a new cluster
        lastCPairID = cPairID;
        for (int i =0;i < energies.length;i ++){
            for (int j =0;j< energies.length; j ++){
                lastEnergies[i][j]=energies[i][j];
                energies[i][j]=0;
            }
        }
        cPairID = thisCPairID;
        IMoleculeList molecules = box.getMoleculeList();
        ArrayList[] pairs = null;
       if(doDebug){
           pairs=new ArrayList[box.getLeafList().size()];
           for(int i=0;i<pairs.length;i++) {
               pairs[i] = new ArrayList();
           }
           NeighborIteratorCellFaster.doDebug=false;
       }

        for (IMolecule m0 :molecules){
            for (IAtom a0:m0.getChildList()){

                IPotential2[] p0 = atomPotentials[a0.getType().getIndex()];
                ArrayList[] finalPairs = pairs;
                iterator.iterUpNeighbors(a0.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom a1, Vector rij, int n) {

                        IPotential2 p2 = p0[a1.getType().getIndex()];
                        if (p2 == null) return;
                        IMolecule m1 = a1.getParentGroup();
                        //if(m1==m0) return;
                        if(doDebug && m1==m0){
                            finalPairs[a0.getLeafIndex()].add(a1);
                        }

                        double u = p2.u(rij.squared());
if (u!=0) pairCount++;
                        //check
                        //if (doDebug && u != 0 && m0 != m1)

                            //System.out.println("WITH_CELL: " + a0.getLeafIndex() + " " + a1.getLeafIndex() + " " + rij.squared() + " " + u);

                        //energies[m0.getIndex()] [m1.getIndex()]+= u;
                       // if(i0==4 ||i1==4)System.out.println(Math.min(i0,i1)+" "+ Math.max(i0,i1)+" "+rij.squared()+" "+ u);
                        energies[m0.getIndex()] [m1.getIndex()]+= u;
                        if(m0!=m1) energies[m1.getIndex()] [m0.getIndex()]+= u;

                    }
                });

            }

        }
        if(doDebug) {
            ArrayList[] finalPairs1 = pairs;
            double Ucheck = pcheck.computeAll(false, new PotentialCallback() {
                @Override
                public void pairCompute(int i, int j, Vector dr, double[] u012) {
                    boolean found = false;
                    for (int k = 0; k < finalPairs1[i].size(); k++) {
                        if (((IAtom) finalPairs1[i].get(k)).getLeafIndex() == j) {
                            found = true;
                        }
                    }
                    for (int k = 0; k < finalPairs1[j].size(); k++) {
                        if (((IAtom) finalPairs1[j].get(k)).getLeafIndex() == i) {
                            found = true;
                        }

                    }
                    if (!found) {
                        throw new RuntimeException("couldn't find " + i + " " + j);
                    }
                }
            });
            System.out.println(Ucheck + " " + (energies[0][0] + energies[1][1]));
            //System.exit(0);
        }





        return energies[molecule1.getIndex()][molecule2.getIndex()];

    }
    public PotentialComputeIntra makePotentialComputeIntra(){
        return new PotentialComputeIntra();
    }
    public class PotentialComputeIntra implements PotentialCompute{
        protected double lastEnergy;
        @Override
        public void init() {

        }

        @Override
        public Vector[] getForces() {
            return new Vector[0];
        }

        @Override
        public double getLastVirial() {
            return 0;
        }

        @Override
        public double getLastEnergy() {
            return lastEnergy;
        }

        @Override
        public void updateAtom(IAtom atom) {

        }

        @Override
        public double computeAll(boolean doForces, PotentialCallback pc) {
            double uintra = 0;
            IMoleculeList molecules = box.getMoleculeList();
            for (IMolecule m:molecules){
                uintra+=energy(m,m);
            }
            lastEnergy=uintra;
            return uintra;
        }

        @Override
        public double computeOneOld(IAtom iAtom) {
            return 0;
        }

        @Override
        public double computeOne(IAtom iAtom) {
            return 0;
        }

        @Override
        public double computeManyAtomsOld(IAtom... atoms) {
            return 0;
        }

        @Override
        public double computeManyAtoms(IAtom... atoms) {
            return 0;
        }

        @Override
        public void processAtomU(double fac) {

        }

        @Override
        public IntegratorListener makeIntegratorListener() {
            return new IntegratorListener() {};
        }
    }


}
