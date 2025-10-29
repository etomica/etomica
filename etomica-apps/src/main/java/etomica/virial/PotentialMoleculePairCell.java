/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.virial;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.NeighborCellManager;
import etomica.potential.IPotential2;
import etomica.potential.IPotentialMolecular;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.species.SpeciesManager;

public class PotentialMoleculePairCell implements IPotentialMolecular {

    protected final IPotential2[][] atomPotentials;
    protected BoxCluster box;
    protected double [][]energies,lastEnergies;
    protected long cPairID = -1, lastCPairID = -1;
    protected NeighborIterator iterator;

    public PotentialMoleculePairCell( SpeciesManager sm) {
        this( makeAtomPotentials(sm));
    }

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

    public double energy(IMolecule molecule1, IMolecule molecule2) {
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
        for (IMolecule m0 :molecules){
            for (IAtom a0:m0.getChildList()){
                IPotential2[] p0 = atomPotentials[a0.getType().getIndex()];
                iterator.iterUpNeighbors(a0.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom a1, Vector rij, int n) {
                        IPotential2 p2 = p0[a1.getType().getIndex()];
                        if (p2 == null) return;
                        IMolecule m1 = a1.getParentGroup();
                        if(m1==m0) return;
                        double u = p2.u(rij.squared());
                        energies[m0.getIndex()] [m1.getIndex()]+= u;
                        energies[m1.getIndex()] [m0.getIndex()]+= u;

                    }
                });
            }

        }

        return energies[molecule1.getIndex()][molecule2.getIndex()];
    }

}
