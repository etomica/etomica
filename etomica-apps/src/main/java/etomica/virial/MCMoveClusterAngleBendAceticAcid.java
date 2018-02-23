/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * An MC Move for cluster simulations that bends the bond angle for 3-atom
 * molecule.  The COM is not moved and both end atoms are moved by the same
 * amount such that the move does not alter the orientation of the molecule.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterAngleBendAceticAcid extends MCMoveBoxStep {

    public MCMoveClusterAngleBendAceticAcid(Simulation sim, PotentialMaster potentialMaster, Space _space) {
    	this(potentialMaster,sim.getRandom(), 1.0, _space);
    }
    
    public MCMoveClusterAngleBendAceticAcid(PotentialMaster potentialMaster,
                                            IRandom random, double stepSize, Space _space) {
        super(potentialMaster);
        this.space = _space;
        this.random = random;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)moleculeSource).setRandomNumberGenerator(random);
        setStepSizeMax(Math.PI/2);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;

        energyMeter = new MeterPotentialEnergy(potential);
        energyMeter.setIncludeLrc(false);
        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        dTheta = new double[p.getMoleculeList().getMoleculeCount()];
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getMoleculeCount(); i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) {
                continue;
            }
            IMolecule molecule = moleculeList.getMolecule(i);
            IAtomList childList = molecule.getChildList();
            int numChildren = childList.size();
            if (numChildren != 5) continue;
            double dt = stepSize * (random.nextDouble() - 0.5);
            dTheta[i] = dt;
            transform(molecule, dt);
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    protected void transform(IMolecule molecule, double dt) {
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();
        if (numChildren != 3) return;
        
        Vector pos0 = childList.get(3).getPosition();//SBO
        Vector pos1 = childList.get(1).getPosition();//C
        Vector pos2 = childList.get(2).getPosition();//DBO
        Vector pos3 = childList.get(4).getPosition();//H
        Vector pos4 = childList.get(0).getPosition();//CH3
        
        work1.Ev1Mv2(pos0, pos1);
        double bondLength01 = Math.sqrt(work1.squared());
        // normalize bond lengths -- we'll scale our vectors back up later
        work1.TE(1.0/bondLength01);
        work2.Ev1Mv2(pos2, pos1);
        double bondLength12 = Math.sqrt(work2.squared());
        work2.TE(1.0/bondLength12);
        double dot = work1.dot(work2);
        work2.PEa1Tv1(-dot, work1);
        work2.TE(1.0/Math.sqrt(work2.squared()));
        work3.Ev1Mv2(pos3, pos1);
        double bondLength13 = Math.sqrt(work3.squared());

        double cdt = Math.cos(dt);
        double sdt = Math.sin(dt);
        
        work3.E(pos0);
        pos0.E(pos1);
        pos0.PEa1Tv1(bondLength01*cdt, work1);
        pos0.PEa1Tv1(bondLength01*sdt, work2);
        work3.ME(pos0);

        // normalize bond lengths -- we'll scaled our vectors back up later
        work2.Ev1Mv2(pos2, pos1);
        work2.TE(1.0/bondLength12);
        dot = work1.dot(work2);
        work1.PEa1Tv1(-dot, work2);
        work1.TE(1.0/Math.sqrt(work1.squared()));

        work3.PE(pos2);
        pos2.E(pos1);
        pos2.PEa1Tv1(bondLength12*cdt, work2);
        pos2.PEa1Tv1(bondLength12*sdt, work1);
        work3.ME(pos2);
        
        work2.Ev1Mv2(pos3, pos1);
        work2.TE(1.0/bondLength13);
        dot = work1.dot(work2);
        work1.PEa1Tv1(-dot, work2);
        work1.TE(1.0/Math.sqrt(work1.squared()));

        work3.PE(pos3);
        pos3.E(pos1);
        pos3.PEa1Tv1(bondLength13*cdt, work2);
        pos3.PEa1Tv1(bondLength13*sdt, work1);
        work3.ME(pos3);
        
        // translate COM back to its original position
        work3.TE(1.0/5.0);//5 is number of atoms(acetic acid)
        pos0.PE(work3);
        pos1.PE(work3);
        pos2.PE(work3);
        pos3.PE(work3);
        pos4.PE(work3);
    }
    
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<box.getMoleculeList().getMoleculeCount(); i++) {
            if (dTheta[i] == 0) continue;
            transform(moleculeList.getMolecule(i), -dTheta[i]);
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    public double energyChange() {return uNew - uOld;}

    public AtomIterator affectedAtoms() {
        affectedAtomIterator.setList(box.getLeafList());
        return affectedAtomIterator;
    }
    
    /**
     * @return Returns the atomSource.
     */
    public MoleculeSource getAtomSource() {
        return moleculeSource;
    }
    /**
     * @param source The atomSource to set.
     */
    public void setMoleculeSource(MoleculeSource source) {
        moleculeSource = source;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomIteratorArrayListSimple affectedAtomIterator = new AtomIteratorArrayListSimple();
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector work1, work2, work3;
    protected double[] dTheta;
    protected double wOld, wNew;
    protected final Space space;
    protected ISpecies species;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected final IRandom random;
    protected MoleculeSource moleculeSource;
}
