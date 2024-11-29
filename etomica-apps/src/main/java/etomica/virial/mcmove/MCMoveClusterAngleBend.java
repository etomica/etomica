/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeSource;
import etomica.molecule.MoleculeSourceRandomMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * An MC Move for cluster simulations that bends the bond angle for 3-atom
 * molecule.  The COM is not moved and both end atoms are moved by the same
 * amount such that the move does not alter the orientation of the molecule.
 *
 * @author Andrew Schultz
 */
public class MCMoveClusterAngleBend extends MCMoveBoxStep {

    protected final PotentialCompute pc;
    protected final Vector work1, work2, work3;
    protected double[] dTheta;
    protected double wOld, wNew;
    protected ISpecies species;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected final IRandom random;
    protected MoleculeSource moleculeSource;

    public MCMoveClusterAngleBend(IRandom random, PotentialCompute pc, Space _space) {
        this(pc, random, 1.0, _space);
    }

    public MCMoveClusterAngleBend(PotentialCompute pc,
                                  IRandom random, double stepSize, Space _space) {
        super();
        this.pc = pc;
        this.random = random;
        moleculeSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule) moleculeSource).setRandomNumberGenerator(random);
        setStepSizeMax(Math.PI / 2);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;

        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        dTheta = new double[p.getMoleculeList().size()];
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    public boolean doTrial() {
        uOld = pc.computeAll(false);
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<moleculeList.size(); i++) {
            if (species != null && moleculeList.get(i).getType() != species) {
                continue;
            }
            //hack for angle bend of root molecules only
            if(i ==0) {
                IMolecule molecule = moleculeList.get(i);
                IAtomList childList = molecule.getChildList();
                int numChildren = childList.size();
                if (numChildren != 3) continue;
                double dt = stepSize * (random.nextDouble() - 0.5);
                dTheta[i] = dt;
                transform(molecule, dt);
            }
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        uNew = pc.computeAll(false);

        return true;
    }
    
    protected void transform(IMolecule molecule, double dt) {
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();
        if (numChildren != 3) return;
        
        Vector pos0 = childList.get(0).getPosition();
        Vector pos1 = childList.get(1).getPosition();
        Vector pos2 = childList.get(2).getPosition();
        
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

        double cdt = Math.cos(dt);
        double sdt = Math.sin(dt);
        double m0 = childList.get(0).getType().getMass();
        double m1 = childList.get(1).getType().getMass();
        double m2 = childList.get(2).getType().getMass();

//        work3.Ea1Tv1(m0, pos0);
//        pos0.E(pos1);
//        pos0.PEa1Tv1(bondLength01*cdt, work1);
//        pos0.PEa1Tv1(bondLength01*sdt, work2);
//        work3.PEa1Tv1(-m0, pos0);

        // normalize bond lengths -- we'll scale our vectors back up later
        work2.Ev1Mv2(pos2, pos1);
        work2.TE(1.0/bondLength12);
        dot = work1.dot(work2);
        work1.PEa1Tv1(-dot, work2);
        work1.TE(1.0/Math.sqrt(work1.squared()));

        work3.PEa1Tv1(m2, pos2);
        pos2.E(pos1);
        pos2.PEa1Tv1(bondLength12*cdt, work2);
        pos2.PEa1Tv1(bondLength12*sdt, work1);
        work3.PEa1Tv1(-m2, pos2);
        
        // translate COM back to its original position
//        work3.TE(1/(m0+m1+m2));
//        pos0.PE(work3);
//        pos1.PE(work3);
//        pos2.PE(work3);
    }
    
    public void acceptNotify() {
//        System.out.println("accepted angle");
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
//        System.out.println("rejected wiggle");
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<box.getMoleculeList().size(); i++) {
            if (dTheta[i] == 0) continue;
            transform(moleculeList.get(i), -dTheta[i]);
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }
	
    public double energyChange() {return uNew - uOld;}

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

}
