/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.IVectorMutable;
import etomica.atom.MoleculeSource;
import etomica.atom.MoleculeSourceRandomMolecule;
import etomica.atom.iterator.AtomIterator;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class MCMoveReptate extends MCMoveBoxStep {
    
    public MCMoveReptate(ISimulation sim, IPotentialMaster potentialMaster, ISpace _space){
        this(potentialMaster, sim.getRandom(), 1.0, 15.0, false, _space);
    }
    
    public MCMoveReptate(IPotentialMaster potentialMaster, IRandom random, 
            double stepSize, double stepSizeMax, boolean fixOverlap, ISpace _space){
        super(potentialMaster);
        this.space = _space;
        this.random = random;
        atomSource = new MoleculeSourceRandomMolecule();
        ((MoleculeSourceRandomMolecule)atomSource).setRandomNumberGenerator(random);
        energyMeter = new MeterPotentialEnergy(potentialMaster);
//        conf = con;
//        if(!(conf instanceof ConformationChainZigZag)){
//            throw new IllegalArgumentException("MCMoveRepate will only work with ConformationZigZag right now.  Sorry about that.");
//        }
        
        setStepSizeMax(stepSizeMax);
        setStepSizeMin(0.0);
        setStepSize(stepSize);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        this.fixOverlap = fixOverlap;

        throw new IllegalStateException("MCMoveReptate is broken.  Abort!  Abort!");
    }
    
    
    
    /**
     * method to perform a trial move.
     */
    public boolean doTrial(){
       atom = atomSource.getMolecule();
       if(atom == null) return false;
       energyMeter.setTarget(atom);
       uOld = energyMeter.getDataAsScalar();
//       if(uOld > 1e10 && !fixOverlap) {
//           throw new RuntimeException(new ConfigurationOverlapException(atom.node.parentBox()));
//       }
       
//       //This is the part where we make the change.  We've selected an atom already.
//       AtomIteratorBasis aim = new AtomIteratorBasis();
//       aim.setBasis(atom.node.parentMolecule());
//       aim.reset();
//       
//       
//       //Decide which direction the molecule is moving.
//       double dir = Simulation.random.nextDouble();
//       Atom at = new Atom(box.space());
//       Vector store = box.space().makeVector();
//       
//       
//       int leng = aim.size();
//       
//       if(dir <= 0.5){
//           
//       } else {
//           //make sure the vectors are pointing in the proper direction (original)
//           if(conf.isReversed()) {conf.reverse();}
//           
//           at = aim.nextAtom();
//           
//       }
       
       //Pick direction & set up list of atoms to iterate
       forward = random.nextInt(2) == 0;
       IAtomList childlist = atom.getChildList();
       int numChildren = childlist.getAtomCount();
       
       if(forward){
           IVectorMutable position = childlist.getAtom(numChildren-1).getPosition();
           positionOld.E(position);
           for (int j = numChildren - 1; j > 0; j--) {
               IVectorMutable position2 = childlist.getAtom(j-1).getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere(random);
           tempV.TE(bondLength);
           childlist.getAtom(0).getPosition().PE(tempV);
       }
       else {
           IVectorMutable position = childlist.getAtom(0).getPosition();
           positionOld.E(position);
           for(int j = 0; j < numChildren-1; j++){
               IVectorMutable position2 = childlist.getAtom(j+1).getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere(random);
           tempV.TE(bondLength);
           childlist.getAtom(numChildren - 1).getPosition().PE(tempV);
           
       }
       
       uNew = energyMeter.getDataAsScalar();
       return true;
    }

    public double getA(){
        return 1.0;
    }
    
    public double getB(){
        uNew = energyMeter.getDataAsScalar();
        return -(uNew - uOld);
    }
    
    public void acceptNotify(){
        //we don't actually need to do anything here.
    }
    
    public void rejectNotify(){
        IAtomList childlist = atom.getChildList();
        int numChildren = childlist.getAtomCount();
        if (!forward) {
            IVectorMutable position = childlist.getAtom(numChildren-1).getPosition();
            for (int j=numChildren-1; j>0; j--) {
                IVectorMutable position2 = childlist.getAtom(j-1).getPosition();
                position.E(position2);
                position = position2;
            }
            childlist.getAtom(0).getPosition().E(positionOld);
        }
        else {
            IVectorMutable position = childlist.getAtom(0).getPosition();
            for (int j=0; j<numChildren-1; j++) {
                IVectorMutable position2 = childlist.getAtom(j+1).getPosition();
                position.E(position2);
                position = position2;
            }
            childlist.getAtom(numChildren-1).getPosition().E(positionOld);
        }
        
        
    }
    
    public AtomIterator affectedAtoms(){
        return null;
    }
    
    public double energyChange(){return uNew - uOld;}
     
    public void setBox(IBox p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
        tempV = (IVectorRandom)space.makeVector();
        positionOld = space.makeVector();
    }
    
    /**
     * @return Returns the atomSource.
     */
    public MoleculeSource getAtomSource() {
        return atomSource;
    }

    /**
     * @param source The atomSource to set.
     */
    public void setAtomSource(MoleculeSource source) {
        atomSource = source;
    }
    
    
    
    
    
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
//    protected final Vector translationVector;
    private IMolecule atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected MoleculeSource atomSource;
    protected boolean fixOverlap;
    private IVectorRandom tempV;
    private IVectorMutable positionOld;
    private boolean forward;
    private double bondLength;
    protected final IRandom random;
    private final ISpace space;
    
    
}
