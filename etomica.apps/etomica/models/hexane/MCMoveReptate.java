package etomica.models.hexane;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMovePhaseStep;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.util.IRandom;

public class MCMoveReptate extends MCMovePhaseStep {
    
    public MCMoveReptate(Simulation sim){
        this(sim.getPotentialMaster(), sim.getRandom(), sim.getDefaults().atomSize, 
                sim.getDefaults().boxSize/2, sim.getDefaults().ignoreOverlap /*, con*/);
    }
    
    public MCMoveReptate(PotentialMaster potentialMaster, IRandom random, 
            double stepSize, double stepSizeMax, boolean fixOverlap){
        super(potentialMaster);
        this.random = random;
        atomSource = new AtomSourceRandomMolecule();
        ((AtomSourceRandomMolecule)atomSource).setRandom(random);
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
        setName("MCMoveReptate");
        this.fixOverlap = fixOverlap;

        
    }
    
    
    
    /**
     * method to perform a trial move.
     */
    public boolean doTrial(){
       atom = atomSource.getAtom();
       if(atom == null) return false;
       energyMeter.setTarget(atom);
       uOld = energyMeter.getDataAsScalar();
//       if(uOld > 1e10 && !fixOverlap) {
//           throw new RuntimeException(new ConfigurationOverlapException(atom.node.parentPhase()));
//       }
       
//       //This is the part where we make the change.  We've selected an atom already.
//       AtomIteratorBasis aim = new AtomIteratorBasis();
//       aim.setBasis(atom.node.parentMolecule());
//       aim.reset();
//       
//       
//       //Decide which direction the molecule is moving.
//       double dir = Simulation.random.nextDouble();
//       Atom at = new Atom(phase.space());
//       Vector store = phase.space().makeVector();
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
       AtomArrayList childlist = ((IAtomGroup)atom).getChildList();
       int numChildren = childlist.getAtomCount();
       
       if(forward){
           IVector position = ((IAtomPositioned)childlist.getAtom(numChildren-1)).getPosition();
           positionOld.E(position);
           for (int j = numChildren - 1; j > 0; j--) {
               IVector position2 = ((IAtomPositioned)childlist.getAtom(j-1)).getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere(random);
           tempV.TE(bondLength);
           ((IAtomPositioned)childlist.getAtom(0)).getPosition().PE(tempV);
       }
       else {
           IVector position = ((IAtomPositioned)childlist.getAtom(0)).getPosition();
           positionOld.E(position);
           for(int j = 0; j < numChildren-1; j++){
               IVector position2 = ((IAtomPositioned)childlist.getAtom(j+1)).getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere(random);
           tempV.TE(bondLength);
           ((IAtomPositioned)childlist.getAtom(numChildren - 1)).getPosition().PE(tempV);
           
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
        AtomArrayList childlist = ((IAtomGroup)atom).getChildList();
        int numChildren = childlist.getAtomCount();
        if (!forward) {
            IVector position = ((IAtomPositioned)childlist.getAtom(numChildren-1)).getPosition();
            for (int j=numChildren-1; j>0; j--) {
                IVector position2 = ((IAtomPositioned)childlist.getAtom(j-1)).getPosition();
                position.E(position2);
                position = position2;
            }
            ((IAtomPositioned)childlist.getAtom(0)).getPosition().E(positionOld);
        }
        else {
            IVector position = ((IAtomPositioned)childlist.getAtom(0)).getPosition();
            for (int j=0; j<numChildren-1; j++) {
                IVector position2 = ((IAtomPositioned)childlist.getAtom(j+1)).getPosition();
                position.E(position2);
                position = position2;
            }
            ((IAtomPositioned)childlist.getAtom(numChildren-1)).getPosition().E(positionOld);
        }
        
        
    }
    
    public AtomIterator affectedAtoms(){
        return null;
    }
    
    public double energyChange(){return uNew - uOld;}
     
    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
        atomSource.setPhase(p);
        tempV = (IVectorRandom)phase.getSpace().makeVector();
        positionOld = phase.getSpace().makeVector();
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param atomSource The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }
    
    
    
    
    
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
//    protected final Vector translationVector;
    private IAtom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    private IVectorRandom tempV;
    private IVector positionOld;
    private boolean forward;
    private double bondLength;
    protected final IRandom random;
    
    
}
