package etomica.models.hexane;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomMolecule;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIterator;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMovePhaseStep;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;

public class MCMoveReptate extends MCMovePhaseStep {
    
    public MCMoveReptate(Simulation sim){
        this(sim.getPotentialMaster(), sim.getDefaults().atomSize, 
                sim.getDefaults().boxSize/2, sim.getDefaults().ignoreOverlap /*, con*/);
    }
    
    public MCMoveReptate(PotentialMaster potentialMaster, double stepSize,
            double stepSizeMax, boolean fixOverlap){
        super(potentialMaster);
        atomSource = new AtomSourceRandomMolecule();
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
       forward = Simulation.random.nextBoolean();
       AtomArrayList childlist = ((AtomTreeNodeGroup)atom.getNode()).getChildList();
       int numChildren = childlist.size();
       
       if(forward){
           IVector position = ((AtomLeaf)childlist.get(numChildren-1)).getCoord().getPosition();
           positionOld.E(position);
           for (int j = numChildren - 1; j > 0; j--) {
               IVector position2 = ((AtomLeaf)childlist.get(j-1)).getCoord().getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere();
           tempV.TE(bondLength);
           ((AtomLeaf)childlist.get(0)).getCoord().getPosition().PE(tempV);
       }
       else {
           IVector position = ((AtomLeaf)childlist.get(0)).getCoord().getPosition();
           positionOld.E(position);
           for(int j = 0; j < numChildren-1; j++){
               IVector position2 = ((AtomLeaf)childlist.get(j+1)).getCoord().getPosition();
               position.E(position2);
               position = position2;
           }
           tempV.setRandomSphere();
           tempV.TE(bondLength);
           ((AtomLeaf)childlist.get(numChildren - 1)).getCoord().getPosition().PE(tempV);
           
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
        AtomArrayList childlist = ((AtomTreeNodeGroup)atom.getNode()).getChildList();
        int numChildren = childlist.size();
        if (!forward) {
            IVector position = ((AtomLeaf)childlist.get(numChildren-1)).getCoord().getPosition();
            for (int j=numChildren-1; j>0; j--) {
                IVector position2 = ((AtomLeaf)childlist.get(j-1)).getCoord().getPosition();
                position.E(position2);
                position = position2;
            }
            ((AtomLeaf)childlist.get(0)).getCoord().getPosition().E(positionOld);
        }
        else {
            IVector position = ((AtomLeaf)childlist.get(0)).getCoord().getPosition();
            for (int j=0; j<numChildren-1; j++) {
                IVector position2 = ((AtomLeaf)childlist.get(j+1)).getCoord().getPosition();
                position.E(position2);
                position = position2;
            }
            ((AtomLeaf)childlist.get(numChildren-1)).getCoord().getPosition().E(positionOld);
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
        tempV = phase.space().makeVector();
        positionOld = phase.space().makeVector();
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
    private Atom atom;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected boolean fixOverlap;
    private IVector tempV;
    private IVector positionOld;
    private boolean forward;
    private double bondLength;
    
    
}
