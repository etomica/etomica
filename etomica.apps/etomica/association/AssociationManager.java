package etomica.association;
import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IEvent;
import etomica.api.IListener;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;

/**
 * Class to define and track atom associations.  Constructed given an iterator
 * that defines the set of atoms that is managed, and an association definition
 * that is used to determine if atoms are associated.  Only pairwise associations
 * are considered.
 * To incorporate an instance of this class in a simulation, add an instance of 
 * this EnergySum inner class to the simulation:<br>
 *       sim.setEnergySum(associationManager.new EnergySum());
 *
 * If a MC simulation, this should be done before any MCMove classes are instantiated.
 * Also, register this as a listener to IntegratorMC:<br>
 *       integrator.addMCMoveListener(associationManager);
 *
 */
public class AssociationManager implements AgentSource,IListener {
    
    private AssociationDefinition associationDefinition;
    private final IBox box;
    private final AtomLeafAgentManager agentManager;
    private final Api1ACell neighborIterator;
    private final AtomArrayList associatedAtoms = new AtomArrayList();
    
    public AssociationManager(IBox box, PotentialMasterCell potentialMaster, AssociationDefinition definition) {
    	this.box = box;
    	agentManager = new AtomLeafAgentManager(this,box);
        associationDefinition = definition;
        this.neighborIterator = new Api1ACell(3,1.0,potentialMaster.getCellAgentManager());
        }
    
    public void initialize() {
        IAtomList atomList = box.getLeafList();
        for (int i=0; i<atomList.getAtomCount();i+=1) {
        	IAtom atomi = atomList.getAtom(i);
        	for (int j=0; j<atomList.getAtomCount();j+=1) {
            	IAtom atomj = atomList.getAtom(j);
            	if(associationDefinition.isAssociated(atomi,atomj)) {
                    ((AtomArrayList)agentManager.getAgent(atomi)).add(atomj);
                    ((AtomArrayList)agentManager.getAgent(atomj)).add(atomi);
                }
        	}
         }
        for (int i=0; i<atomList.getAtomCount();i+=1) {
        	IAtom atomi = atomList.getAtom(i);
        	if(((AtomArrayList)agentManager.getAgent(atomi)).getAtomCount() > 0){
        		associatedAtoms.add(atomi);
        	}
        }
                
    }
    
    public IAtomList getAssociatedAtoms() {return associatedAtoms;}
    
    //need also to handle associatedAtoms list
    public void actionPerformed(IEvent evt) {
    	MCMoveEvent mcEvent = (MCMoveEvent)evt;
        if(mcEvent instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent) mcEvent).isAccepted()) {
        	return;
            }
        
        AtomIterator iterator = mcEvent.getMCMove().affectedAtoms(box);
        iterator.reset();
        neighborIterator.setBox(box);
        for (IAtom atomi = iterator.nextAtom();atomi != null; atomi =iterator.nextAtom()){
        	AtomArrayList listi = (AtomArrayList)agentManager.getAgent(atomi);
        	neighborIterator.setTarget(atomi);
        	neighborIterator.reset();
        	for (IAtomList atomij = neighborIterator.next();atomij != null; atomij =neighborIterator.next()){
            	IAtom atomj = atomij.getAtom(0);
            	if (atomj == atomi){
            		atomj = atomij.getAtom(1);
            	}
//            	if (atomi.getLeafIndex() == 500 ||atomj.getLeafIndex() == 500 ||atomi.getLeafIndex() == 501 ||atomj.getLeafIndex() == 501){
//            		System.out.println("atomij= "+atomij);
//            	}
            	if (associationDefinition.isAssociated(atomi, atomj)){ //they are associated
            		if (listi.indexOf(atomj) == -1) { // they were not associated
            			if (listi.getAtomCount() == 0) {
            				associatedAtoms.add(atomi);
            			}
            			listi.add(atomj);
            			AtomArrayList listj = (AtomArrayList)agentManager.getAgent(atomj);
            			if (listj.getAtomCount() == 0) {
            				associatedAtoms.add(atomj);
            			}
            			listj.add(atomi);
            		}
            	}
            	else { 
            		int indexj = listi.indexOf(atomj);
            		if ( indexj > -1) {//they were associated
            			listi.remove(indexj);//remove atom j from the listi
            			if ( listi.getAtomCount() == 0) {
            				associatedAtoms.remove(associatedAtoms.indexOf(atomi));
            			}
            			AtomArrayList listj = (AtomArrayList)agentManager.getAgent(atomj);
            			listj.remove(listj.indexOf(atomi));//remove atom i from the listj
            			if ( listj.getAtomCount() == 0) {
            				associatedAtoms.remove(associatedAtoms.indexOf(atomj));
            			}
            		}
            	}
        	}
        }
    }
    
        
    /**
     * Returns the number of atoms on the list of associations of the given atom.
     */
    public IAtomList getAssociatedAtoms(IAtom atom) {
       return (AtomArrayList)agentManager.getAgent(atom);
    }

	public Class getAgentClass() {
		return AtomArrayList.class;
	}

	public Object makeAgent(IAtom a) {
		return new AtomArrayList();
	}

	public void releaseAgent(Object agent, IAtom atom) {
		
	}
}
    