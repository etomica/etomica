package etomica.association;
import etomica.*;

/**
 * Manages information about the association of atoms with each other.  
 * Causes a list to be kept with each atom that holds references to other
 * atoms associated with it, and keeps a master list of all atoms that have 
 * at least one associate.  Definition of "association" is given by the 
 * AssociationDefinition class passed to the constructor.
 */
 //needs to be further developed to handle MCMoveEvents when affected atoms
 //are not among those managed by this class
public class AssociationManager implements MCMoveListener, Atom.AgentSource {
    
    private final AssociationDefinition associationDefinition;
    private final int index = Atom.requestAgentIndex(this);
    private AtomListRestorable associationList = new AtomListRestorable();
    private AtomTreeNodeGroup basisNode1, basisNode2;
    private Phase phase;
    private final AtomIteratorListSimple listIterator = new AtomIteratorListSimple();
    private AtomPairIterator pairIterator1, pairIteratorA;
    private final Simulation simulation;
    private final IteratorDirective directive = new IteratorDirective(IteratorDirective.BOTH);
   
    public AssociationManager(Simulation sim, AssociationDefinition definition) {
        simulation = sim;
        associationDefinition = definition;
    }
    
    /**
     * Identifies the parent group(s) of the atoms to be managed by this object.
     */
    public void setBasis(AtomTreeNodeGroup basisNode1, AtomTreeNodeGroup basisNode2) {
        phase = basisNode1.parentPhase();
        if(phase != basisNode2.parentPhase()) throw new IllegalArgumentException("AssociationManager.setBasis:  basis nodes must be in the same phase");
        this.basisNode1 = basisNode1;
        this.basisNode2 = basisNode2;
        if(basisNode1 == basisNode2) {//intragroup association
            pairIterator1 = new ApiIntragroup1A(simulation);
            pairIteratorA = new ApiIntragroupAA(simulation);
        } else {
            pairIterator1 = new ApiIntergroup1A(simulation);
            pairIteratorA = new ApiIntergroupAA(simulation);
        }
        pairIterator1.setBasis(basisNode1.atom(), basisNode2.atom());
        pairIteratorA.setBasis(basisNode1.atom(), basisNode2.atom());
        
        clearLists();
    }
 
    /**
     * Clears association lists of all managed atoms, and clears list
     * of associated atoms.  Also clears all memories.
     */
    public void clearLists() {
        associationList.clear();
        associationList.clearMemory();
        listIterator.setBasis(basisNode1.childList);
        listIterator.reset();
        while(listIterator.hasNext()) {
            Atom atom = listIterator.next();
            ((AtomListRestorable)atom.allatomAgents[index]).clear();
            ((AtomListRestorable)atom.allatomAgents[index]).clearMemory();
        }
        if(basisNode1.equals(basisNode2)) return;
        listIterator.setBasis(basisNode2.childList);
        listIterator.reset();
        while(listIterator.hasNext()) {
            Atom atom = listIterator.next();
            ((AtomListRestorable)atom.allatomAgents[index]).clear();
            ((AtomListRestorable)atom.allatomAgents[index]).clearMemory();
        }
    }//end clearLists        
    
    /**
     * Clears only memory of changes for association lists of all managed atoms and
     * for list of associated atoms.  Lists themselves are unaltered.
     */
    public void clearMemories() {
        associationList.clearMemory();
        listIterator.setBasis(basisNode1.childList);
        listIterator.reset();
        while(listIterator.hasNext()) {
            Atom atom = listIterator.next();
            ((AtomListRestorable)atom.allatomAgents[index]).clearMemory();
        }
        if(basisNode1.equals(basisNode2)) return;
        listIterator.setBasis(basisNode2.childList);
        listIterator.reset();
        while(listIterator.hasNext()) {
            Atom atom = listIterator.next();
            ((AtomListRestorable)atom.allatomAgents[index]).clearMemory();
        }
    }        
    
    /**
     * Initializes all lists based on current configuration.
     */
    public void initialize() {
        clearLists();
        
        //loop over all managed pairs
        pairIteratorA.reset();
        while(pairIteratorA.hasNext()){
            AtomPair pair = pairIteratorA.next();
            Atom atom1 = pair.atom1();
            Atom atom2 = pair.atom2();
            if(atom1 != atom2) {//this check shouldn't be needed
                //see if atom1 and atom2 are associated and update accordingly
                AtomListRestorable atom1List = (AtomListRestorable)atom1.allatomAgents[index];
                AtomListRestorable atom2List = (AtomListRestorable)atom2.allatomAgents[index];
                if(associationDefinition.isAssociated(atom1, atom2)) {
                    if(!(atom1List.contains(atom2))){//this check shouldn't be needed since atom1 and atom2 are encountered only once in iteration
                        atom1List.add(atom2);
                    }
                    if(!(atom2List.contains(atom1))){
                        atom2List.add(atom1);
                    }
                    //add to association list if this is the first addition to each atom's list
                    if(atom1List.size() == 1) associationList.add(atom1);
                    if(atom2List.size() == 1) associationList.add(atom2);
                }
            } else throw new RuntimeException("AssociationManager.initialize encountered pair formed from same atom");
        }//end while
        
        clearMemories();
     }//end initialize
    
    /**
     * Updates lists in response to a MCMove trial event.
     */
    public void actionPerformed(MCMoveEvent evt) {
       
        AtomIterator ai = evt.mcMove.affectedAtoms(phase);
        if(evt.isTrialNotify) {
           //update neighbors
            while(ai.hasNext()) updateList(ai.next());
            
        } else { //notification of acceptance/rejection -- restore or save association lists
            if(evt.wasAccepted) {//accepted, so clear memory of changes
                clearMemories();
            
            } else {//rejected, restore relevant lists
                while(ai.hasNext()) restoreLists(ai.next());
                associationList.restore();
            }
        }
    }//end actionPerformed
    
    /**
     * Restores associations list of the given atom, as well as those of that
     * atoms that were on its association list and are on it before/after
     * its list is restored.
     */
    private void restoreLists(Atom atom){
        //restore list for atoms on atom's current list
        listIterator.setBasis((AtomList)atom.allatomAgents[index]);
        while(listIterator.hasNext()) {
            ((AtomListRestorable)listIterator.next().allatomAgents[index]).restore();
        }
        
        //restore atom's list
        ((AtomListRestorable)atom.allatomAgents[index]).restore();
        
        //restore lists for atoms on atom's restored list
        listIterator.setBasis((AtomList)atom.allatomAgents[index]);
        while(listIterator.hasNext()) {
            ((AtomListRestorable)listIterator.next().allatomAgents[index]).restore();
        }
    }//end restoreLists
    
    /**
     * Returns the total number of atoms with at least one associate.
     */
    public int associatedAtomCount() {
       return associationList.size();
    }
    
    /**
     * Returns the number of atoms in the list of associations of the given atom.
     */
    public int associationCount(Atom atom) {
       return ((AtomList)atom.allatomAgents[index]).size();
    }
    
    /**
     * Returns a randomly selected atom from the list of all associated atoms.
     */
    public Atom randomAssociatedAtom() {
        return associationList.getRandom();
    }
    
    /**
     * Updates the association list of the given atom.  First clears the list while removing
     * atom from other atoms' association lists, then loops through neighbors and places
     * them on list as appropriate.
     */
    public void updateList(Atom atom){
        //clear association list of atom and remove it from associates' list
        AtomList atomList = (AtomList)atom.allatomAgents[index];
        if(atomList.size() > 0) {
            listIterator.setBasis(atomList);
            listIterator.reset();
            while(listIterator.hasNext()) {
                Atom atom2 = listIterator.next();
                AtomListRestorable atom2List = (AtomListRestorable)atom2.allatomAgents[index];
                atomList.remove(atom2);//iterator can handle atom2 removal even in middle of list iteration
                atom2List.remove(atom);
                if(atom2List.size() == 0) associationList.remove(atom2);
            }//end while
            associationList.remove(atom);
        }//end if
        
        //look at neighbors to see if they belong
        pairIterator1.reset(directive.set(atom));//maybe should work with neighbor AtomIterator, but complicated by different bases
        while(pairIterator1.hasNext()){
            AtomPair pair = pairIterator1.next();
            Atom atomA = pair.atom1();
            Atom atomB = pair.atom2();
            AtomListRestorable atomAList = (AtomListRestorable)atomA.allatomAgents[index];
            AtomListRestorable atomBList = (AtomListRestorable)atomB.allatomAgents[index];
            
            if(associationDefinition.isAssociated(atomA, atomB)) {
                if(atomAList.size() == 0) associationList.add(atomA);
                if(atomBList.size() == 0) associationList.add(atomB);
                atomAList.add(atomB);
                atomBList.add(atomA);
            }//end if
        }//end while
    }//end updateList
    
    
    /**
     * Returns a new AtomList as an agent for placement in each atom as it is constructed.
     * Implementation of the Atom.AgentSource interface.
     */
    public Object makeAgent(Atom a) {return new AtomListRestorable();}
    
}//end of AssociationManager

