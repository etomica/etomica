package etomica;

import etomica.action.AtomsetAction;

public class ApiInnerVariableBonded extends ApiInnerVariable {
    
    /**
     * Construct a bond iterator using the given atom iterators.  Requires
     * call to reset() before beginning iteration.
     */
    public ApiInnerVariableBonded(AtomIterator aiOuter, 
                            AtomIteratorAtomDependent aiInner) {
        super(aiOuter,aiInner);
    }
    
    /**
     * Indicates whether the given atom pair will be returned by the
     * iterator during its iteration. The order of the atoms in the pair
     * is significant (this means that a value of true is returned only if
     * one of the pairs returned by the iterator will have the same two 
     * atoms in the same atom1/atom2 position as the given pair). Not dependent
     * on state of hasNext.
     */
    public boolean contains(Atom[] pair) {
        if(aiOuter.contains(new Atom[] {pair[0]})) {
            aiInner.setAtom(pair[0]);
            return aiInner.peek()[0] == pair[1];
        }
        return false;
    }

    /**
     * Returns the number of pairs given by this iterator.  Independent
     * of state of hasNext. Clobbers the iteration state (i.e., status
     * of hasNext/next) but does not recondition iterator (i.e., does not 
     * change set of iterates that would be given on iteration after reset).
     * Must perform reset if attempting iteration after using size() method.
     */
    public int size() {
        aiOuter.reset();
        return aiOuter.size()-1;
    }        
    
    /**
     * Returns the next pair of atoms. The same AtomPair instance
     * is returned every time, but the Atoms it holds are (of course)
     * different for each iterate. 
     */
    public Atom[] next() {
        if(!hasNext) return null;
        //we use this update flag to indicate that atom1 in pair needs to be set to a new value.
        //it is not done directly in the while-loop because pair must first return with the old atom1 intact
        if(needUpdate1) {pair[0] = atom1; needUpdate1 = false;}  //aiOuter was advanced
        pair[1] = aiInner.nextAtom();
        if(aiOuter.hasNext()) {     //Outer has another atom1...
            atom1 = aiOuter.nextAtom();           //...get it
            aiInner.setAtom(atom1);
            aiInner.reset();
            needUpdate1 = true;           //...flag update of pair.atom1 for next time
            while(!aiInner.hasNext()) {     //Inner is done for this atom1, loop until it is prepared for next
                if(aiOuter.hasNext()) {     //Outer has another atom1...
                    atom1 = aiOuter.nextAtom();           //...get it
                    aiInner.setAtom(atom1);
                    aiInner.reset();
                    needUpdate1 = true;           //...flag update of pair.atom1 for next time
                }
                else {hasNext = false; break;} //Outer has no more; all done with pairs
            }//end while
        }
        else {hasNext = false;} //Outer has no more; all done with pairs
        return pair;
    }

    /**
     * Performs the given action on all pairs returned by this iterator.
     */
    public void allAtoms(AtomsetAction act) {  
        aiOuter.reset();
        while(aiOuter.hasNext()) {
            pair[0] = aiOuter.nextAtom();
            aiInner.setAtom(pair[0]);
            aiInner.reset();
            if(aiInner.hasNext()){
                pair[1] = aiInner.nextAtom();
                act.actionPerformed(pair);
            }
        }
    }
    
}
    
