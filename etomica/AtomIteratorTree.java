package etomica;

/**
 * Atom iterator that traverses all atoms at a specified depth below a
 * specified root atom in the atom tree hierarchy.  The depth may be
 * specified as any non-negative integer.  If the bottom of the hierarchy
 * is reached before the specified depth, the leaf atoms encountered
 * there are the iterates.
 *
 * @author David Kofke
 * 02.02.16
 */
 
 /* History of changes
  * 07/26/02 started history of changes
  * 07/26/02 (DAK/DW) fixed error in setBasis that left iterator null if given atom is null.
  * 08/24/02 (DAK) revised with overhaul of iterators
  */
  
public class AtomIteratorTree implements AtomIterator {
    
    public AtomIteratorTree() {
    	this(Integer.MAX_VALUE);
    }

	public AtomIteratorTree(int d) {
		setIterationDepth(d);
	}
    /**
     * Returns a new iterator to loop over all the atoms at the
     * given depth below the given atom.  Must call reset() before using.
     */
    public AtomIteratorTree(Atom atom, int d) {
        this();
        setIterationDepth(d);
        setRoot(atom);
    }

    /**
     * Performs action on all iterates.
     */
    public void allAtoms(AtomsetActive act) {
//        if(rootNode == null) return;
//        if(iterationDepth == 0 || rootNode.isLeaf()) act.actionPerformed(rootNode.atom);
//        else if(!rootNode.childrenAreGroups() || iterationDepth == 1) {
//			AtomIteratorListSimple.allAtoms(act, rootNode.childList);
    	if(!doTreeIteration) {
    		iterator.allAtoms(act);
		} else {
			AtomIteratorTree localIterator = new AtomIteratorTree(iterationDepth-1);
			AtomLinker.Tab header = rootNode.childList.header;
			for(AtomLinker e=header.next; e!=header; e=e.next) {
				if(e.atom != null) {
					localIterator.setRoot(e.atom);
					localIterator.allAtoms(act);
				}
			}//end for
		}//end else
    }

    /**
     * Indicates if the iterate has another atom to give.
     */
    public boolean hasNext() {return next != null;}
    
    
    /**
     * Puts iterator in state in which hasNext is false.
     */
    public void unset() {next = null;}
    
    /**
     * Returns true if the given atom is among the iterates as
     * currently configured.
     */
    public boolean contains(Atom[] atom) {
        if(rootNode == null || atom == null || atom[0] == null) return false;
        if(iterationDepth == 0 || rootNode.isLeaf()) return rootNode.atom.equals(atom);
        if(!doTreeIteration) return listIterator.contains(atom);
        
        //eliminated simpler cases; look for atom via explicit iteration
		AtomsetActiveDetect detector = new AtomsetActiveDetect(atom);
		allAtoms(detector);
		return detector.detectedAtom();
    }
    
    
    /**
     * Reinitializes the iterator according to the most recently specified basis
     * and iteration depth.
     */
    public void reset() {
        if(!doTreeIteration) {
        	iterator.reset();
        } else {
            listIterator.reset();
            treeIterator.unset();
            while(listIterator.hasNext() && !treeIterator.hasNext()) {
                treeIterator.setRoot(listIterator.nextAtom()); 
                treeIterator.reset();
            }
        }
        next = iterator.hasNext() ? iterator.nextAtom() : null;
    }
        
    /**
     * Returns the next atom in the iteration sequence.
     */
    public Atom nextAtom() {
    	if(next == null) return null;
        Atom nextAtom = next;
        next = null;
        if(iterator.hasNext()) next = iterator.nextAtom();
        else if(doTreeIteration) {
            while(listIterator.hasNext() && !treeIterator.hasNext()) {
                treeIterator.setRoot(listIterator.nextAtom()); 
                treeIterator.reset();
            }
            if(treeIterator.hasNext()) next = treeIterator.nextAtom();
        }
        return nextAtom;
    }
    
    public Atom[] next() {
    	atoms[0] = nextAtom();
    	return atoms;
    }
    
    public Atom[] peek() {
    	atoms[0] = next;
    	return atoms;
    }
        
    /**
     * Defines the root of the tree under which the iteration is performed.
     * If atom is null, hasNext will report false.
     * User must perform a subsequent call to reset() before beginning iteration.
     */
    public void setRoot(Atom atom) {
        if(atom == null) {
        	rootNode = null;
        	doTreeIteration = false;
        	iterator = AtomIterator.NULL;
        } else if(iterationDepth == 0 || atom.node.isLeaf()) {//singlet iteration of basis atom
	        rootNode = (AtomTreeNodeGroup)atom.node;
            doTreeIteration = false;
            iterator = new AtomIteratorSinglet(atom);
            return;
        } else {
	        rootNode = (AtomTreeNodeGroup)atom.node;
	        listIterator.setList(rootNode.childList);
	        doTreeIteration = (iterationDepth > 1 &&
	                            rootNode.childrenAreGroups());
	        if(doTreeIteration) {
	            if(treeIterator == null) treeIterator = new AtomIteratorTree(iterationDepth-1);
	            iterator = treeIterator;
	        } else {
	            iterator = listIterator;
	        }
        }
        unset();
    }
        
    /**
     * Returns the number of iterates given by a full cycle of this iterator.
     */
    public int size() {
        if(rootNode == null) return 0;
        if(iterationDepth == 0 || rootNode.isLeaf()) return 1;
        if(!doTreeIteration) return listIterator.size();
        
		AtomsetActiveCount counter = new AtomsetActiveCount();
		allAtoms(counter);
		return counter.callCount();
    }
    
    public final int nBody() {return 1;}
    
    /**
     * Sets the depth below the current basis for which the iteration will occur.
     * Any non-negative value is permitted.  A value of zero causes singlet iteration
     * returning just the root atom. A value of 1 returns all children of the basis
     * atom, a value of 2 returns all children of all children of the basis, etc.
     * Iterator returns only atoms at the specified level, and not those above it.
     * Returns atoms at bottom of hierarchy (i.e., leafs) if specified depth
     * exceeds depth of hierarchy.  Default is Integer.MAX_VALUE, which causes
     * all leaf atoms to be iterated.
     */
    public void setIterationDepth(int depth) {
        if (iterationDepth == depth) return;
        if(depth < 0) throw new IllegalArgumentException("Error: attempt to set iteration depth to negative value in AtomIteratorTree");
        iterationDepth = depth;
        if(rootNode != null) setRoot(rootNode.atom);
    }
    /**
     * Returns the currently set value of iteration depth.
     */
    public int getIterationDepth() {return iterationDepth;}
    
    /**
     * Sets iterator to iterate over all leaf atoms below the basis atom.
     * Equivalent to setIterationDepth(Integer.MAX_VALUE).
     */
    public void setAsLeafIterator() {
        setIterationDepth(Integer.MAX_VALUE);
    }
            
    private AtomTreeNodeGroup rootNode;
    private boolean doTreeIteration;
    private AtomIteratorTree treeIterator;
    private final AtomIteratorListSimple listIterator = new AtomIteratorListSimple();
    private AtomIterator iterator;
    private int iterationDepth = Integer.MAX_VALUE;
    private Atom next;
    private Atom[] atoms = new Atom[1];
    
    /**
     * main method to test and demonstrate use of this class.
     */
    public static void main(String args[]) {
        
        Simulation sim = new Simulation();
        Species species2 = new SpeciesSpheresMono();
        Species species1 = new SpeciesSpheres(3,3);
        Species species0 = new SpeciesSpheres(3,2);
        species0.setNMolecules(3);
        species1.setNMolecules(2);
        species2.setNMolecules(2);
        Phase phase = new Phase();
        sim.elementCoordinator.go();
        
        int k = 0;
        AtomIteratorTree treeIterator = new AtomIteratorTree();
 //       treeIterator.setIterationDepth(2);
        treeIterator.setRoot(phase.speciesMaster);
        treeIterator.reset(); k = 0;
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println(treeIterator.size() + "  " + (treeIterator.size() == k));
        System.out.println();
        
        treeIterator.setIterationDepth(2);
        treeIterator.reset(); k = 0;
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println(treeIterator.size() + "  " + (treeIterator.size() == k));
        System.out.println();
        
        treeIterator.setIterationDepth(1);
        treeIterator.reset(); k = 0;
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println(treeIterator.size() + "  " + (treeIterator.size() == k));
        System.out.println();
        
        treeIterator.setIterationDepth(0);
        treeIterator.reset(); k = 0;
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println(treeIterator.size() + "  " + (treeIterator.size() == k));
        System.out.println();
        
        treeIterator.setRoot(((AtomTreeNodeGroup)phase.speciesMaster.node).firstChildAtom());
        treeIterator.setAsLeafIterator();
        treeIterator.setIterationDepth(1);
        treeIterator.reset(); k = 0;
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println(treeIterator.size() + "  " + (treeIterator.size() == k));
        System.out.println();
        
        System.out.print("null-basis test ");
        treeIterator.setRoot(null);
        treeIterator.setAsLeafIterator();
        treeIterator.reset();
        while(treeIterator.hasNext()) System.out.println(k++ + "  " + treeIterator.next().toString());
        System.out.println("ok");
        
    }//end of main
    
}//end of AtomIteratorTree