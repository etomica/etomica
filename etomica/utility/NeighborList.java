package etomica.utility;
import etomica.*;

import java.util.Arrays;

/**
 * Holds a static method to create an array of nearest-neighbor atoms to each atom
 * given by an iterator.
 * Uses Java 1.2 
 */
 
public class NeighborList {
        
/**
 * Sets up atomList elements in each atom to point to the nearest atoms
 * to it in the current configuration.
 */
    public static void setup(AtomIterator iterator, 
                                Space.Boundary boundary, 
                                int neighborCount,
                                int iUp, int iDn) {
            
        if(neighborCount > iterator.size() || neighborCount < 0) {//need an exception
            System.out.println("Bad argument in NeighborList.setup");
            System.exit(1);
        }
        //make copy of atom iterator for double looping
        Atom[] atoms = new Atom[iterator.size()];
        Atom[] sortArray = new Atom[iterator.size()];
        for(int i=0; i<atoms.length; i++) {
            atoms[i] = iterator.next();
            sortArray[i] = atoms[i];
        }
            
        R2Comparator comparator = new R2Comparator(boundary); 
            
        for(int i=0; i<atoms.length; i++) {//loop over atoms, setting up their neighbor lists
            Atom root = atoms[i];
            comparator.setOriginAtom(root);
            Arrays.sort(sortArray, comparator);
            for(int k=1; k<neighborCount; k++) {//loop over neighborCount nearest atoms (skip first, which is root atom)
                for(int j=0; j<atoms.length; j++) {//look for neighbor downlist of root
                    if(atoms[j] == root) {//uplisted to root without finding neighbor; neighbor is uplist
                        root.atomList[iUp].add(sortArray[k]);
                        break; //out of j loop
                    }
                    if(atoms[j] == sortArray[k]) {//found it downlist of root
                        root.atomList[iDn].add(sortArray[k]);
                        break; //out of j loop
                    }
                }//end for j
            }//end for k
        }//end for i
    }//end of setup
            
    private static class R2Comparator implements java.util.Comparator {
            
        public final Space.Boundary boundary;
        private Atom originAtom;
        private Space.Vector origin;
            
        public R2Comparator(Space.Boundary boundary) {
            this.boundary = boundary;
        }
            
        public void setOriginAtom(Atom a) {
            originAtom = a;
            origin = a.coord.position();
        }
            
        public int compare(Object obj1, Object obj2) {
            Atom a1 = (Atom)obj1;
            Atom a2 = (Atom)obj2;
            double r12 = Space.r2(a1.coord.position(), origin, boundary);
            double r22 = Space.r2(a2.coord.position(), origin, boundary);
            if(r12 < r22) return -1;
            else if(r12 > r22) return +1;
            else return 0;
        }
    }//end of R2Comparator
    
    public static void main(String[] args) {
        
    Simulation sim = new Simulation(new etomica.Space3D());
    Simulation.instance = sim;
    etomica.Phase phase  = new etomica.Phase();
    phase.setConfiguration(new ConfigurationFcc(sim.space));
    etomica.SpeciesSpheresMono speciesSpheres  = new etomica.SpeciesSpheresMono();
    speciesSpheres.setNMolecules(32);
    
    int iUp = Atom.requestAtomListIndex();
    int iDn = Atom.requestAtomListIndex();
    
    sim.mediator().go(); 
    
    AtomIteratorSequential iterator = new AtomIteratorSequential();
    iterator.setBasis(phase.speciesMaster);
    
    NeighborList.setup(iterator, phase.boundary(), 12, iUp, iDn);
    
    //write out something to check that it worked ok
    
  }//end of main
        
   
} //end of NeighborList


