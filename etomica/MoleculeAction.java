package simulate;

/**
 * 
 * @author David Kofke
 * @see Molecule.Iterator (when it is written)
 */

public abstract class MoleculeAction extends simulate.Action {
    
    protected Molecule molecule;
    public void setMolecule(Molecule m) {molecule = m;}
    public Molecule getMolecule() {return molecule;}
                
    public void actionPerformed() {if(molecule != null) actionPerformed(molecule);}
    
    /**
        * Method called by the iterator, once for each of its molecules.
        * 
        * @param m Molecule passed to method by iterator
        */
    public abstract void actionPerformed(Molecule m);
        
    //***** end of Action methods; begin definition of subclasses *****//

                
    public static class Translate extends MoleculeAction {
        protected Space.Vector displacement;
            
        public Translate(Space space) {
            super();
            displacement = space.makeVector();
        }
            
        public final void actionPerformed(Molecule m) {m.r.PE(displacement);}
        public void actionPerformed(Molecule m, Space.Vector d) {m.r.PE(d);}
        public final void setDisplacement(Space.Vector d) {displacement.E(d);}
    }//end of Translate
    

} //end of MoleculeAction   