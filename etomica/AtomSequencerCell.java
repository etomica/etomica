package etomica;

public class AtomSequencerCell implements AtomSequencer, AbstractLattice.Occupant {
    
    Atom next, previous;  //next and previous coordinates in neighbor list
    public AtomCell cell;             //cell currently occupied by this coordinate
    public Lattice lattice;               //cell lattice in the phase occupied by this coordinate
    private final Atom atom;
    
    public AtomSequencerCell(Atom a) {
        atom = a;
    }
    public void setNextAtom(Atom a) {
        next = a;
        if(a != null) a.seq.setPreviousAtom(this.atom);
    }
    
    public void setPreviousAtom(Atom a) {
        previous = a;
    }
    public Site site() {return cell;}   //Lattice.Occupant interface method
    public void clearPreviousAtom() {previous = null;}
    public Atom nextAtom() {return next;}
    public Atom previousAtom() {return previous;}

//Determines appropriate cell and assigns it
    public void assignCell() {
        AtomCell newCell = (AtomCell)lattice.nearestSite(atom.position(), dimensions);
        if(newCell != cell) {assignCell(newCell);}
    }
//Assigns atom to given cell; if removed from another cell, repairs tear in list
    public void assignCell(AtomCell newCell) {
        if(previous != null) {previous.setNextAtom(next);}
        else {//removing first atom in cell
            if(cell != null) cell.setFirst(next); 
            if(next != null) next.clearPreviousAtom();
        }   
        cell = newCell;
        if(cell == null) {setNextAtom(null);}
        else {
            setNextAtom(cell.first);
            cell.first = this;
        }
        clearPreviousAtom();
    }//end of assignCell
}//end of AtomSequencerCell
    
