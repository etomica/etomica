package simulate;
import java.awt.Graphics;
import java.util.Random;

public class Space2DCell extends Space2D implements Iterator.Maker {
    
    public Iterator makeIterator(Phase p) {return new NeighborIterator(p);}

    public Space.Coordinate makeCoordinate(Space.Occupant o) {return new Coordinate(o);}

    public Potential makePotential() {return new CellPotential();}
       
    public class Coordinate extends Space2D.Coordinate implements Lattice.Occupant {
        Coordinate nextNeighbor, previousNeighbor;  //next and previous coordinates in neighbor list
        public LatticeSquare.Site cell;             //cell currently occupied by this coordinate; ok to work with Site rather than Cell (Cell was needed only to make up neighbor-cell list)
        public Coordinate(Space.Occupant o) {super(o);}
        public final void setNextNeighbor(Coordinate c) {
            nextNeighbor = c;
            if(c != null) {c.previousNeighbor = this;}
        }
        public final Lattice.Site site() {return cell;}   //Lattice.Occupant interface method
        public final void clearPreviousNeighbor() {previousNeighbor = null;}
        public final Coordinate nextNeighbor() {return nextNeighbor;}
        public final Coordinate previousNeighbor() {return previousNeighbor;}

   //Determines appropriate cell and assigns it
        public void assignCell() {             
            LatticeSquare cells = ((NeighborIterator)parentPhase().iterator()).cells;
            LatticeSquare.Site newCell = cells.nearestSite(this.r, (Space2D.Vector)parentPhase().dimensions());
            if(newCell != cell) {assignCell(newCell);}
        }
   //Assigns atom to given cell; if removed from another cell, repairs tear in list
        public void assignCell(LatticeSquare.Site newCell) {
            if(previousNeighbor != null) {previousNeighbor.setNextNeighbor(nextNeighbor);}
            else {if(cell != null) cell.setFirst(nextNeighbor); if(nextNeighbor != null) nextNeighbor.clearPreviousNeighbor();}   //removing first atom in cell
            cell = newCell;
            if(newCell == null) return;
            setNextNeighbor((Space2DCell.Coordinate)cell.first());
            cell.setFirst(this);
            clearPreviousNeighbor();
//            ((Atom)parent).setColor(cell.color);  //move this to a colorscheme
        }
    } 
    
    /**
     * Iterator of atoms based on a cell-list organization of neighbors
     */
    public static final class NeighborIterator extends Iterator {
        public LatticeSquare cells;  //want to declare final, but won't compile
        private int xCells, yCells;
        private double neighborRadius;
        
        public NeighborIterator(Phase p) {
            super(p);
            xCells = yCells = 10;
            cells = new LatticeSquare(LatticeSquare.Cell.class, new int[] {xCells,yCells});
        }
        
        public void moveNotify(Atom a) {
            ((Coordinate)a.coordinate).assignCell();
        }
        
        public void addMolecule(Molecule m) {
            m.atomIterator.reset();
            while(m.atomIterator.hasNext()) {
                ((Coordinate)m.atomIterator.next().coordinate).assignCell();
            }
        }
        
        public void deleteMolecule(Molecule m) {
            m.atomIterator.reset();
            while(m.atomIterator.hasNext()) {
                ((Coordinate)m.atomIterator.next().coordinate).assignCell(null);
            }
        }
        
        /**
         * Clears cells of all atoms and reassigns all atoms to their cells
         */
        public void reset() {
            cells.clearOccupants();
            for(Atom a=phase.firstAtom(); a!=null; a=a.nextAtom()) {
                Coordinate c = (Coordinate)a.coordinate();
                c.cell = null;
                c.clearPreviousNeighbor();
                c.setNextNeighbor(null);
                ((Coordinate)a.coordinate()).assignCell();
            }
        }
        
        public void setNeighborRadius(double radius) {
            neighborRadius = radius; 
            cells.setNeighborIndexCutoff(radius*radius);  //probably should do assignCell for each atom
            reset();
        }
        public double getNeighborRadius() {return neighborRadius;}

        public final AtomPair.Iterator.A makeAtomPairIteratorUp() {return new AtomPairIteratorUp(phase, cells);}
        public final AtomPair.Iterator.A makeAtomPairIteratorDown() {return new AtomPairIteratorDown(phase, cells);}
        public final Atom.Iterator makeAtomIteratorUp() {return new AtomIteratorUp(cells);}
        public final Atom.Iterator makeAtomIteratorDown() {return new AtomIteratorDown(cells);}

        /**
         * Iterator for atom pairs, going up cell-based neighbor list
         * Takes a given atom and generates pairs from it and uplist neighbors
         * Handles (untested as of yet) case where given atom is not in phase
         */
        private static final class AtomPairIteratorUp implements simulate.AtomPair.Iterator.A {
            private final Phase phase;
            private final LatticeSquare cells;
            final AtomPair pair;
            final Space2D.CoordinatePair cPair;
            Atom atom;
            private boolean hasNext;
            private Coordinate neighborCoordinate;
            private LatticeSquare.Site neighborCell;
            private LatticeSquare.Linker nextLinker;
            public AtomPairIteratorUp(Phase p, LatticeSquare c) {
                phase = p;
                cells = c;
                pair = new AtomPair(p);
                cPair = (Space2D.CoordinatePair)pair.cPair;
                hasNext = false;
            }
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {reset(atom,true);}   //change if INTRA is used
            public void reset(Atom a, boolean intra) {  //presently ignores intra, acting as if it were true
                atom = a;
                if(a.parentPhase() == phase) {
                    Coordinate c = (Coordinate)a.coordinate;
                    cPair.c1 = c;
                    nextLinker = c.cell.firstUpNeighbor();  //this points to the cell after the current neighbor cell
                    neighborCoordinate = c.nextNeighbor();   //this is the next atom to be paired with the fixed atom
                }
                else {  //this atom is not currently in this phase
                    Coordinate c = (Coordinate)a.coordinate;        
                    cPair.c1 = c;
                    LatticeSquare.Site cell = cells.nearestSite(c.r,(Space2D.Vector)phase.dimensions());  //find cell containing atom
                    nextLinker = cell.firstUpNeighbor();               //next cell up list
                    neighborCoordinate = (Coordinate)cell.first();     //"inserted" atom at beginning of list; all atoms in cell are uplist from it
                }    
                hasNext = true;
                if(neighborCoordinate == null) {advanceCell();}  //sets hasNext to false if can't find neighbor
            }
            // Finds first neighbor of next occupied cell
            private void advanceCell() {
                while(neighborCoordinate == null) {
                    if(nextLinker == null) {hasNext = false; return;} //no more neighbor cells; no upNeighbors for atom
                    neighborCell = (LatticeSquare.Site)nextLinker.site();  //don't need to cast all the way to cell
                    neighborCoordinate = (Coordinate)neighborCell.first();   //get first atom of another neighbor cell; this is null if cell is empty
                    nextLinker = nextLinker.next();
                }
            }
            public AtomPair next() {
                cPair.c2 = neighborCoordinate;
                cPair.reset();
                pair.atom1 = atom;
                pair.atom2 = (Atom)neighborCoordinate.parent();
//                pair.reset(atom,(Atom)neighborCoordinate.parent(),cPair);  //maybe put an Atom in Coordinate to avoid cast
                neighborCoordinate = neighborCoordinate.nextNeighbor;
                if(neighborCoordinate == null) {advanceCell();}
                return pair;
            }
        } //end of AtomPairIteratorUp

        /**
         * Iterator for atom pairs, going down cell-based neighbor list
         * Takes a given atom and generates pairs from it and downlist neighbors
         * Handles (untested as of yet) case where given atom is not in phase
         */
        private static final class AtomPairIteratorDown implements simulate.AtomPair.Iterator.A {
            final AtomPair pair;
            private final Phase phase;
            private final LatticeSquare cells;
            final Space2D.CoordinatePair cPair;
            Atom atom;
            private boolean hasNext;
            private Coordinate coordinate, neighborCoordinate;
            private LatticeSquare.Site neighborCell;
            private LatticeSquare.Linker nextLinker;
            private boolean firstCell;
            public AtomPairIteratorDown(Phase p, LatticeSquare c) {
                phase = p; 
                cells = c;
                pair = new AtomPair(p);
                cPair = (Space2D.CoordinatePair)pair.cPair;
                hasNext = false;
            }
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {reset(atom,true);}  //change if intra is used
            public void reset(Atom a, boolean intra) {  //presently ignores intra, acting as if it were true
                atom = a;
                Coordinate c = (Coordinate)a.coordinate;
                cPair.c1 = c;
                if(a.parentPhase() == phase) {  //this atom is currently in this phase
                    nextLinker = c.cell.firstDownNeighbor();  //this points to the cell after the current neighbor cell
                    neighborCoordinate = c.previousNeighbor();   //this is the next atom to be paired with the fixed atom
                }
                else {  //this atom is not currently in this phase
                    LatticeSquare.Site cell = cells.nearestSite(c.r, (Space2D.Vector)phase.dimensions());  //find cell containing atom
                    nextLinker = cell.firstDownNeighbor();             //next cell up list
                    neighborCoordinate = null;     //"inserted" atom at beginning of cell list; none in cell are downlist from it
                }    
               firstCell = true;
               hasNext = true;
               if(neighborCoordinate == null) {advanceCell();}  //sets hasNext to false if can't find neighbor
            }
            // Finds first neighbor of next occupied cell      *** CHANGES (not?) NEEDED SOMEWHERE IN HERE ***
            private void advanceCell() {
                firstCell = false;
                while(neighborCoordinate == null) {
                    if(nextLinker == null) {hasNext = false; return;} //no more neighbor cells; no downNeighbors for atom
                    neighborCell = (LatticeSquare.Site)nextLinker.site();  //don't need to cast all the way to cell
                    neighborCoordinate = (Coordinate)neighborCell.first;   //get first atom of another neighbor cell; this is null if cell is empty
                    nextLinker = nextLinker.next();
                }
            }
            public AtomPair next() {
                cPair.c2 = neighborCoordinate;
                cPair.reset();
                pair.atom1 = atom;
                pair.atom2 = (Atom)neighborCoordinate.parent();
//                pair.reset(atom,(Atom)neighborCoordinate.parent(),cPair);  //maybe put an Atom in Coordinate to avoid cast
                neighborCoordinate = firstCell ? neighborCoordinate.previousNeighbor : neighborCoordinate.nextNeighbor;  //first cell advances backwards, subsequent cells (before the first) can advance forward
                if(neighborCoordinate == null) {advanceCell();}
                return pair;
            }
        } //end of AtomPairIteratorDown

        private static final class AtomIteratorUp implements Atom.Iterator {
            private Coordinate coordinate;
            private Atom nextAtom, firstAtom;
            private final LatticeSquare cells;
            private boolean hasNext;
            private LatticeSquare.Site cell;
            public AtomIteratorUp(LatticeSquare c) {hasNext = false;  cells = c;}
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {  //resets to first atom in list
                hasNext = true;
                cell = cells.origin;  //rewrite this to use cells.firstOccupant
                coordinate = (Coordinate)cell.first();
                if(coordinate == null) {advanceCell();}
            }
            public void reset(Atom a) {
                firstAtom = a;
                if(a == null) {allDone(); return;}
                coordinate = (Coordinate)a.coordinate();
                cell = coordinate.cell;
                hasNext = true;
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                while(coordinate == null) {
                    cell = cell.nextSite();
                    if(cell == null) {        //no more cells;
                        allDone();
                        return;
                    }
                    coordinate = (Coordinate)cell.first;
                }
            }
            public Atom next() {
                nextAtom = (Atom)coordinate.parent();
                coordinate = coordinate.nextNeighbor;      //next atom in cell
                if(coordinate == null) {advanceCell();}
                return nextAtom;
            }
        } //end of AtomIteratorUp
        
        
        private static final class AtomIteratorDown implements Atom.Iterator {
            private Coordinate coordinate;
            private Atom nextAtom, firstAtom;
            private final LatticeSquare cells;  
            private boolean hasNext;
            private LatticeSquare.Site cell, neighborCell;
            private LatticeSquare.Linker nextLinker;
            private boolean firstCell;
            public AtomIteratorDown(LatticeSquare c) {hasNext = false; cells = c;}
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {System.out.println("Error:  should not be calling reset() in Space2DCell.NeighborIterator.AtomIteratorDown"); 
                                    reset(firstAtom);}
            public void reset(Atom a) {
                firstAtom = a;
                if(a == null) {allDone(); return;}
                coordinate = (Coordinate)a.coordinate();
                cell = coordinate.cell;
                hasNext = true;
                firstCell = true;
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                firstCell = false;
                while(coordinate == null) {
                    cell = cell.previousSite();
                    if(cell == null) {        //no more cells;
                        allDone();
                        return;
                    }
                    coordinate = (Coordinate)cell.first();
                }
            }
            public Atom next() {
                nextAtom = (Atom)coordinate.parent();
                coordinate = (firstCell) ? coordinate.previousNeighbor : coordinate.nextNeighbor; //advance down in atom's cell, advance up in all cells below it (because they begin with first atom in cell, not last)
                if(coordinate == null) {advanceCell();}
                return nextAtom;
            }
        } //end of AtomIteratorDown
    }
    
    
    public static class CellPotential extends Potential implements PotentialHard {
        
        public void bump(AtomPair pair) {
            Atom atom = pair.atom1();
            Coordinate atomCoordinate = (Coordinate)atom.coordinate();
            LatticeSquare.Point cell = (LatticeSquare.Point)atomCoordinate.cell;
            double dx = atomCoordinate.r.x - cell.position()[0];
            double dy = atomCoordinate.r.y - cell.position()[1];
            if(Math.abs(dx) > Math.abs(dy)) { //collision with left or right wall
                if(dx > 0) {atomCoordinate.r.x = ((LatticeSquare.Cell)cell.E()).vertices()[2][0]; atomCoordinate.assignCell(cell.E());}
                else {atomCoordinate.r.x = ((LatticeSquare.Cell)cell.W()).vertices()[0][0]; atomCoordinate.assignCell(cell.W());}
            }
            else {
                if(dy > 0) {atomCoordinate.r.y = ((LatticeSquare.Cell)cell.S()).vertices()[2][1]; atomCoordinate.assignCell(cell.S());}
                else {atomCoordinate.r.y = ((LatticeSquare.Cell)cell.N()).vertices()[0][1]; atomCoordinate.assignCell(cell.N());}
            }
        }
        
        public double collisionTime(AtomPair pair) {
            Coordinate atomCoordinate = (Coordinate)pair.atom1().coordinate;
            LatticeSquare.Cell cell = (LatticeSquare.Cell)atomCoordinate.cell;
            double dx = (atomCoordinate.p.x > 0) ? cell.vertices()[0][0] - atomCoordinate.r.x : cell.vertices()[2][0] - atomCoordinate.r.x; 
            double tx = Math.abs(dx/atomCoordinate.p.x);
            double dy = (atomCoordinate.p.y > 0) ? cell.vertices()[0][1] - atomCoordinate.r.y : cell.vertices()[2][1] - atomCoordinate.r.y;
            double ty = Math.abs(dy/atomCoordinate.p.y);
            return (tx > ty) ? ty*atomCoordinate.parent().mass() : tx*atomCoordinate.parent().mass();
//        double tnew = Math.abs((0.5*atom.phase().parentSimulation.space.getNeighborRadius()-0.5*1.0001*((AtomType.Disk)atom.type).diameter())/Math.sqrt(atom.coordinate.momentum().squared()));  //assumes range of potential is .le. diameter
        }
    }
}