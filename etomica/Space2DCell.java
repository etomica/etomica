package simulate;
import java.awt.Graphics;
import java.util.Random;

public class Space2DCell extends Space2D implements Space.NeighborIterator {
    
    public Iterator makeIterator(Phase p) {return new NeighborIterator(p);}

    public Space.Coordinate makeCoordinate(Space.Occupant o) {return new Coordinate(o);}

    public Space.Vector makeVector() {return new Vector();}
    public Potential makePotential() {return new CellPotential();}
    public Space.Boundary makeBoundary(int b) {
        switch(b) {
            case(Space2D.Boundary.NONE):            return new Space2D.BoundaryNone();
            case(Space2D.Boundary.PERIODIC_SQUARE): return new Space2D.BoundaryPeriodicSquare();
            default:                                return null;
        }
    }
       
    public static class Vector extends Space2D.Vector {
    }
    
    public class Coordinate extends Space2D.Coordinate implements Lattice.Occupant {
        Coordinate nextNeighbor, previousNeighbor;
        public LatticeSquare.Site cell;        //ok to work with Site rather than Cell (Cell was needed only to make up neighbor-cell list)
        public Coordinate(Space.Occupant o) {super(o);}
        public final void setNextNeighbor(Coordinate c) {
            nextNeighbor = c;
            if(c != null) {c.previousNeighbor = this;}
        }
        public final Lattice.Site site() {return cell;}
        public final void clearPreviousNeighbor() {previousNeighbor = null;}
        public final Coordinate nextNeighbor() {return nextNeighbor;}
        public final Coordinate previousNeighbor() {return previousNeighbor;}
   //Determines appropriate cell and assigns it
        public void assignCell() {             //assumes coordinates ranges [0,1)
            LatticeSquare cells = ((NeighborIterator)parentPhase().iterator()).cells;
            int ix = (int)Math.floor(r.x * cells.dimensions()[0]);
            int iy = (int)Math.floor(r.y * cells.dimensions()[1]);
//            System.out.println(ix+"  "+iy);
            LatticeSquare.Site newCell = cells.sites()[ix][iy];
            if(newCell != cell) {assignCell(newCell);}
            //debugging code
//            LatticeSquare.Cell cellCell = (LatticeSquare.Cell)newCell;
//            System.out.println(cellCell.vertices()[2][0] + "  " + r.x + "  " + cellCell.vertices()[0][0]);
//            System.out.println(cellCell.vertices()[2][1] + "  " + r.y + "  " + cellCell.vertices()[0][1]);
//            System.out.println();
        }
   //Assigns atom to given cell; if removed from another cell, repairs tear in list
        public void assignCell(LatticeSquare.Site newCell) {
            if(previousNeighbor != null) {previousNeighbor.setNextNeighbor(nextNeighbor);}
            else {if(cell != null) cell.setFirst(nextNeighbor); if(nextNeighbor != null) nextNeighbor.clearPreviousNeighbor();}   //removing first atom in cell
            cell = newCell;
            setNextNeighbor((Space2DCell.Coordinate)cell.first());
            cell.setFirst(this);
            clearPreviousNeighbor();
            ((Atom)parent).setColor(cell.color);
        }
    } 
    
    static class Link {
        public Link next;
        public Coordinate coordinate;
    }

    public static final class NeighborIterator extends Iterator {
        public LatticeSquare cells;  //want to declare final, but won't compile
        private int xCells, yCells;
        private double neighborRadius;
        
        public NeighborIterator(Phase p) {
            super(p);
            xCells = yCells = 10;
            cells = new LatticeSquare(LatticeSquare.Cell.class, new int[] {xCells,yCells});
        }
        
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

//        public final AtomPair.Iterator.A makeAtomPairIteratorFull() {return new AtomPairIteratorFull(phase, cells.origin);}
        public final AtomPair.Iterator.A makeAtomPairIteratorUp() {return new AtomPairIteratorUp(phase, cells.origin);}
        public final AtomPair.Iterator.A makeAtomPairIteratorDown() {return new AtomPairIteratorDown(phase, cells.origin);}
        public final Atom.Iterator makeAtomIteratorUp() {return new AtomIteratorUp(cells.origin);}
        public final Atom.Iterator makeAtomIteratorDown() {return new AtomIteratorDown(cells.origin);}

        private static final class AtomPairIteratorUp implements simulate.AtomPair.Iterator.A {
            private final Phase phase;
            final AtomPair pair;
            final Space2D.CoordinatePair cPair;
            Atom atom;
            private boolean hasNext;
            private final LatticeSquare.Site origin;
            private Coordinate coordinate, neighborCoordinate;
            private LatticeSquare.Site neighborCell;
            private LatticeSquare.Linker nextLinker;
            public AtomPairIteratorUp(Phase p, LatticeSquare.Site o) {
                phase = p;
                origin = o;
                pair = new AtomPair(p);
                cPair = (Space2D.CoordinatePair)phase.space().makeCoordinatePair(p.boundary());
                hasNext = false;
            }
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {reset(atom,null,null,null);}
            public void reset(Atom a, Atom d1, Atom d2) {reset(a, d1, d2, null);}
            public void reset(Atom a, Atom d1, Atom d2, Atom d3) {}
            public void reset(Atom a, boolean intra) {  //presently ignores intra, acting as if it were true
                atom = a;
                coordinate = (Coordinate)a.coordinate;
                cPair.c1 = coordinate;
                nextLinker = coordinate.cell.firstUpNeighbor();  //this points to the cell after the current neighbor cell
    //            neighborAtom = coordinate;
    //            hasNext = true;
               neighborCoordinate = coordinate.nextNeighbor();   //this is the next atom to be paired with the fixed atom
               hasNext = true;
               if(neighborCoordinate == null) {advanceCell();}  //sets hasNext to false if can't find neighbor
            }
            // Finds first neighbor of next occupied cell
            private void advanceCell() {
                while(neighborCoordinate == null) {
                    if(nextLinker == null) {hasNext = false; return;} //no more neighbor cells; no upNeighbors for atom
                    neighborCell = (LatticeSquare.Site)nextLinker.site();  //don't need to cast all the way to cell
                    neighborCoordinate = (Coordinate)neighborCell.first;   //get first atom of another neighbor cell; this is null if cell is empty
                    nextLinker = nextLinker.next();
                }
            }
            public AtomPair next() {
                cPair.c2 = neighborCoordinate;
                cPair.reset();
                pair.reset(atom,(Atom)neighborCoordinate.parent(),cPair);  //maybe put an Atom in Coordinate to avoid cast
                neighborCoordinate = neighborCoordinate.nextNeighbor;
                if(neighborCoordinate == null) {advanceCell();}
                return pair;
            }
        } //end of AtomPairIteratorUp

        private static final class AtomPairIteratorDown implements simulate.AtomPair.Iterator.A {
            final AtomPair pair;
            private final Phase phase;
            private final LatticeSquare.Site origin;  //not used correctly down
            final Space2D.CoordinatePair cPair;
            Atom atom;
            private boolean hasNext;
            private Coordinate coordinate, neighborCoordinate;
            private LatticeSquare.Site neighborCell;
            private LatticeSquare.Linker nextLinker;
            private boolean firstCell;
            public AtomPairIteratorDown(Phase p, LatticeSquare.Site o) {
                phase = p; 
                origin = o;
                pair = new AtomPair(p);
                cPair = (Space2D.CoordinatePair)phase.space().makeCoordinatePair(p.boundary());
                hasNext = false;
            }
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {reset(atom,null,null,null);}
            public void reset(Atom a, Atom d1, Atom d2) {reset(a, d1, d2, null);}
            public void reset(Atom a, Atom d1, Atom d2, Atom d3) {}
            public void reset(Atom a, boolean intra) {  //presently ignores intra, acting as if it were true
                atom = a;
                coordinate = (Coordinate)a.coordinate;
                cPair.c1 = coordinate;
                // changes from AtomPairIteratorUp
                nextLinker = coordinate.cell.firstDownNeighbor();  //this points to the cell after the current neighbor cell
                neighborCoordinate = coordinate.previousNeighbor();   //this is the next atom to be paired with the fixed atom
                //////
    //            neighborAtom = coordinate;
    //            hasNext = true;
               firstCell = true;
               hasNext = true;
               if(neighborCoordinate == null) {advanceCell();}  //sets hasNext to false if can't find neighbor
            }
            // Finds first neighbor of next occupied cell      *** CHANGES NEEDED SOMEWHERE IN HERE ***
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
                pair.reset(atom,(Atom)neighborCoordinate.parent(),cPair);  //maybe put an Atom in Coordinate to avoid cast
                neighborCoordinate = firstCell ? neighborCoordinate.previousNeighbor : neighborCoordinate.nextNeighbor;  //first cell advances backwards, subsequent cells (before the first) can advance forward
                if(neighborCoordinate == null) {advanceCell();}
                return pair;
            }
        } //end of AtomPairIteratorDown

        private static final class AtomIteratorUp implements Atom.Iterator {
            private Coordinate coordinate;
            private Atom nextAtom, firstAtom;
            private final LatticeSquare.Site origin;
            private boolean hasNext;
            private LatticeSquare.Site cell;
            public AtomIteratorUp(LatticeSquare.Site o) {hasNext = false;  origin = o;}
            public boolean hasNext() {return hasNext;}
            public void allDone() {hasNext = false;}
            public void reset() {  //resets to first atom in list
                hasNext = true;
                cell = origin;  //rewrite this to use cells.firstOccupant
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
            private final LatticeSquare.Site origin;  //not used correctly down
            private boolean hasNext;
            private LatticeSquare.Site cell, neighborCell;
            private LatticeSquare.Linker nextLinker;
            private boolean firstCell;
            public AtomIteratorDown(LatticeSquare.Site o) {hasNext = false; origin = o;}
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

        
        //These iterators are identical in every Space class; they are repeated in each
        //because they make direct use of the Coordinate type in the class; otherwise casting would be needed
        // Perhaps interitance would work, but haven't tried it
            
        //"Full" --> Each iteration of inner loop begins with same first atom
 /*       private static final class AtomPairIteratorFull implements simulate.AtomPair.Iterator.A {
            final AtomPair pair;
            Atom outer, inner;
            private Atom iFirst, iLast, oFirst, oLast;
            private boolean hasNext;
            public AtomPairIteratorFull(Space.Boundary b) {  //null constructor
                pair = new AtomPair((Boundary)b);
                hasNext = false;
            }  
            public AtomPairIteratorFull(Space.Boundary b, Atom iF, Atom iL, Atom oF, Atom oL) {  //constructor
                pair = new AtomPair((Boundary)b);
                reset(iF,iL,oF,oL);
            }
            public void reset(Atom iL, Atom oF, Atom oL) {reset(oF,iL,oF,oL);}  //take inner and outer first atoms as same
            public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {
                if(iF == null || oF == null) {hasNext = false; return;}
                iFirst = iF; 
                iLast =  iL; 
                oFirst = oF; 
                oLast =  oL;
                inner = iFirst;
                outer = oFirst;
                hasNext = true;
            }
            public simulate.AtomPair next() {
                if(!hasNext) {return null;}
                pair.reset(outer, inner);  //atom1 should always be outer
                if(inner == iLast) {                                     //end of inner loop
                    if(outer == oLast) {hasNext = false;}                //all done
                    else {outer = outer.nextAtom(); inner = iFirst;} //advance outer, reset inner
                }
                else {inner = inner.nextAtom();}
                return pair;
            }
            public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
            public boolean hasNext() {return hasNext;}
            public void reset() {reset(iFirst, iLast, oFirst, oLast);}
        }
            
        //"Half" --> Each iteration of inner loop begins with atom after outer loop atom
/*        private static final class PairIteratorHalf implements simulate.AtomPair.Iterator.A {
            final AtomPair pair;
            AtomCoordinate outer, inner;
            private AtomCoordinate iFirst, iLast, oFirst, oLast;
            private boolean hasNext;
            public PairIteratorHalf(Space.Boundary b) {
                pair = new AtomPair((Boundary)b);
                hasNext = false;
            }
            public PairIteratorHalf(Space.Boundary b, Atom iL, Atom oF, Atom oL) {  //constructor
                pair = new AtomPair((Boundary)b);
                reset(iL,oF,oL);
            }
            public void reset(Atom iF, Atom iL, Atom oF, Atom oL) {reset(iL,oF,oL);} //ignore first argument
            public void reset(Atom iL, Atom oF, Atom oL) {
                if(oF == null) {hasNext = false; return;}
                iLast =  iL; 
                oFirst =  oF; 
                oLast =  oL;
                outer = oFirst;
                inner = outer.nextAtom();
                hasNext = (inner != null);
            }
            public simulate.AtomPair next() {
                pair.c1 = outer;   //c1 should always be outer
                pair.c2 = inner;
                pair.reset();
                if(inner == iLast) {                                     //end of inner loop
                    if(outer == oLast) {hasNext = false;}                //all done
                    else {outer = outer.nextCoordinate; inner = outer.nextCoordinate;} //advance outer, reset inner
                }
                else {inner = inner.nextCoordinate;}
                return pair;
            }
            public final void allDone() {hasNext = false;}   //for forcing iterator to indicate it has no more pairs
            public boolean hasNext() {return hasNext;}
            public void reset() {reset(iLast.atom(), oFirst.atom(), oLast.atom());}
        }
    }    
*/
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