package simulate;
import simulate.lattice.*;
import java.awt.Graphics;
import java.awt.Color;
import java.util.Random;
import java.util.Observer;
import java.util.Observable;

public class Space2DCell extends Space2D implements IteratorFactory.Maker {
    
    public IteratorFactory makeIteratorFactory(Phase p) {return new CellListIteratorFactory(p);}

    public Space.Coordinate makeCoordinate(Space.Occupant o) {return new Coordinate(o);}

//    public Potential makePotential(Phase p) {return new CellPotential(p);}
    
    public Space.Boundary makeBoundary(Space.Boundary.Type t) {
        if(t == Space2D.Boundary.PERIODIC_SQUARE) return new BoundaryPeriodicSquare();  //this space's version overrides centralImage methods
        else return super.makeBoundary(t);
    }

    public final class Coordinate extends Space2D.Coordinate implements AbstractLattice.Occupant {
        Coordinate nextNeighbor, previousNeighbor;  //next and previous coordinates in neighbor list
        public AtomCell cell;             //cell currently occupied by this coordinate
        public SquareLattice lattice;               //cell lattice in the phase occupied by this coordinate
        public Coordinate(Space.Occupant o) {
            super(o);
                //lattice can be reassigned by add/delete molecule of CellListIteratorFactory            
            if(o.parentPhase() != null) lattice = ((CellListIteratorFactory)o.parentPhase().iteratorFactory()).lattice(); 
        }
        public void setNextNeighbor(Coordinate c) {
            nextNeighbor = c;
            if(c != null) {c.previousNeighbor = this;}
        }
        public Site site() {return cell;}   //Lattice.Occupant interface method
        public void clearPreviousNeighbor() {previousNeighbor = null;}
        public Coordinate nextNeighbor() {return nextNeighbor;}
        public Coordinate previousNeighbor() {return previousNeighbor;}
        
        public void setLattice(SquareLattice squareLattice) {
            lattice = squareLattice;
            if(lattice != null) assignCell();
        }

   //Determines appropriate cell and assigns it
        public void assignCell() {             
            AtomCell newCell = (AtomCell)lattice.nearestSite(this.r, (Space2D.Vector)parent.parentPhase().dimensions());
            if(newCell != cell) {assignCell(newCell);}
        }
   //Assigns atom to given cell; if removed from another cell, repairs tear in list
        public void assignCell(AtomCell newCell) {
            if(previousNeighbor != null) {previousNeighbor.setNextNeighbor(nextNeighbor);}
            else {//removing first atom in cell
                if(cell != null) cell.setFirst(nextNeighbor); 
                if(nextNeighbor != null) nextNeighbor.clearPreviousNeighbor();
            }   
            cell = newCell;
            if(cell == null) {setNextNeighbor(null);}
            else {
                setNextNeighbor(cell.first);
                cell.first = this;
            }
            clearPreviousNeighbor();
        }//end of assignCell
    }//end of Space2DCell.Coordinate
    
    //////////// end of Space2DCell methods ////////////
    
    /**
     * Iterator of atoms based on a cell-list organization of neighbors
     */
    public static final class CellListIteratorFactory extends IteratorFactory implements PotentialField.Maker {
        private SquareLattice lattice; 
        private SiteIterator.List latticeIterator;
        private int xCells, yCells;
        private double neighborDistance, absoluteNeighborDistance;  //look at this more carefully if boundary size fluctuates
        
        public CellListIteratorFactory(Phase p) {
            super(p);
            absoluteNeighborDistance = 1.5*Default.ATOM_SIZE;
            setNCells(3,3);
        }
        
        public PotentialField makePotentialField(Phase p) {
            return new CellPotential(p);
        }
        
        public SquareLattice lattice() {return lattice;}
        
        // don't make public until rectangular lattice class is written
        public final void setNCells(int nx, int ny) {
            xCells = nx;
            yCells = ny;
            double cellWidthX = 1.0/(double)nx;
            double cellWidthY = 1.0/(double)ny;
            double halfWidthX = 0.5*cellWidthX;
            double halfWidthY = 0.5*cellWidthY;
            //Construct lattice
            lattice = new SquareLattice(new int[] {xCells, yCells}, new AtomCellFactory(), cellWidthX); 
            lattice.translateBy(new Space2D.Vector(halfWidthX, halfWidthY));
            for(int i=0; i<xCells; i++) {
                for(int j=0; j<yCells; j++) {
                    AtomCell centralCell = (AtomCell)lattice.sites[i][j];
                    centralCell.position = (Space2D.Vector)((BravaisLattice.Coordinate)centralCell.coordinate()).position();
                    centralCell.E = (i+1==xCells) ? (AtomCell)lattice.sites[0][j] : (AtomCell)lattice.sites[i+1][j];
                    centralCell.W = (i==0) ? (AtomCell)lattice.sites[xCells-1][j] : (AtomCell)lattice.sites[i-1][j];
                    centralCell.S = (j+1==yCells) ? (AtomCell)lattice.sites[i][0] : (AtomCell)lattice.sites[i][j+1];
                    centralCell.N = (j==0) ? (AtomCell)lattice.sites[i][yCells-1] : (AtomCell)lattice.sites[i][j-1];
                    centralCell.xE = centralCell.position.x + halfWidthX;
                    centralCell.xW = centralCell.position.x - halfWidthX;
                    centralCell.yN = centralCell.position.y - halfWidthY;//y increases downward
                    centralCell.yS = centralCell.position.y + halfWidthY;
                //    System.out.println("Cell: "+centralCell.coordinate().toString());
                //    System.out.println("E "+centralCell.E.coordinate().toString());
                //    System.out.println("W "+centralCell.W.coordinate().toString());
                //    System.out.println("S "+centralCell.S.coordinate().toString());
                //    System.out.println("N "+centralCell.N.coordinate().toString());
                //    System.out.println();
                }
            }
            //Construct a List iterator that can make Cursors
            latticeIterator = new SiteIterator.List(lattice.iterator());
            //Set up neighbors of each site
            setupNeighbors();
            reset();
        }
        public final void setNCells(int n) {
            setNCells(n,n);
        }
       
        public void setNeighborDistance(double distance) {
            absoluteNeighborDistance = distance;
            if(phase().boundary() != null) {
                neighborDistance = distance/phase().boundary().dimensions().component(0); 
                setupNeighbors();
            }
        }
        public double getNeighborDistance() {return neighborDistance;}
        
        //anonymous member class
        SiteIterator.Neighbor.Criterion neighborCriterion = new SiteIterator.Neighbor.Criterion() {
            Space2D.BoundaryPeriodicSquare boundary = new Space2D.BoundaryPeriodicSquare(1.0,1.0);
            public boolean areNeighbors(Site s1, Site s2) {
                return ((AbstractCell)s1).r2NearestVertex((AbstractCell)s2, boundary) < neighborDistance;
            }
        };
        
        private void setupNeighbors() {
            SiteIterator.Cursor cursor = latticeIterator.makeCursor();
            latticeIterator.reset();
            while(latticeIterator.hasNext()) {
                SiteIterator.Neighbor nbrIterator = latticeIterator.next().adjacentIterator();
                nbrIterator.clearAll();
                nbrIterator.setNeighbors(cursor, neighborCriterion);
         //       System.out.println(nbrIterator.neighborCount());
            }
         //   while(
         //       System.out.println("Cell: "+centralCell.coordinate().toString());
        }

        public void moveNotify(Atom a) {
            ((Coordinate)a.coordinate).assignCell();
        }
        
        public void addMolecule(Molecule m) {
            m.atomIterator.reset();
            while(m.atomIterator.hasNext()) {
                ((Coordinate)m.atomIterator.next().coordinate).setLattice(lattice);
            }
        }
        
        public void deleteMolecule(Molecule m) {
            m.atomIterator.reset();
            while(m.atomIterator.hasNext()) {
                ((Coordinate)m.atomIterator.next().coordinate).setLattice(null);
            }
        }
        
        /**
         * Clears cells of all atoms and reassigns all atoms to their cells
         */
        public void reset() {
            latticeIterator.reset();
            while(latticeIterator.hasNext()) {((AtomCell)latticeIterator.next()).setFirst(null);}
            for(Atom a=phase().firstAtom(); a!=null; a=a.nextAtom()) {
                ((Coordinate)a.coordinate()).setLattice(lattice);
            }
        }
        public final Atom.Iterator makeAtomIteratorUpNeighbor() {return new AtomIteratorUpNeighbor();}
        public final Atom.Iterator makeAtomIteratorDownNeighbor() {return new AtomIteratorDownNeighbor();}
        public final Atom.Iterator makeAtomIteratorUp() {return new AtomIteratorUp();}
        public final Atom.Iterator makeAtomIteratorDown() {return null;} //not yet implemented new AtomIteratorDown();}

        /**
         * Iterator for atom pairs, going up cell-based neighbor list
         * Takes a given atom and generates pairs from it and uplist neighbors
         * Handles (untested as of yet) case where given atom is not in phase
         */
        //member class of Space2DCell
        private final class AtomIteratorUpNeighbor implements Atom.Iterator {
            private Coordinate nextCoordinate;
            private boolean hasNext;
            private SiteIterator.Neighbor cellIterator;
            private Atom referenceAtom;
            public AtomIteratorUpNeighbor() {
                reset();
            }
            public boolean hasNext() {return hasNext;}
            public void reset() {  //resets to first atom in list
                reset(referenceAtom);
            }
            public void reset(Atom a) { //resets iterator to begin with atom first uplist from given one
                referenceAtom = a;
                if(a == null) {hasNext = false; return;}
                nextCoordinate = (Coordinate)a.coordinate();
                cellIterator = nextCoordinate.cell.adjacentIterator(); //cell-neighbor iterator for cell containing this atom
                cellIterator.resetUp();  //next cell returned by iterator is first up the neighbor list
                //need to finish current cell before advancing to next one
                nextCoordinate = nextCoordinate.nextNeighbor; //next atom returned is first up list from given atom
                hasNext = true;
                advanceCell();//if nextCoordinate is non-null, this will do nothing
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                while(nextCoordinate == null) {
                    AtomCell cell = cellIterator.hasNext() ? (AtomCell)cellIterator.next() : null;
                    if(cell == null) {        //no more cells;
                        hasNext = false;
                        return;
                    }
                    nextCoordinate = cell.first;
                }
            }
            public Atom next() {
                Atom nextAtom = (Atom)nextCoordinate.parent();
                nextCoordinate = nextCoordinate.nextNeighbor;      //next atom in cell
                if(nextCoordinate == null) {advanceCell();}
                return nextAtom;
            }
            //we cheat a bit in setting up this allAtoms method
            public void allAtoms(AtomAction action) {
                reset();
                while(hasNext) {action.actionPerformed(next());}
            }
        } //end of AtomIteratorUpNeighbor


        /**
         * Iterator for atom pairs, going up cell-based neighbor list
         * Takes a given atom and generates pairs from it and uplist neighbors
         * Handles (untested as of yet) case where given atom is not in phase
         */
        //member class of Space2DCell
        private final class AtomIteratorDownNeighbor implements Atom.Iterator {
            private Coordinate nextCoordinate;
            private boolean hasNext;
            private SiteIterator.Neighbor cellIterator;
            private Atom referenceAtom;
            public AtomIteratorDownNeighbor() {
                reset();
            }
            public boolean hasNext() {return hasNext;}
            public void reset() {  //resets to first atom in list
                reset(referenceAtom);
            }
            public void reset(Atom a) { //resets iterator to begin with atom first downList from given one
                referenceAtom = a;
                if(a == null) {hasNext = false; return;}
                nextCoordinate = (Coordinate)a.coordinate();
                cellIterator = nextCoordinate.cell.adjacentIterator(); //cell-neighbor iterator for cell containing this atom
                cellIterator.resetDown();  //next cell returned by iterator is first down the neighbor list
                //need to finish current cell before advancing to next one
                nextCoordinate = nextCoordinate.previousNeighbor; //next atom returned is first down list from given atom
                hasNext = true;
                advanceCell();//if nextCoordinate is non-null, this will do nothing
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                while(nextCoordinate == null) {
                    AtomCell cell = cellIterator.hasNext() ? (AtomCell)cellIterator.next() : null;
                    if(cell == null) {        //no more cells;
                        hasNext = false;
                        return;
                    }
                    nextCoordinate = cell.first;
                }
            }
            public Atom next() {
                Atom nextAtom = (Atom)nextCoordinate.parent();
                nextCoordinate = nextCoordinate.previousNeighbor;      //next atom in cell
                if(nextCoordinate == null) {advanceCell();}
                return nextAtom;
            }
            //we cheat a bit in setting up this allAtoms method
            public void allAtoms(AtomAction action) {
                reset();
                while(hasNext) {action.actionPerformed(next());}
            }
        } //end of AtomIteratorDownNeighbor


        //member class of Space2DCell
        private final class AtomIteratorUp implements Atom.Iterator {
            private Coordinate nextCoordinate;
            private boolean hasNext;
            private SiteIterator.List.Cursor cellIterator;
            public AtomIteratorUp() {
                cellIterator = latticeIterator.makeCursor(); //latticeIterator is field in CellListIteratorFactory outer class
                reset();
            }
            public boolean hasNext() {return hasNext;}
            private void allDone() {hasNext = false;}
            public void reset() {  //resets to first atom in list
                hasNext = true;
                cellIterator.reset(); 
                advanceCell();
            }
            public void reset(Atom a) { //resets to the point where the given atom is next
                if(a == null) {allDone(); return;}
                nextCoordinate = (Coordinate)a.coordinate();
                //set cellIterator as if it just returned the cell containing this atom
                hasNext = false;  //if cell is not found, atom is not in phase, and want hasNext false
                cellIterator.reset();   
                while(cellIterator.hasNext()) {
                    if(cellIterator.next() == nextCoordinate.cell) {hasNext = true; break;}  //found the cell
                }
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                while(nextCoordinate == null) {
                    AtomCell cell = cellIterator.hasNext() ? (AtomCell)cellIterator.next() : null;
                    if(cell == null) {        //no more cells;
                        allDone();
                        return;
                    }
                    nextCoordinate = cell.first;
                }
            }
            public Atom next() {
                Atom nextAtom = (Atom)nextCoordinate.parent();
                nextCoordinate = nextCoordinate.nextNeighbor;      //next atom in cell
                if(nextCoordinate == null) {advanceCell();}
                return nextAtom;
            }
            public void allAtoms(AtomAction action) {
                cellIterator.reset();
                while(cellIterator.hasNext()) {
                    AtomCell cell = (AtomCell)cellIterator.next();
                    nextCoordinate = cell.first;
                    while(nextCoordinate != null) {
                        action.actionPerformed((Atom)nextCoordinate.parent());
                        nextCoordinate = nextCoordinate.nextNeighbor;
                    }
                }
            }
                    
        } //end of AtomIteratorUp
        
        //member class of Space2DCell
        
        //under development; presently this code is identical to AtomIteratorUp
        private final class AtomIteratorDown implements Atom.Iterator {
            private Coordinate nextCoordinate;
            private boolean hasNext;
            private SiteIterator.List.Cursor cellIterator;
            public AtomIteratorDown() {
                cellIterator = latticeIterator.makeCursor(); //latticeIterator is field in CellListIteratorFactory outer class
                reset();
            }
            public boolean hasNext() {return hasNext;}
            private void allDone() {hasNext = false;}
            public void reset() {  //resets to first atom in list
                hasNext = true;
                cellIterator.reset(); 
                advanceCell();
            }
            public void reset(Atom a) { //resets to the point where the given atom is next
                if(a == null) {allDone(); return;}
                nextCoordinate = (Coordinate)a.coordinate();
                //set cellIterator as if it just returned the cell containing this atom
                hasNext = false;  //if cell is not found, atom is not in phase, and want hasNext false
                cellIterator.reset();   
                while(cellIterator.hasNext()) {
                    if(cellIterator.next() == nextCoordinate.cell) {hasNext = true; break;}  //found the cell
                }
            }
            // Finds first atom of next occupied cell
            private void advanceCell() {
                while(nextCoordinate == null) {
                    AtomCell cell = cellIterator.hasNext() ? (AtomCell)cellIterator.next() : null;
                    if(cell == null) {        //no more cells;
                        allDone();
                        return;
                    }
                    nextCoordinate = cell.first;
                }
            }
            public Atom next() {
                Atom nextAtom = (Atom)nextCoordinate.parent();
                nextCoordinate = nextCoordinate.nextNeighbor;      //next atom in cell
                if(nextCoordinate == null) {advanceCell();}
                return nextAtom;
            }
            public void allAtoms(AtomAction action) {
                cellIterator.reset();
                while(cellIterator.hasNext()) {
                    AtomCell cell = (AtomCell)cellIterator.next();
                    nextCoordinate = cell.first;
                    while(nextCoordinate != null) {
                        action.actionPerformed((Atom)nextCoordinate.parent());
                        nextCoordinate = nextCoordinate.nextNeighbor;
                    }
                }
            }
        } //end of AtomIteratorDown
    }//end of CellListIteratorFactory
    
    
    // might be better to define this potential's affected atoms iterator to return no atoms?
    public static class CellPotential extends PotentialField implements PotentialField.Hard {
        
   //     Phase phase;
        Space2D.Vector d;  //dimensions of phase
        private static final double eps = 1.e-8;
        
        public CellPotential(Phase p) {
            super(p);
            //need dimensions of phase boundary to normalize atom coordinates to unity
            setPhaseBoundary(p.boundary());
            p.boundaryMonitor.addObserver(boundaryObserver());
        }
        
	    private Observer boundaryObserver() {
	        return new Observer() {
	            //This is the action that is to happen if phase takes a new boundary
	            public void update(Observable o, Object arg) {
	                setPhaseBoundary((Space.Boundary)arg);
	            }
	        };
	    }
	    void setPhaseBoundary(Space.Boundary b) {
	        if(b == null) return;
	        d = (Space2D.Vector)b.dimensions();
	    }
	    
	    public double energy(Atom a) {return 0.0;}
	    
        public void bump(Atom atom) {
            Coordinate atomCoordinate = (Coordinate)atom.coordinate();
            AtomCell cell = atomCoordinate.cell;
            //System.out.println(atomCoordinate.r.x+"  "+atomCoordinate.r.y);
            double dx = atomCoordinate.r.x - d.x*cell.position.x;  //position is a field in AtomCell
            double dy = atomCoordinate.r.y - d.y*cell.position.y;
            //System.out.println("dx, dy: "+dx + "  " + dy);
            //System.out.println("Old cell "+cell.coordinate().toString());
            if(Math.abs(dx) > Math.abs(dy)) { //collision with left or right wall
              //  System.out.println("x "+atomCoordinate.r.x+"  "+atomCoordinate.p.x);
                if(dx > 0) {atomCoordinate.r.x = d.x*cell.E.xW/*+eps*/; atomCoordinate.assignCell(cell.E);}
                else       {atomCoordinate.r.x = d.x*cell.W.xE/*-eps*/; atomCoordinate.assignCell(cell.W);}
              //  System.out.println("x "+atomCoordinate.r.x+"  "+atomCoordinate.p.x);
            }
            else {
              //  System.out.println("y "+atomCoordinate.r.y+"  "+atomCoordinate.p.y);
                if(dy > 0) {atomCoordinate.r.y = d.y*cell.S.yN/*+eps*/; atomCoordinate.assignCell(cell.S);}
                else       {atomCoordinate.r.y = d.y*cell.N.yS/*-eps*/; atomCoordinate.assignCell(cell.N);}
              //  System.out.println("y "+atomCoordinate.r.y+"  "+atomCoordinate.p.y);
            }
            //System.out.println("New cell "+((Coordinate)atom.coordinate()).cell.coordinate().toString());
            //System.out.println();
        }
        
        public double collisionTime(Atom atom) {
            Coordinate atomCoordinate = (Coordinate)atom.coordinate;
            AtomCell cell = atomCoordinate.cell;
            double dx = (atomCoordinate.p.x > 0) ? d.x*cell.xE - atomCoordinate.r.x : d.x*cell.xW - atomCoordinate.r.x; 
            double tx = Math.abs(dx/atomCoordinate.p.x);
            double dy = (atomCoordinate.p.y > 0) ? d.y*cell.yS - atomCoordinate.r.y : d.y*cell.yN - atomCoordinate.r.y;
            double ty = Math.abs(dy/atomCoordinate.p.y);  //remember, down (South) is positive
            return (tx > ty) ? ty*atomCoordinate.parent().mass() : tx*atomCoordinate.parent().mass();
//        double tnew = Math.abs((0.5*atom.phase().parentSimulation.space.getNeighborDistance()-0.5*1.0001*((AtomType.Disk)atom.type).diameter())/Math.sqrt(atom.coordinate.momentum().squared()));  //assumes range of potential is .le. diameter
        }
        
  //      public double lastCollisionVirial() {return 0.0;}
  //      public Space.Tensor lastCollisionVirialTensor() {return Space2D.Tensor.ZERO;}
  //      public double energy(AtomPair pair) {return 0.0;}
  //      public boolean overlap(AtomPair pair) {return false;}
  //      public double energyLRC(int i, int j, double x) {return 0.0;}
    }//end of Space2DCell.CellPotential
    
    protected static class BoundaryPeriodicSquare extends Space2D.BoundaryPeriodicSquare {
        public BoundaryPeriodicSquare() {super();}
        public BoundaryPeriodicSquare(Phase p) {super(p);}
        public void centralImage(Coordinate c) {}  //central image enforcement is done with cells
        public void centralImage(Space.Vector r) {}
        public void centralImage(Vector r) {}
        public double collisionTime(AtomPair pair) {return Double.MAX_VALUE;}
        public void bump() {}
    }//end of BoundaryPeriodicSquare
    
    /**
     * A factory that makes Sites of type AtomCell
     */
    private static class AtomCellFactory implements SiteFactory {
        public Site makeSite(AbstractLattice parent, SiteIterator.Neighbor iterator, AbstractLattice.Coordinate coord) {
            if(!(parent instanceof SquareLattice)) {
                throw new IllegalArgumentException("Space2DCell.AtomCellFactory: parent lattice must be of type SquareLattice");
            }
            if(!(coord instanceof BravaisLattice.Coordinate)) {
                throw new IllegalArgumentException("Space2DCell.AtomCell.Factory: coordinate must be of type BravaisLattice.Coordinate");
            }
            return (new AtomCell((SquareLattice)parent, (BravaisLattice.Coordinate)coord));
//            return ((SquareLattice)parent).new Cell(coord);
        }
    }//end of SquareLattice.CellFactory
    
    /**
     * A lattice cell that holds a reference to an atom coordinate, which marks the beginning of list of atoms in that cell.
     */
    private static class AtomCell extends SquareLattice.Cell {
        public Coordinate first;
        public Space2D.Vector position;
        public Color color;
        public AtomCell N, E, S, W;   //immediately adjacent cells (North, East, South, West)
        public double yN, xE, yS, xW; //x, y coordinates of boundaries of cell
        public AtomCell(SquareLattice parent, BravaisLattice.Coordinate coord) {
            super(parent, coord);
            color = Constants.RandomColor();
//            position = (Space2D.Vector)coord.position();
        }
        public Coordinate first() {return first;}
        public void setFirst(Coordinate f) {first = f;}
    }//end of AtomCell
    
    public static class ColorByCell extends ColorScheme {
        public void colorAtom(Atom a) {
            a.setColor(((Coordinate)a.coordinate()).cell.color);
        }
    }
    
    
    /**
     * Method for testing and demonstrating use of class
     */
    public static void main(String[] args) {
        java.awt.Frame f = new java.awt.Frame();   //create a window
        f.setSize(600,350);

        Simulation.instance = new Simulation(new Space2DCell());
        
	    IntegratorHard integratorHard1 = new IntegratorHard();
	    SpeciesDisks speciesDisks1 = new SpeciesDisks();
	    speciesDisks1.setAtomsPerMolecule(1);
	    speciesDisks1.setNMolecules(20);
	    Phase phase1 = new Phase();
	    P2HardDisk P2HardDisk1 = new P2HardDisk();  //P2 must follow species until setSpeciesCount is fixed
	    P1TetherHardDisk P1TetherHardDisk1 = new P1TetherHardDisk();
	    Controller controller1 = new Controller();
	    controller1.setMakeButton(true);
	    DisplayPhase displayPhase1 = new DisplayPhase();
	    IntegratorMD.Timer timer = integratorHard1.new Timer(integratorHard1.chronoMeter());
	    timer.setUpdateInterval(10);
		Simulation.instance.setBackground(Color.yellow);
		
		displayPhase1.setColorScheme(new Space2DCell.ColorByCell());
		integratorHard1.setIsothermal(true);
		integratorHard1.setTemperature(50.0);
		
		((CellListIteratorFactory)phase1.iteratorFactory()).setNeighborDistance(1.1*Default.ATOM_SIZE);
	//	((CellListIteratorFactory)phase1.iteratorFactory()).setNCells(10,10);
                                            
		Simulation.instance.elementCoordinator.go(); //invoke this method only after all elements are in place
		                                    //calling it a second time has no effect
		                                    
        f.add(Simulation.instance);         //access the static instance of the simulation to
                                            //display the graphical components
        f.pack();
        f.show();
        f.addWindowListener(new java.awt.event.WindowAdapter() {   //anonymous class to handle window closing
            public void windowClosing(java.awt.event.WindowEvent e) {System.exit(0);}
        });
    }
    
}//end of Space2DCell