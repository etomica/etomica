package etomica;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorListCompound;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Fills phase with molecules on a lattice, taking each molecule in successive order
 * from the linked list of molecules.  Takes no special action when list moves from
 * one species to the next.
 * Wall "molecules" are ignored becuase the (super) add method will not add them.
 *
 * Need to improve this to handle different dimensions more elegantly
 */

/* History
 * 01/04/03 (SKK/DAK) added hexagonal lattice for 2D configuration
 * 01/14/03 (DAK) fixed typo in name of getSquareConfig method
 */
 
public class ConfigurationSequential extends Configuration {

	private boolean fill;
	private boolean squareConfig;
    private final AtomIteratorListCompound atomIterator;
    private final AtomActionTranslateTo atomActionTranslateTo;

	public ConfigurationSequential(Space space) {
		super(space);
		setFillVertical(true);
		setSquareConfig(false); // hexagonalLattice is Default!!
        atomIterator = new AtomIteratorListCompound();
        atomActionTranslateTo = new AtomActionTranslateTo(space);
	}
    
	public void setFillVertical(boolean b) {fill = b;}
	public boolean getFillVertical() {return fill;}
    
	public void setSquareConfig(boolean b){ squareConfig = b;}
	public boolean getSquareConfig() {return squareConfig;}
    
    public void initializePositions(AtomList[] lists) {
        atomIterator.setLists(lists);

        double Lx = dimensions[0];
        double Ly = 0.0;
        if(dimensions.length>1)  Ly = dimensions[1];

        int sumOfMolecules = atomIterator.size();
         
        if(sumOfMolecules == 0) return;
 //       System.out.println("ConfigurationSequential sumOfMolecules = "+sumOfMolecules);
        
        Vector[] rLat;
        
        switch(space.D()) {
            case 1:
                rLat = lineLattice(sumOfMolecules, Lx);
                break;
            default:
            case 2:
//skkwak
                if(squareConfig){rLat = squareLattice(sumOfMolecules, Lx, Ly, fill);
                } else {rLat = hexagonalLattice(sumOfMolecules,Lx,Ly,fill);}
         
                break;
            case 3:
                rLat = null;
///                rLat = new etomica.lattice.LatticeFCC(sumOfMolecules, Default.BOX_SIZE).positions();//ConfigurationFcc.lattice(sumOfMolecules);
                break;
        }
        
   // Place molecules     
        int i = 0;
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            //initialize coordinates of child atoms
            if (!a.node.isLeaf()) {
                Conformation config = a.type.creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup)a.node).childList);
                atomActionTranslateTo.setDestination(rLat[i]);
                atomActionTranslateTo.actionPerformed(a);
            }
            else {
                a.coord.position().E(rLat[i]);
            }
            i++;
        }
    }

    public static etomica.space1d.Vector1D[] lineLattice(int n, double Lx) {
        etomica.space1d.Vector1D[] r = new etomica.space1d.Vector1D[n];
        double delta = Lx/n;
        for(int i=0; i<n; i++) {
            r[i] = new etomica.space1d.Vector1D();
            r[i].setX(0, (i+0.5)*delta);
        }
        return r;
    }
        
    /**
     * Returns a set of n coordinates filling a square lattice of sides Lx and Ly
     * If n is not suitable for square lattice, then last sites are left unfilled
     * Lattice is centered between (0,0) and (Lx, Ly).
     * The final argument should be passed one of the class variables VERTICAL or HORIZONTAL, indicating
     *   whether successive points fill the lattice across or down.
     */
    public final static etomica.space2d.Vector2D[] squareLattice(int n, double Lx, double Ly, boolean fillVertical) {
        etomica.space2d.Vector2D[] r = new etomica.space2d.Vector2D[n];
        for(int i=0; i<n; i++) {r[i] = new etomica.space2d.Vector2D();}

        int moleculeColumns, moleculeRows;
        double moleculeInitialSpacingX, moleculeInitialSpacingY;

    //  Number of molecules per row (moleculeColumns) and number of rows (moleculeRows)
    //  in initial configuration
        moleculeColumns = (int)Math.sqrt(Lx/Ly*n);
        moleculeRows = n/moleculeColumns;

        if(moleculeRows*moleculeColumns < n) moleculeRows++;

    //moleculeColumns may be greater than the actual number of columns drawn
    //Need to center columns in the initial position.
        int columnsDrawn = (int)((double)n/(double)moleculeRows - 1.0E-10) + 1;
        
    //moleculeColumnsShift used to center initial coordinates
        double moleculeColumnsShift = Lx/columnsDrawn/2;
        double moleculeRowsShift = Ly/moleculeRows/2;
  
    //assign distance between molecule centers
        moleculeInitialSpacingX = Lx/columnsDrawn;
        moleculeInitialSpacingY = Ly/moleculeRows;
        int i = 0;
        int ix = 0;
        int iy = 0;
        while(i < n) {
            r[i].setX(0, ix*moleculeInitialSpacingX + moleculeColumnsShift);
            r[i].setX(1, iy*moleculeInitialSpacingY + moleculeRowsShift);
            i++;
            if(fillVertical) {
                iy++;
                if(iy >= moleculeRows) {
                    iy = 0;
                    ix++;
                }}
            else {
                ix++;
                if(ix >= columnsDrawn) {
                    ix = 0;
                    iy++;
                }
            }
        }
        return r;
    }//end of squareLattice
 
    public final static etomica.space2d.Vector2D[] hexagonalLattice(int n, double Lx, double Ly, boolean fillVertical) {
        etomica.space2d.Vector2D[] r = new etomica.space2d.Vector2D[n];
        if(n == 0) return r;
                        
        for(int i=0; i<n; i++) {r[i] = new etomica.space2d.Vector2D();}

        int[] moleculeL = new int[2];
        double[] L = new double[2];
        L[0] = Lx;
        L[1] = Ly;
 
        int fillDim = fillVertical ? 1 : 0;
        // ratio of moleculeRows/moleculeColumns (or vica versa)
        double ratio = L[fillDim]/L[1-fillDim]*Math.sqrt(3.0)*0.5;
        // round off lower of the two and round up the other, ensuring both are even
        if (ratio>1.0) {
            moleculeL[1-fillDim] = 2*Math.max(1,(int)Math.round(Math.sqrt(n/ratio)/2.0));
            moleculeL[fillDim] = ((n + 2*moleculeL[1-fillDim] - 1) / (2*moleculeL[1-fillDim])) * 2;
        }
        else {
            moleculeL[fillDim] = 2*Math.max(1,(int)Math.round(Math.sqrt(n*ratio)/2.0));
            moleculeL[1-fillDim] = ((n + 2*moleculeL[fillDim]- 1) / (2*moleculeL[fillDim])) * 2;
        }
        // pixel spacing between molecules
        double[] moleculeInitialSpacing = new double[2];
        moleculeInitialSpacing[0] = L[0] / moleculeL[0];
        moleculeInitialSpacing[1] = L[1] / moleculeL[1];
        double half = 0.0;
        // counters for each direction
        int[] j = new int[2];
        for (int i=0; i<n; i++) {
            r[i].setX(1-fillDim, (0.5 + j[1-fillDim]) * moleculeInitialSpacing[1-fillDim]);
            r[i].setX(fillDim,   (0.25 + j[fillDim] + half) * moleculeInitialSpacing[fillDim]);
            j[fillDim]++;
            if (j[fillDim] == moleculeL[fillDim]) {
                j[fillDim] = 0;
                if (half == 0.0) {
                    half = 0.5;
                }
                else {
                    half = 0.0;
                }
                j[1-fillDim]++;
            }
        }
           
        return r;
    }   
    
}
