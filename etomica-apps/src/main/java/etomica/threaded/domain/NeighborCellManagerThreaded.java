/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded.domain;

import etomica.box.Box;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.nbr.cell.NeighborCellManager;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

public class NeighborCellManagerThreaded extends NeighborCellManager {

    public int totalCells;
    
    public NeighborCellManagerThreaded(Simulation sim, Box box, double potentialRange, Space space) {
        super(sim, box, potentialRange, space);
        // TODO Auto-generated constructor stub
    }

    public NeighborCellManagerThreaded(Simulation sim, Box box, double potentialRange,
                                       IMoleculePositionDefinition positionDefinition, Space space) {
        super(sim, box, potentialRange, positionDefinition, space);
        // TODO Auto-generated constructor stub
    }
    
    
    
    public boolean checkDimensions(){
        setNumCells(totalCells);
        return true;
    }
    
    //Divides up simulation box for threads
    public void setNumCells(int totalCells){
        this.totalCells = totalCells;
        if(totalCells==0){return;}
        Vector dimensions = box.getBoundary().getBoxSize();
        lattice.setDimensions(dimensions);
       
        int [] nCells = calculateLatticeDimensions(totalCells, dimensions);
     
        //only update the lattice (expensive) if the number of cells changed
        int[] oldSize = lattice.getSize();
        for (int i=0; i<nCells.length; i++) {
            if (oldSize[i] != nCells[i]) {
                lattice.setSize(nCells);
                break;
            }
        }
        
        
        
    }
    
    protected int[] calculateLatticeDimensions(int nCells, Vector shape) {
        int dimLeft = shape.getD();
        int nCellsLeft = nCells;
        int[] latticeDimensions = new int[shape.getD()];
        while (dimLeft > 0) {
            double smin = Double.POSITIVE_INFINITY;
            int dmin = 0;
            double product = 1.0;
            for (int idim = 0; idim < shape.getD(); idim++) {
                if (latticeDimensions[idim] > 0)
                    continue;
                if (shape.getX(idim) < smin) {
                    smin = shape.getX(idim);
                    dmin = idim;
                }
                product *= shape.getX(idim);
            }
            // round off except for last dimension (then round up)
            if (dimLeft > 1) {
                latticeDimensions[dmin] = (int) Math.round(shape.getX(dmin)
                        * Math.pow((nCellsLeft / product), 1.0 / dimLeft));
                if(latticeDimensions[dmin] == 0){
                    latticeDimensions[dmin] = 1;
                }
            } else {
                latticeDimensions[dmin] = (int) Math.ceil(shape.getX(dmin)
                        * nCellsLeft / product);
            }
            
            nCellsLeft = (nCellsLeft + latticeDimensions[dmin] - 1)
                    / latticeDimensions[dmin];
            dimLeft--;
        }
        return latticeDimensions;
    }

}
