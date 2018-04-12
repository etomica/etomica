/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;
import java.util.HashMap;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.box.BoxAgentManager;
import etomica.lattice.FiniteLattice;
import etomica.nbr.PotentialMasterNbr;
import etomica.nbr.cell.NeighborCellManager;

public class ColorSchemeCell extends ColorScheme {
    
    public ColorSchemeCell(PotentialMasterNbr potentialMaster, IRandom random, Box box) {
    	super();
        cellManager = (NeighborCellManager)potentialMaster.getBoxCellManager(box);
        this.random = random;
    }
    
    public void setLattice(FiniteLattice lattice) {
        Object[] sites = lattice.sites();
        for(int i=0; i<sites.length; i++) {
            hash.put(sites[i], new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble()));
        }
    }
    
    public Color getAtomColor(IAtom a) {
        return hash.get(cellManager.getCell(a));
    }
    
    private static final long serialVersionUID = 1L;
    private final HashMap<Object,Color> hash = new HashMap<Object,Color>();
    private final NeighborCellManager cellManager;
    private final IRandom random;
}
