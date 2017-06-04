/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.atom;

import java.util.HashMap;

import etomica.api.IElement;

/**
 * This class hashes atomic diameters based on the element.  A static method is
 * provided that populates the hash with some default diameters based on VWD
 * radii.
 * 
 * @author Andrew Schultz
 */
public class DiameterHashByElement implements DiameterHash {

    public DiameterHashByElement() {
        elementDiameterHash = new HashMap<String,Double>();
    }
    
    public double getDiameter(IAtom atom) {
        return getDiameter(atom.getType().getElement());
    }
    
    public double getDiameter(IElement element) {
        Double dbl = elementDiameterHash.get(element.getSymbol());
        if (dbl != null) return dbl;
        return -1;
    }

    public void setDiameter(String elementSymbol, double newDiameter) {
        elementDiameterHash.put(elementSymbol, newDiameter);
    }
    
    public void clear() {
        elementDiameterHash.clear();
    }

    protected final HashMap<String,Double> elementDiameterHash;

    /**
     * Populates the given instance with diameters based on atomic radii
     * (taken from various sources).
     */
    public static void populateVDWDiameters(DiameterHashByElement manager) {
        manager.clear();
        // VDW radii taken from
        //  M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer, and D. G. Truhlar,
        //  "Consistent van der Waals Radii for the Whole Main Group",
        //  J. Phys. Chem. A 2009, 113, 5806­5812
        // (who in turn found some new values and reported those in combination with others)
        manager.setDiameter("H", 2*1.1);
        manager.setDiameter("He", 2*1.4);
        manager.setDiameter("Li", 2*1.81);
        manager.setDiameter("B", 2*1.92);
        manager.setDiameter("C", 2*1.7);
        manager.setDiameter("N", 2*1.55);
        manager.setDiameter("O", 2*1.52);
        manager.setDiameter("F", 2*1.47);
        manager.setDiameter("Ne", 2*1.54);  // 10
        manager.setDiameter("Na", 2*2.27);
        manager.setDiameter("Mg", 2*1.73);
        manager.setDiameter("Al", 2*1.84);
        manager.setDiameter("Si", 2*2.1);
        manager.setDiameter("P", 2*1.8);
        manager.setDiameter("S", 2*1.8);
        manager.setDiameter("Cl", 2*1.75);
        manager.setDiameter("Ar", 2*1.88);
        manager.setDiameter("K", 2*2.75);
        manager.setDiameter("Ca", 2*2.31);  // 20
        //manager.setDiameter("Cu", 2*1.4);
        manager.setDiameter("Ga", 2*1.87);  // 31
        manager.setDiameter("Ge", 2*2.11);
        manager.setDiameter("As", 2*1.85);
        manager.setDiameter("Se", 2*1.9);
        manager.setDiameter("Br", 2*1.83);
        manager.setDiameter("Kr", 2*2.02);
        manager.setDiameter("Rb", 2*3.03);
        manager.setDiameter("Sr", 2*2.49);  //38
        manager.setDiameter("In", 2*1.93);
        //manager.setDiameter("Sn", 2*2.17); // 50
        manager.setDiameter("Sb", 2*2.06);
        manager.setDiameter("Te", 2*2.06);
        manager.setDiameter("I", 2*1.98);
        manager.setDiameter("Xe", 2*2.16);
        manager.setDiameter("Cs", 2*3.43);
        manager.setDiameter("Ba", 2*2.68);
        manager.setDiameter("Tl", 2*1.96);
        manager.setDiameter("Pb", 2*2.02);
        manager.setDiameter("Bi", 2*2.07);
        manager.setDiameter("Po", 2*1.97);
        manager.setDiameter("At", 2*2.02);
        manager.setDiameter("Rn", 2*2.20);
        manager.setDiameter("Fr", 2*3.48);
        manager.setDiameter("Ra", 2*2.83);

        // msellers' diameters
        // The following values come from either the ASM Handbook or Cullity & Stock's 
        // "Elements of X-Ray Diffraction" (2001)
        
        manager.setDiameter("Cu", 2.5561);
        manager.setDiameter("Ag", 2.8895);
        manager.setDiameter("Sn", 3.022);
    }
}
