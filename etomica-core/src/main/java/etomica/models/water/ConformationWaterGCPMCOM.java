/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.units.Electron;

/**
 * Conformation for 4-point water molecule (plus 1 site at the COM).
 */
public class ConformationWaterGCPMCOM extends ConformationWaterGCPM {

    public ConformationWaterGCPMCOM(Space space) {
        super(space);
    }
    
    public void initializePositions(IAtomList list){
        super.initializePositions(list);
        
        IAtom o = list.getAtom(SpeciesWater4PCOM.indexO);
               
        IAtom h1 = list.getAtom(SpeciesWater4PCOM.indexH1);
                
        IAtom h2 = list.getAtom(SpeciesWater4PCOM.indexH2);
        
        Vector cr = list.getAtom(SpeciesWater4PCOM.indexC).getPosition();
        cr.Ea1Tv1(Oxygen.INSTANCE.getMass(), o.getPosition());
        cr.PEa1Tv1(Hydrogen.INSTANCE.getMass(), h1.getPosition());
        cr.PEa1Tv1(Hydrogen.INSTANCE.getMass(), h2.getPosition());
        cr.TE(1.0/(Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
    }
    
    public final static double [] Echarge = new double [5];
    static {
        ConformationWaterGCPMCOM.Echarge[SpeciesWater4PCOM.indexH1] = Electron.UNIT.toSim( 0.52);
        ConformationWaterGCPMCOM.Echarge[SpeciesWater4PCOM.indexH2] = Electron.UNIT.toSim( 0.52);
        ConformationWaterGCPMCOM.Echarge[SpeciesWater4PCOM.indexO] = Electron.UNIT.toSim( 0.00);
        ConformationWaterGCPMCOM.Echarge[SpeciesWater4PCOM.indexM] = Electron.UNIT.toSim(-1.04);
        ConformationWaterGCPMCOM.Echarge[SpeciesWater4PCOM.indexC] = Electron.UNIT.toSim( 0.00);
    }
    
}
