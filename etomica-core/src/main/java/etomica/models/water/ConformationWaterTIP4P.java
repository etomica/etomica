/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;
import etomica.units.Electron;

/**
 * Conformation for 4-point water molecule.
 */
public class ConformationWaterTIP4P implements IConformation, java.io.Serializable {

    public ConformationWaterTIP4P(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
        
        IAtom o = list.get(SpeciesWater4P.indexO);
        o.getPosition().E(new double[] {0, 0, 0.0});

        double x = bondLengthOH*Math.sin(0.5*angleHOH);
        double y = bondLengthOH*Math.cos(0.5*angleHOH);
        
        IAtom h1 = list.get(SpeciesWater4P.indexH1);
        h1.getPosition().E(new double[] {-x, y, 0.0});
                
        IAtom h2 = list.get(SpeciesWater4P.indexH2);
        h2.getPosition().E(new double[] {+x, y, 0.0});
        
        IAtom m = list.get(SpeciesWater4P.indexM);
        m.getPosition().E(new double[] {0, rOM, 0.0});
    }
    
    public final static double [] Echarge = new double [4];
    static {
        ConformationWaterTIP4P.Echarge[SpeciesWater4P.indexH1] = Electron.UNIT.toSim( 0.52);
        ConformationWaterTIP4P.Echarge[SpeciesWater4P.indexH2] = Electron.UNIT.toSim( 0.52);
        ConformationWaterTIP4P.Echarge[SpeciesWater4P.indexO] = Electron.UNIT.toSim( 0.00);
        ConformationWaterTIP4P.Echarge[SpeciesWater4P.indexM] = Electron.UNIT.toSim(-1.04);
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    public static final double bondLengthOH = 0.9572;
    public static final double angleHOH = 104.52*Math.PI/180.;
    public static final double rOM=0.15;
}
