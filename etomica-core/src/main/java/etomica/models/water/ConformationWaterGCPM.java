/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.water;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

public class ConformationWaterGCPM implements IConformation, java.io.Serializable {

    private static final long serialVersionUID = 1L;
    private double bondLengthOH = 0.9572;
    private double angleHOH = 104.52*Math.PI/180.;
    private double rOM=0.27;

    public ConformationWaterGCPM(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
        
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        
        
        IAtom o = list.get(SpeciesWater4P.indexO);
        o.getPosition().E(new double[] {x, y, z});
               
        IAtom h1 = list.get(SpeciesWater4P.indexH1);
        h1.getPosition().E(new double[] {x-bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y, z+bondLengthOH*Math.sin((Math.PI-angleHOH)/2)});
                
        IAtom h2 = list.get(SpeciesWater4P.indexH2);
        h2.getPosition().E(new double[] {x+bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y, z+bondLengthOH*Math.sin((Math.PI-angleHOH)/2)});
        
        IAtom m = list.get(SpeciesWater4P.indexM);
        m.getPosition().E(new double[] {x, y, z+rOM});

        
/*        Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x-bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y+bondLengthOH*Math.sin((Math.PI-angleHOH)/2), 0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos((Math.PI-angleHOH)/2), y+bondLengthOH*Math.sin((Math.PI-angleHOH)/2), 0});
        
        Atom m = iterator.nextAtom();
        m.coord.position().E(new double[] {x, y+rOM, 0});
*/
        
 /*       Atom o = iterator.nextAtom();
        o.coord.position().E(new double[] {x, y, 0.0});
               
        Atom h1 = iterator.nextAtom();
        h1.coord.position().E(new double[] {x+bondLengthOH, y, 0.0});
                
        Atom h2 = iterator.nextAtom();
        h2.coord.position().E(new double[] {x+bondLengthOH*Math.cos(angleHOH), y+bondLengthOH*Math.sin(angleHOH), 0.0});
        
        Atom m = iterator.nextAtom();
        m.coord.position().E(new double[] {x+rOM*Math.cos(angleHOH/2.0), y+rOM*Math.sin(angleHOH/2.0), 0.0});
*/
    }//end of initializePositions

    protected final Space space;
}
