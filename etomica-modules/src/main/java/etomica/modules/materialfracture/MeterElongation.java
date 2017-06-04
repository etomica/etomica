/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.materialfracture;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.units.Null;


public class MeterElongation extends DataSourceScalar {
    
    public MeterElongation() {
        super("Elongation", Null.DIMENSION);
        setHexagonal(true);
    }


    public void setBox(Box newBox) {
        box = newBox;
        originalGageLength = calcGageLength();
    }

    public Box getBox() {
        return box;
    }

    public void setAtomNumber(int n){
        aNumber = n;
        if (box != null) {
            originalGageLength = calcGageLength();
        }
    }

    public void setHexagonal(boolean b){
        hex = b;
        if (box != null) {
            originalGageLength = calcGageLength();
        }
    }
    
    protected double calcGageLength() {
        IAtomList leafList = box.getLeafList();
        if(hex){
            double sum = 0;
            for (int i=aNumber; i<aNumber+18; i+=2) {
                sum -= leafList.getAtom(i).getPosition().getX(0);
            }
            for (int i=197-aNumber; i>197-aNumber-18; i-=2) {
                sum += leafList.getAtom(i).getPosition().getX(0);
            }
            return sum/9.0;
        }
        double firstLine = 0, secondLine = 0;
        for(int i=0; i<10; i++){
//                firstLine =firstLine+phase.speciesMaster.atomList.get((int)(190-aNumber+i)).coord.position().get(0);
//                secondLine=secondLine+phase.speciesMaster.atomList.get(aNumber+i).coord.position().get(0);
        }
        return (firstLine-secondLine)/10.0;
    }
    
    public double getDataAsScalar(){
        System.out.println("l "+calcGageLength()+" "+originalGageLength);
        return calcGageLength()-originalGageLength; 
    }

    protected int aNumber;
    protected double originalGageLength;
    protected boolean hex;
    protected Box box;
}
