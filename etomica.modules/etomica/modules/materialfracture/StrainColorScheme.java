package etomica.modules.materialfracture;
import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.graphics.ColorScheme;

public class StrainColorScheme extends ColorScheme {
    protected int atomNumber;
    protected IBox box;
    protected boolean hex;
    protected int atomNumberMod;

    public StrainColorScheme() {
        super();
        setHexagonal(true);
    }
    public java.awt.Color centerAtomColor = java.awt.Color.red;
    public java.awt.Color parentColor = java.awt.Color.black;
    public void setBox(IBox newBox){
        box = newBox;
    }
    public IBox getBox(){ return box;}
    
    public void setNumber(double n){
        atomNumber = (int)n;
        atomNumberMod = atomNumber%2;
    }
    public int getNumber(){ return atomNumber;}       


    public void setHexagonal(boolean b){
        hex = b;
    }

    public Color getAtomColor(IAtom a) {
        int idx = a.getParentGroup().getIndex();
        if(hex){
            if( (idx>=atomNumber && idx<atomNumber+18 && idx%2 == atomNumber%2)
                    || (idx>197-atomNumber-18 && idx<=197-atomNumber && idx%2==(197-atomNumber-18)%2) ) return centerAtomColor;
        }
        else {
            if( (idx>=atomNumber && idx<atomNumber+10)
                || (idx>=190-atomNumber && idx<200-atomNumber) ) return centerAtomColor;
        }
        return parentColor;
    }
}