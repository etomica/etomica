package etomica.gui;

import etomica.*;
import java.awt.*;

public class P2DefinedPotential extends Potential2 {

    private Class[][] atomPairClassArray;
    public Potential[][] atomPairPotArray;
    
    public P2DefinedPotential(Class[][] array) {
        this(Simulation.instance, array);
    }
    public P2DefinedPotential(Simulation sim, Class[][] array) {
        super(sim);
        atomPairClassArray = (Class[][])array.clone();
        atomPairPotArray = new Potential[atomPairClassArray.length][atomPairClassArray[0].length];
        for (int i = 0; i < atomPairClassArray.length; i++){
            for (int j = 0; j < atomPairClassArray[i].length; j++){
                try { 
                    atomPairPotArray[i][j] = (Potential)atomPairClassArray[i][j].newInstance();
                    String name = ((String)atomPairPotArray[i][j].toString()).substring(9,((String)atomPairPotArray[i][j].toString()).indexOf("@"));
                }
                catch(NullPointerException exc) {}
	            catch(IllegalAccessException exc) {}
	            catch(InstantiationException exc) {}
            }
        }
    }
    
    public final Potential getPotential(Atom a1, Atom a2) {
        return atomPairPotArray[a1.atomIndex()][a2.atomIndex()];
    }
    public final Potential getOnlyPotential() {return atomPairPotArray[0][0];}

}


