package simulate;

import simulate.*;
import java.beans.*;

public class MeterEnergy extends simulate.Meter
{
    public MeterEnergy()
    {
        super();
        setLabel("Energy");
    }

    public double currentValue()
    {
        double ke = 0.0;
        double pe = 0.0;
        for(AtomC a=(AtomC)phase.firstAtom(); a!=null; a=(AtomC)a.getNextAtom()) {
            ke += a.kineticEnergy();
            Atom nextMoleculeAtom = a.getMolecule().lastAtom.getNextAtom();  //first atom on next molecule
            
            Potential1 p1 = phase.potential1[a.getSpeciesIndex()];            
            for(Atom b=a.getNextAtom(); b!=nextMoleculeAtom; b=b.getNextAtom()) { //intramolecular
                pe += p1.getPotential(a,b).energy(a,b);
            }
            
            Potential2[] p2 = phase.potential2[a.getSpeciesIndex()];
            for(Atom b=nextMoleculeAtom; b!=null; b=b.getNextAtom()) {       //intermolecular
                pe += p2[b.getSpeciesIndex()].getPotential(a,b).energy(a,b);
            }
        }
        
        return ke + pe;
    }
}