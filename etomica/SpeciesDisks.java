package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesDisks extends Species {

///  double mass;
  double diameter;
  double radius;
  
  public SpeciesDisks() {
    super(new AtomHardDisk(null));
    this.add(new ConfigurationMoleculeLinear());
    initializeMolecules(0.1,1.0,Color.black);
  }
  
  //initializeMolecules is called by constructor via setNMolecules
  void initializeMolecules() {
    initializeMolecules(diameter, mass, colorScheme.getBaseColor());
  }
  void initializeMolecules(double d, double m, Color c) {
    setDiameter(d);    //call set methods to pass diameter and mass to atoms
    setMass(m);
    setColor(c);
  }
  
  // Exposed Properties

    public final double getMass() {return mass;}
    public final void setMass(double mass) {
        this.mass = mass;
        if(firstAtom() == null) {return;}  //return if atoms have not yet been ordered
        for(AtomC a=(AtomC)firstAtom(); a!=lastAtom().getNextAtom(); a=(AtomC)a.getNextAtom()) {a.setMass(mass);}
        for(Molecule m=firstMolecule; m!=lastMolecule.getNextMolecule(); m=m.getNextMolecule()) {m.updateMass();}        
    }
    
    public final double getDiameter() {return diameter;}
    public void setDiameter(double d) {
        diameter = d;
        radius = 0.5*d;
        if(firstAtom() == null) {return;}
        for(Atom a=firstAtom(); a!=lastAtom().getNextAtom(); a=a.getNextAtom()) {((AtomDisk)a).setDiameter(d);}
    }
        
    public final Color getColor() {return colorScheme.getBaseColor();}
    public final void setColor(Color c) {
        colorScheme.setBaseColor(c);
        if(firstAtom() == null) {return;}
        for(Atom a=firstAtom(); a!=lastAtom().getNextAtom(); a=a.getNextAtom()) {a.setColor(c);}
    }
}


