package simulate;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class SpeciesWalls extends Species {

    double diameter;
    int thickness;
    AtomWall firstAtom, lastAtom;
    MoleculeWall[] molecule;

    public SpeciesWalls() {
        super(4);
        setSpeciesIndex(1);
        name = "Wall";
        setThickness(4);
        setFillVolume(true);
    }

  public void setDefaults() {
  }

  //makeMolecules is called by constructor via setNMolecules
  void makeMolecules() {
    molecule = new MoleculeWall[nMolecules];
    for(int i=0; i<nMolecules; i++) {molecule[i] = new MoleculeWall(this);}
  }
  
  //initializeMolecules is called by constructor via setNMolecules
  void initializeMolecules() {
    initializeMolecules(Color.orange);
  }
  void initializeMolecules(Color c) {
    setColor(c);
  }
    
    public final Color getColor() {return firstAtom.getColor();}
    public final void setColor(Color c) {
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setColor(c);}
    }
        
    public final int getThickness() {return firstAtom.getThickness();}
    public final void setThickness(int t) {
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {((AtomWall)a).setThickness(t);}
    }
        
  public void initializeSpecies(Phase phase) {
    parentPhase = phase;
    int x = getLocation().x;
    int y = getLocation().y;
    setBounds(0,0,phase.getSize().width, phase.getSize().height);
    double xmax = phase.space.dimensions[0];
    double ymax = phase.space.dimensions[1];
    AtomWall a = molecule[0].firstAtom;
    a.setHorizontal(true); a.r[0] = 0.0; a.r[1] = 0.0; a.setDiameter(xmax);  //top
    a = molecule[1].firstAtom;
    a.setVertical(true); a.r[0] = xmax-getThickness(); a.r[1] = 0; a.setDiameter(ymax);  //right
    a = molecule[2].firstAtom;
    a.setHorizontal(true); a.r[0] = ymax-getThickness(); a.r[1] = xmax; a.setDiameter(xmax); //bottom
    a = molecule[3].firstAtom;
    a.setVertical(true); a.r[0] = 0.0; a.r[1] = 0.0; a.setDiameter(ymax);  //left
  }

    public void draw(Graphics g, int[] origin, double scale){
      Atom nextSpeciesAtom = lastAtom.getNextAtom();
      for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        a.draw(g, origin, scale);
      }
    }
}
