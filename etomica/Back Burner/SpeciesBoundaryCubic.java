package simulate;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class SpeciesBoundaryCubic extends Species {

    double diameter;
    int thickness;
//    AtomWall firstAtom, lastAtom;
//    MoleculeWall[] molecule;

    public SpeciesBoundaryCubic() {
        super(4);
        setSpeciesIndex(1);
        name = "Boundary";
        nAtomsPerMolecule = 1;
    }

  public void setDefaults() {;}

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
        
    public final int getThickness() {return ((AtomWall)firstAtom).getThickness();}
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
    AtomWall a = (AtomWall)molecule[0].firstAtom;
    a.setHorizontal(true); a.r[0] = 0.0; a.r[1] = 0.0; a.setDiameter(xmax);  //top
    a = (AtomWall)molecule[1].firstAtom;
    a.setVertical(true); a.r[0] = xmax-(double)getThickness()/Phase.TO_PIXELS; a.r[1] = 0.0; a.setDiameter(ymax);  //right
    a = (AtomWall)molecule[2].firstAtom;
    a.setHorizontal(true); a.r[0] = 0.0; a.r[1] = ymax-(double)getThickness()/Phase.TO_PIXELS; a.setDiameter(xmax); //bottom
    a = (AtomWall)molecule[3].firstAtom;
    a.setVertical(true); a.r[0] = 0.0; a.r[1] = 0.0; a.setDiameter(ymax);  //left
  }

    public void draw(Graphics g, int[] origin, double scale){
      Atom nextSpeciesAtom = lastAtom.getNextAtom();
      for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        a.draw(g, origin, scale);
      }
    }
}
