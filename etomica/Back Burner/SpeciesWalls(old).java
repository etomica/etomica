package simulate;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class SpeciesWalls extends Species {

    int borderTol;
    boolean boundary, vertical, horizontal;
    double diameter;
    int thickness;
    Color color;

    public SpeciesWalls() {
        super(4);
        setSpeciesIndex(1);
        name = "Wall";
        borderTol = 20;
        setThickness(4);
        boundary = vertical = false;
        setHorizontal(true);
        setFillVolume(true);
    }

  public void setDefaults() {
    setDiameter(horizontal ? parentPhase.space.dimensions[0] : parentPhase.space.dimensions[1]);
  }

  //makeMolecules is called by constructor via setNMolecules
  void makeMolecules() {
    molecule = new Molecule[nMolecules];
    for(int i=0; i<nMolecules; i++) {molecule[i] = new MoleculeWall(this);}
  }
  
  //initializeMolecules is called by constructor via setNMolecules
  void initializeMolecules() {
    initializeMolecules(diameter, color);
  }
  void initializeMolecules(double d, Color c) {
    setDiameter(d);    //call set methods to pass diameter and mass to atoms
    setColor(c);
  }

    public final double getDiameter() {return diameter;}
    public final void setDiameter(double d) {
        diameter = d;
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setDiameter(d);}
    }
    
    public final Color getColor() {return color;}
    public final void setColor(Color c) {
        color = c;
        for(Atom a=firstAtom; a!=lastAtom.getNextAtom(); a=a.getNextAtom()) {a.setColor(c);}
    }
        

    public double getMass() {return Double.MAX_VALUE;}

    public boolean isBoundary() {return boundary;}
    public void setBoundary(boolean b) {boundary = b;}

    public boolean isVertical() {return vertical;}
    public void setVertical(boolean b) {
        vertical = b;
        if(vertical) {setHorizontal(false);}
    }

    public boolean isHorizontal() {return horizontal;}
    public void setHorizontal(boolean b) {
        horizontal = b;
        if(horizontal) {setVertical(false);}
    }

    public int getBorderTol() {return borderTol;}
    public void setBorderTol(int borderTol) {this.borderTol = borderTol;}

    public int getThickness() {return thickness;}
    public void setThickness(int thickness) {this.thickness = thickness;}

  public void initializeSpecies(Phase phase) {
    parentPhase = phase;
    int x = getLocation().x;
    int y = getLocation().y;
 //   double scale = phase.getNominalScale()/(2.0*phase.getImageShells()+1);   //temporary workaround to get scale
 //   double scale = Beans.isDesignTime() ? 1.0 : phase.space.getScale(); //needs some work
    double scale = 1.0;
    double toPixels = scale*Phase.TO_PIXELS;
    if(vertical) {
        d[0] = 0.0;
        r[1] = (double)y/toPixels;
        if(fillVolume) {r[1] = 0.0; d[1] = designTimeYDim;}
        else {d[1] = (double)getSize().width/toPixels;}
        if(boundary) {
            if(x < phase.getSize().width/2) {r[0] = 0.0;}
            else {r[0] = designTimeXDim;}
        }
        else {r[0] = (double)x/toPixels;}
        if(Beans.isDesignTime()) setBounds((int)(toPixels*r[0]),(int)(toPixels*r[1]),thickness,(int)(toPixels*d[1]));
    }
    else if(horizontal) {
        d[1] = 0.0;
        r[0] = (double)x/toPixels;
        if(fillVolume) {r[0] = 0.0; d[0] = designTimeXDim;}
        else {d[0] = (double)getSize().height/toPixels;}
        if(boundary) {
            if(y < phase.getSize().height/2) {r[1] = 0.0;}
            else {r[1] = designTimeYDim;}
        }
        else {r[1] = (double)y/toPixels;}
        if(Beans.isDesignTime()) setBounds((int)(toPixels*a.r[0]),(int)(toPixels*a.r[1]),(int)(toPixels*a.d[0]),thickness);
    }
  }

    public void draw(Graphics g, int[] origin, double scale){
      Atom nextSpeciesAtom = lastAtom.getNextAtom();
      for(Atom a=firstAtom; a!=nextSpeciesAtom; a=a.getNextAtom()) {
        a.draw(g, origin, scale);
      }
    }
}
