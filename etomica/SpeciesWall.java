package simulate;
import java.awt.*;//for Graphics
import java.util.*; //for neighborList
import java.beans.Beans;

public class SpeciesWall extends Species {

    int borderTol;
    boolean boundary, vertical, horizontal;

    public SpeciesWall() {
        super(1);
        setSpeciesIndex(1);
        name = "Wall";
        borderTol = 20;
        setThickness(4);
        boundary = vertical = false;
        setHorizontal(true);
        firstElement = lastElement = this;
        neighbors = new Vector();
        setFillVolume(true);
    }

  //makeMolecules is called by constructor via setNMolecules
  void makeMolecules() {
    molecule = new Molecule[nMolecules];
    for(int i=0; i<nMolecules; i++) {molecule[i] = new MoleculeWalls(1);}
    linkMolecules();
    setDiameter(diameter);    //call set methods to pass diameter and mass to atoms
    setMass(mass);
    setColor(color);
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

  public void initializeSpecies(double scale) {
    r[0] = getLocation().x/Phase.TO_PIXELS;
    r[1] = getLocation().y/Phase.TO_PIXELS;
  }

  public void initializeSpecies(Phase phase) {
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
        if(Beans.isDesignTime()) setBounds((int)(toPixels*r[0]),(int)(toPixels*r[1]),(int)(toPixels*d[0]),thickness);
    }
  }

    public void draw(Graphics g, int[] origin, double scale){
        double toPixels = scale*Phase.TO_PIXELS;
        g.setColor(Color.orange);
        int xP = origin[0] + (int)(toPixels*r[0]);
        int yP = origin[1] + (int)(toPixels*r[1]);
        int wP = vertical ? thickness : (int)(toPixels*d[0]);
        int hP = horizontal ? thickness : (int)(toPixels*d[1]);
        g.fillRect(xP,yP,wP,hP);
    }

// SpeciesElement methods

  public double getCollisionTime() {return collisionTime;}
  public void setCollisionTime(double collisionTime) {this.collisionTime = collisionTime;}
  public void decrementCollisionTime(double interval) {this.collisionTime -= interval;}

  public SpeciesElement getCollisionPartner() {return collisionPartner;}
  public void setCollisionPartner(SpeciesElement partner) {this.collisionPartner = partner;}

  public SpeciesElement getNext() {return next;}
  public void setNext(SpeciesElement element) {this.next = element;}

  public SpeciesElement getPrevious() {return previous;}
  public void setPrevious(SpeciesElement element) {this.previous = element;}

//  public void setSpeciesIndex(int index) {this.speciesIndex = index;}
  // getSpeciesIndex is defined in the Species superclass

  public void setRm(double rm) {this.rm = rm;}

  public void zeroForce() {
    Space.uEa1(f,0.0);
  }
  public void addForce(double[] force) {
    Space.uPEv1(f,force);
  }
  public void subtractForce(double[] force) {
    Space.uMEv1(f,force);
  }

  public double getKineticEnergy() {return 0.0;}

  public void addNeighbor(SpeciesElement e) {neighbors.addElement(e);}
  public void clearNeighborList() {neighbors.removeAllElements();}
  public Enumeration getNeighborList() {return neighbors.elements();}

}
