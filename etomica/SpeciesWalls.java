package simulate;
import java.awt.*;//for Graphics
import java.beans.Beans;

public class SpeciesWalls extends Species {

    int borderTol;
    boolean boundary;

    public SpeciesWalls() {
        this(1,1);
    }
    
    public SpeciesWalls(int n, int na) {
        super(n,na);
        initializeMolecules();
    }

  public void initializeMolecules() {
    setThickness(4);
    setAngle(0);
    setColor(Constants.DARK_RED);
  }

  public Molecule makeMolecule() {
    return new MoleculeWall(this,nAtomsPerMolecule);
  }

    public final boolean isBoundary() {return boundary;}
    public final void setBoundary(boolean b) {boundary = b;}

    public final int getBorderTol() {return borderTol;}
    public final void setBorderTol(int borderTol) {this.borderTol = borderTol;}

    public final int getAngle() {return ((AtomWall)firstAtom()).getAngle();}
    public final void setAngle(int theta) {
        ((AtomWall)firstAtom()).setAngle(theta);
    }

    public final int getThickness() {return ((AtomWall)firstAtom()).getThickness();}
    public final void setThickness(int t) {((AtomWall)firstAtom()).setThickness(t);}
    
    public final Color getColor() {return firstAtom().getColor();}
    public final void setColor(Color c) {firstAtom().setColor(c);}

  public void initializeSpecies(Phase phase) {
    parentPhase = phase;
    double toPixels = Phase.TO_PIXELS;
    int x = getLocation().x;
    int y = getLocation().y;
    if(((AtomWall)firstAtom()).isVertical()) {
        firstAtom().r[1] = (double)y/Phase.TO_PIXELS;
        if(fillVolume) {firstAtom().r[1] = 0.0; firstAtom().setDiameter(designTimeYDim);}
        else {firstAtom().setDiameter((double)getSize().width/Phase.TO_PIXELS);}
        if(boundary) {
            if(x < phase.getSize().width/2) {firstAtom().r[0] = 0.0;}
            else {firstAtom().r[0] = designTimeXDim;}
        }
        else {firstAtom().r[0] = (double)x/toPixels;}
        if(Beans.isDesignTime()) setBounds((int)(Phase.TO_PIXELS*firstAtom().r[0]),(int)(Phase.TO_PIXELS*firstAtom().r[1]),getThickness(),(int)(toPixels*firstAtom().getDiameter()));
    }
    else if(((AtomWall)firstAtom()).isHorizontal()) {
        firstAtom().r[0] = (double)x/toPixels;
        if(fillVolume) {firstAtom().r[0] = 0.0; firstAtom().setDiameter(designTimeXDim);}
        else {firstAtom().setDiameter((double)getSize().height/toPixels);}
        if(boundary) {
            if(y < phase.getSize().height/2) {firstAtom().r[1] = 0.0;}
            else {firstAtom().r[1] = designTimeYDim;}
        }
        else {firstAtom().r[1] = (double)y/toPixels;}
        if(Beans.isDesignTime()) setBounds((int)(Phase.TO_PIXELS*firstAtom().r[0]),(int)(Phase.TO_PIXELS*firstAtom().r[1]),(int)(Phase.TO_PIXELS*firstAtom().getDiameter()),getThickness());
    }
  }

/*    public void draw(Graphics g, int[] origin, double scale){
        firstAtom().draw(g, origin, scale);
    }
*/
}
