package simulate;
import java.util.*;
import java.awt.*;
import java.beans.*;

public class SpeciesBoundaryCubic extends Species {
  SpeciesWall[] wall;
  int thickness;
  
  public SpeciesBoundaryCubic() {
    super();
    name = "BoundaryCubic";
    nElements = 4;
    thickness = 4;
    speciesIndex = 1;
    wall = new SpeciesWall[nElements];
    wall[0] = new SpeciesWall(speciesIndex);
    for(int i = 1; i < nElements; i++) {
        wall[i] = new SpeciesWall(speciesIndex);
        wall[i-1].setNext(wall[i]);
        wall[i].setPrevious(wall[i-1]);
    }
    firstElement = wall[0];
    lastElement = wall[nElements-1];
  }

  // Exposed Properties
  public void setThickness(int t) {thickness = t;}
  public int getThickness() {return thickness;}
    
  public void setSpeciesIndex(int index) {
    this.speciesIndex = index;
    for(int i = 0; i < nElements; i++) {
        wall[i].setSpeciesIndex(index);
    }
  }
  // overrides super class' set method
  public void setNElements(int nElements) {;}
    
  // Class methods

  public void draw(Graphics g, int[] origin, double scale) {
    for(int i = 0; i < nElements; i++) {
        wall[i].draw(g, origin, scale);
    }
  }
  
  public void initializeSpecies(Phase phase) {
    setBounds(0,0,phase.getSize().width, phase.getSize().height);
    int xmax = getSize().width;
    int ymax = getSize().height;
    wall[0].setHorizontal(true);
    wall[1].setVertical(true);
    wall[2].setHorizontal(true);
    wall[3].setVertical(true);
    wall[0].setBounds(0,0,xmax,thickness);  //top
    wall[1].setBounds(xmax-thickness,0,thickness,ymax);  //right
    wall[2].setBounds(0,ymax-thickness,xmax,thickness); //bottom
    wall[3].setBounds(0,0,thickness,ymax);  //left
    for(int i=0; i<nElements; i++) {
        wall[i].initializeSpecies(phase);
    }
  }

}


