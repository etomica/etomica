package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;

public class Controller extends Container {

  public Integrator integrator;
  public Button button;
  public Phase phase;
  public Simulation parentSimulation;

  public Controller() {
    setSize(100,40);
    button = new Button("Start");
    button.setBounds(0,0,60,40);
	button.setBackground(new Color(12632256));
    add(button);
  }

  public void add(Integrator i) {
    super.add(i);
    this.integrator = i;
    parentSimulation.registerPhases(i);
    button.addMouseListener(i);
  }
  
  public void paint(Graphics g) {
    if(Beans.isDesignTime()) {
      g.setColor(Color.red);
      g.drawRect(0,0,getSize().width-1,getSize().height-1);
      g.drawRect(1,1,getSize().width-3,getSize().height-3);
    } 
    g.setColor(getBackground());
    paintComponents(g);
  }
}


