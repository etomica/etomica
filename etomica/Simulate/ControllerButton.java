package simulate;

import java.awt.*;
import java.awt.event.*;
import java.beans.Beans;

public class ControllerButton extends Controller {

  public Button button;
  public boolean autoStart;

  public ControllerButton() {
    super();
    button = new Button("Start");
    button.setBounds(0,0,60,40);
	button.setBackground(new Color(12632256));
    add(button);
    autoStart = false;
  }

  public void add(Integrator i) {
    super.add(i);
    button.addMouseListener(i);
    if(autoStart) {i.mouseClicked(new MouseEvent(button,MouseEvent.MOUSE_CLICKED,1L,0,0,0,0,false));}
  }
  
  public void setAutoStart(boolean b) {autoStart = b;}
  public boolean getAutoStart() {return autoStart;}
  
  public void setButtonColor(Color c) {button.setBackground(c);}
  public Color getButtonColor() {return button.getBackground();}
  
  public void setButtonTextColor(Color c) {button.setForeground(c);}
  public Color getButtonTextColor() {return button.getForeground();}
  
  public void run() {;}
}


