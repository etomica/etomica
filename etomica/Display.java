package simulate;
import java.awt.*;
import java.beans.Beans;
    
    public class Display extends Container implements simulate.IntegrationIntervalListener{

	View view;
	Phase phase;
    int pixels = 300;
    Image offScreen;
    Graphics osg;

    public Display () {
        setSize(pixels, pixels);
        setBackground(Color.white);
    }

	public void add(View v) {
//	    super.add(v);
	    view = v;
	    view.parentDisplay = this;
	}
	
    public void createOffScreen () {
        if (offScreen == null) {
            offScreen = createImage(pixels, pixels);
            osg = offScreen.getGraphics();
        }
    }
    
    public void paint(Graphics g) {
      if(Beans.isDesignTime()) {
        g.setColor(Color.red);
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
      } 
      createOffScreen();
      view.paint(osg);
      g.drawImage(offScreen, 0, 0, null);
    }

    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
        view.updateView();
        repaint();
    }

}
