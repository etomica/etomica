package simulate;
import java.awt.*;
import java.beans.Beans;

    public abstract class Display extends Canvas implements simulate.IntegrationIntervalListener {

	Simulation parentSimulation;
    int pixels = 200;
    Image offScreen;
    Graphics osg;
    int updateInterval;
    int iieCount;
    Component displayTool = null;

    public Display () {
        setSize(pixels, pixels);
        setBackground(Color.white);
	    setUpdateInterval(1);
    }

    public void createOffScreen () {
        if (offScreen == null) {
            offScreen = createImage(pixels, pixels);
            osg = offScreen.getGraphics();
        }
    }
    
    public void setMeter(Meter m) {;}
    
    public void update(Graphics g) {paint(g);}
    
    public void paint(Graphics g) {
      if(Beans.isDesignTime()) {
        g.setColor(Color.red);
        g.drawRect(0,0,getSize().width-1,getSize().height-1);
        g.drawRect(1,1,getSize().width-3,getSize().height-3);
      } 
      createOffScreen();
      doPaint(osg);
      g.drawImage(offScreen, 0, 0, null);
    }

    public void integrationIntervalAction(IntegrationIntervalEvent evt) {
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        doUpdate();
            repaint();
	    }
    }

    public abstract void doUpdate();
    
    public abstract void doPaint(Graphics g);
    
	public void setParentSimulation(Simulation s) {parentSimulation = s;}

    public final int getUpdateInterval() {return updateInterval;}
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
}
