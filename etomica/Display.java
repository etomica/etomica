package simulate;
import java.awt.*;
    
    public class Display extends Container {

        int pixels = 300;

        public Display () {
            setSize(pixels, pixels);
            setBackground(Color.white);
        }

        Image offScreen;
        Graphics osg;

        public void createOffScreen () {
            if (offScreen == null) {
                offScreen = createImage(pixels, pixels);
                osg = offScreen.getGraphics();
            }
        }        
    }
