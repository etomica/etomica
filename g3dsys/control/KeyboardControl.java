package g3dsys.control;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

/**
 *	Capture keystrokes to control rotation and translation of the display 
 */

class KeyboardControl implements KeyListener {
	
  /*
   * Keyboard controls:
   * arrow keys - translation
   * =/- - zoom up/down
   * h - home position; removes rotation, translation, and zoom, and rescales 
   * x/X - x axis rotation
   * y/Y - y axis rotation
   * z/Z - z axis rotation
   * r - redraw display
   * p - toggle perspective on/off
   * i - toggle image shell on/off
   * 1-9 - set the number of shells to displays
   * b - cycle through boundary draw styles
   */
  
  private G3DSys master;
	
  public KeyboardControl(G3DSys m) { master = m; }

  public void keyPressed(KeyEvent e) {
    if( e.getKeyCode() == KeyEvent.VK_DOWN ) {
      master.xlateDown();
      master.fastRefresh();
    }
    if( e.getKeyCode() == KeyEvent.VK_UP ) {
      master.xlateUp();
      master.fastRefresh();
    }
    if( e.getKeyCode() == KeyEvent.VK_LEFT ) {
      master.xlateLeft();
      master.fastRefresh();
    }
    if( e.getKeyCode() == KeyEvent.VK_RIGHT ) {
      master.xlateRight();
      master.fastRefresh();
    }
  }
	
  public void keyReleased(KeyEvent e) {}
	
  public void keyTyped(KeyEvent e) {
    if( e.getKeyChar() == 'r' ) {
      master.refresh();
    }
    if( e.getKeyChar() == '=' ) {
      master.zoomUp(10);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '-' ) {
      master.zoomDown(10);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'z' ) {
      master.rotateByZ(-10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'Z' ) {
      master.rotateByZ(10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'x' ) {
      master.rotateByX(-10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'X' ) {
      master.rotateByX(10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'y' ) {
      master.rotateByY(-10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'Y' ) {
      master.rotateByY(10f);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'h' ) {
      master.rotateToHome();
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'p' ) {
      master.setPerspectiveDepth(!master.getPerspectiveDepth());
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'i' ) {
      master.setEnableImages(!master.isEnableImages());
      master.fastRefresh();
    }
    if( e.getKeyChar() == '1' ) {
      master.setLayers(1);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '2' ) {
      master.setLayers(2);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '3' ) {
      master.setLayers(3);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '4' ) {
      master.setLayers(4);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '5' ) {
      master.setLayers(5);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '6' ) {
      master.setLayers(6);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '7' ) {
      master.setLayers(7);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '8' ) {
      master.setLayers(8);
      master.fastRefresh();
    }
    if( e.getKeyChar() == '9' ) {
      master.setLayers(9);
      master.fastRefresh();
    }
    if( e.getKeyChar() == 'b' ) {
      master.cycleDrawBoundaryType();
      master.fastRefresh();
    }
  }

}
