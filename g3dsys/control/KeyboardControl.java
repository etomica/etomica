package g3dsys.control;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

/**
 *	Capture keystrokes to control rotation and translation of the display 
 */

class KeyboardControl implements KeyListener {
	
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
	}

}
