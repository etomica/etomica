/* Copyright (C) 2002-2005  The Jmol Development Team
 *
 * Contact: jmol-developers@lists.sf.net
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package g3dsys.control;

import java.awt.Component;
import java.awt.Event;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

class MouseManager implements MouseListener, MouseMotionListener {

  private G3DSys gsys;
  private Component c; // for giving focus on click

  private int previousDragX, previousDragY;
  
  public MouseManager(Component component, G3DSys gs) {
    component.addMouseListener(this);
    component.addMouseMotionListener(this);
    c = component;
    gsys = gs;
  }

  final static int LEFT = 16;
  final static int RIGHT = Event.META_MASK;  // 4

  private void mouseSinglePressDrag(int deltaX, int deltaY, int modifiers) {
    switch (modifiers) {
    case LEFT:
      gsys.rotateByY((float)deltaX);
      gsys.rotateByX((float)deltaY);
      gsys.fastRefresh();
      break;
    case RIGHT:
    	gsys.xlateX(deltaX);
    	gsys.xlateY(deltaY);
    	gsys.fastRefresh();
    	break;
    }
  }

  public void mousePressed(MouseEvent e) {
    c.requestFocus(); // Sun says this is not foolproof
    previousDragX = e.getX();
    previousDragY = e.getY();
  }
  
  public void mouseDragged(MouseEvent e) {
    int deltaX = e.getX() - previousDragX;
    int deltaY = e.getY() - previousDragY;
    previousDragX = e.getX(); previousDragY = e.getY();
    mouseSinglePressDrag(deltaX/2, deltaY/2, e.getModifiers());
  }

  public void mouseMoved(MouseEvent e) {}
  public void mouseEntered(MouseEvent e) {}
  public void mouseExited(MouseEvent e) {}
  public void mouseClicked(MouseEvent e) {}
  public void mouseReleased(MouseEvent e) {}
  
}
