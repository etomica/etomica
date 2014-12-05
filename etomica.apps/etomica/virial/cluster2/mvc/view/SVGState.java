/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.Color;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.swing.SwingUtilities;

import org.apache.batik.bridge.UpdateManager;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.swing.JSVGCanvas;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.svg.SVGDocument;


/**
 * This class controls the state of a SVGCanvas. It determines which
 * rectangles are selected and hovered using a list of 2D rectangles.
 * 
 * @author Demian Lessa
 */
public class SVGState {

  final boolean CONTROL_HOVER = true;
  final SVGDocument document;
  Node HOVERING_RECTANGLE = null;
  int HOVERING_RECTANGLE_ID = -1;
  Node SELECTED_RECTANGLE = null;
  int SELECTED_RECTANGLE_ID = -1;
  final List<Rectangle2D> tiles;
  private Runnable hovering;
  private Runnable selected;
  private SVGClusterMap svgMap;

  public SVGState(final SVGClusterMap map) {

    svgMap = map;
    document = map.doc;
    tiles = map.tiles;
    hovering = new HoveringRectangle(this);
    selected = new SelectedRectangle(this);
  }

  public MouseListener getMouseClickedListener(final SVGContextMenu menu,
      final JSVGCanvas canvas) {

    return new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

        int newSelected = svgMap.getTile(e.getPoint());
        if (newSelected >= tiles.size()) {
          menu.setVisible(false);
          return;
        }
        Point p = new Point(e.getPoint());
        SwingUtilities.convertPointToScreen(p, canvas);
        if (menu.isVisible()) {
          if (e.getButton() == MouseEvent.BUTTON3) {
            menu.setLocation(p);
            return;
          }
          else {
            menu.setVisible(false);
          }
        }
        if (e.getButton() == MouseEvent.BUTTON3) {
          menu.setLocation(p);
          menu.setVisible(true);
        }
        else if (e.getButton() == MouseEvent.BUTTON1) {
          if (SELECTED_RECTANGLE_ID == newSelected) {
            SELECTED_RECTANGLE_ID = -1;
          }
          else {
            SELECTED_RECTANGLE_ID = newSelected;
          }
          HOVERING_RECTANGLE_ID = -1;
          canvas.getUpdateManager().getUpdateRunnableQueue().invokeLater(
              selected);
        }
      }
    };
  }

  public MouseMotionListener getMouseMotionListener(final SVGContextMenu menu,
      final JSVGCanvas canvas) {

    return new MouseMotionListener() {

      public void mouseMoved(MouseEvent e) {

        if (CONTROL_HOVER && (e.getButton() == MouseEvent.NOBUTTON)) {
          menu.setVisible(false);
          int newHover = svgMap.getTile(e.getPoint());
          if (newHover >= tiles.size()) {
            newHover = -1;
          }
          if (SELECTED_RECTANGLE_ID == newHover) {
            newHover = -1;
          }
          if (newHover != HOVERING_RECTANGLE_ID) {
            if (newHover == -1) {
              HOVERING_RECTANGLE_ID = -1;
            }
            else {
              HOVERING_RECTANGLE_ID = newHover;
            }
            UpdateManager um = canvas.getUpdateManager();
            if (um != null) {
              um.getUpdateRunnableQueue().invokeLater(hovering);
            }
          }
        }
      }

      public void mouseDragged(MouseEvent e) {

        // no-op
      }
    };
  }
}

/**
 * The classes below are Runnable implementations of painting routines for the
 * SVGCanvas. The SVGCanvas canvas follows the Swing approach to painting in
 * that all painting operations must happen on the same thread. Therefore, these
 * runnable objects are queued for execution through the SVGCanvas UpdateManager
 * Queue, which schedules them for execution on the Swing thread. All styles for 
 * painting and drawing are defined as constants at the SVGDraw class.
 * 
 * @author Demian Lessa
 */
abstract class PaintRunnable implements Runnable {

  protected Map<Integer, Node> rectangles = new HashMap<Integer, Node>();
  protected Color color;
  protected Stroke stroke;
  protected SVGState state;

  PaintRunnable(SVGState state, Stroke stroke, Color color) {

    this.state = state;
    this.color = color;
    this.stroke = stroke;
  }

  protected void clearSelection() {

    if (state.SELECTED_RECTANGLE != null) {
      Element root = state.document.getDocumentElement();
      root.removeChild(state.SELECTED_RECTANGLE);
      state.SELECTED_RECTANGLE = null;
    }
  }

  protected void clearHovering() {

    if (state.HOVERING_RECTANGLE != null) {
      Element root = state.document.getDocumentElement();
      root.removeChild(state.HOVERING_RECTANGLE);
      state.HOVERING_RECTANGLE = null;
    }
  }

  protected Node drawRectangle(int i) {

    if (!rectangles.containsKey(i)) {
      // draw the cluster on an SVG canvas
      SVGGraphics2D g = new SVGGraphics2D(state.document);
      g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
          RenderingHints.VALUE_ANTIALIAS_ON);
      g.setRenderingHint(RenderingHints.KEY_RENDERING,
          RenderingHints.VALUE_RENDER_QUALITY);
      g.setStroke(stroke);
      g.setPaint(color);
      g.draw(state.tiles.get(i));
      Node child = g.getTopLevelGroup().getFirstChild();
      rectangles.put(i, child);
    }
    Element root = state.document.getDocumentElement();
    return root.appendChild(rectangles.get(i));
  }
}

/**
 * The purpose of this class is to paint a selected rectangle around a single
 * cluster, clearing up any other existing selected rectangles.
 * 
 * @author Demian Lessa
 */
class SelectedRectangle extends PaintRunnable {

  SelectedRectangle(SVGState state) {

    super(state, SVGDraw.SELECTED_STROKE, SVGDraw.SELECTED_COLOR);
  }

  public void run() {

    clearHovering();
    clearSelection();
    if (state.SELECTED_RECTANGLE_ID != -1) {
      state.SELECTED_RECTANGLE = drawRectangle(state.SELECTED_RECTANGLE_ID);
    }
  }
}

/**
 * The purpose of this class is to paint a hovering rectangle around a single
 * cluster, clearing up any other existing hovering rectangles. No hovering
 * rectangle is painted if the rectangle under the mouse is selected.
 * 
 * @author Demian Lessa
 */
class HoveringRectangle extends PaintRunnable {

  HoveringRectangle(SVGState state) {

    super(state, SVGDraw.HOVERING_STROKE, SVGDraw.HOVERING_COLOR);
  }

  public void run() {

    clearHovering();
    if (state.HOVERING_RECTANGLE_ID != -1) {
      state.HOVERING_RECTANGLE = drawRectangle(state.HOVERING_RECTANGLE_ID);
    }
  }
}