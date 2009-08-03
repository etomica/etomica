package etomica.virial.cluster2.ui;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.*;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;

import javax.swing.*;

import org.apache.batik.bridge.UpdateManager;
import org.apache.batik.swing.*;
import org.apache.batik.swing.svg.JSVGComponent;
import org.apache.batik.svggen.*;
import org.apache.batik.transcoder.Transcoder;
import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.svg2svg.SVGTranscoder;
import org.apache.batik.dom.svg.SVGDOMImplementation;

import org.w3c.dom.*;
import org.w3c.dom.svg.*;

import etomica.virial.cluster2.ui.rasterizer.DestinationType;
import etomica.virial.cluster2.ui.rasterizer.Main;
import etomica.virial.cluster2.ui.rasterizer.SVGConverter;
import etomica.virial.cluster2.ui.rasterizer.SVGConverterController;
import etomica.virial.cluster2.ui.rasterizer.SVGConverterFileSource;

public class SVGDraw {

  SVGState state;
  JSVGCanvas canvas;
  SVGContextMenu menu;

  public static void main(String[] args) {

    SVGDraw d = new SVGDraw();
    d.display();
  }

  public SVGDraw() {

    state = SVGDocumentBuilder.getState();
  }

  private String getExtension(String fileName) {

    String ext = null;
    String s = fileName;
    int i = s.lastIndexOf('.');
    if (i > 0 && i < s.length() - 1) {
      ext = s.substring(i).toLowerCase();
    }
    return ext;
  }

  private String getMimeType(String fileName) {

    String ext = getExtension(fileName);
    if (ext.equals(DestinationType.JPEG_EXTENSION)) {
      return DestinationType.JPEG_STR;
    }
    if (ext.equals(DestinationType.PDF_EXTENSION)) {
      return DestinationType.PDF_STR;
    }
    if (ext.equals(DestinationType.PNG_EXTENSION)) {
      return DestinationType.PNG_STR;
    }
    if (ext.equals(DestinationType.TIFF_EXTENSION)) {
      return DestinationType.TIFF_STR;
    }
    return null;
  }

  public void imageOut(String fileName) {

    File f = null;
    try {
      f = File.createTempFile("tmp_", "svg");
      FileOutputStream fs = new FileOutputStream(f);
      streamOut(fs);
      fs.flush();
      fs.close();
      String mimeType = getMimeType(fileName);
      String[] options = new String[13];
      // destination
      options[0] = Main.CL_OPTION_OUTPUT;
      options[1] = fileName;
      // destination type
      options[2] = Main.CL_OPTION_MIME_TYPE;
      options[3] = mimeType;
      // quality (for JPEG)
      options[4] = Main.CL_OPTION_QUALITY;
      options[5] = "0.8";
      // resolution 
      options[6] = Main.CL_OPTION_DPI;
      options[7] = "300";
      // resolution 
      options[8] = Main.CL_OPTION_WIDTH;
      options[9] = "1024";
      // resolution 
      options[10] = Main.CL_OPTION_BACKGROUND_COLOR;
      options[11] = "255.255.255.255";
      // input file
      options[12] = f.getAbsolutePath();
      // convert
      Main c = new Main(options);
      c.execute();
      // remove temporary file
      f.delete();
    }
    catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  public void streamOut(OutputStream svg) {

    // Stream out SVG to the standard output using UTF-8 encoding.
    try {
      TranscoderInput input = new TranscoderInput(state.document);
      TranscoderOutput output = new TranscoderOutput(new BufferedWriter(
          new OutputStreamWriter(svg, "UTF-8")));
      Transcoder t = new SVGTranscoder();
      t.transcode(input, output);
    }
    catch (TranscoderException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch (UnsupportedEncodingException uee) {
      // TODO Auto-generated catch block
      uee.printStackTrace();
    }
  }

  public void display() {

    SVGGraphics2D g2 = new SVGGraphics2D(state.document);
    g2.getRoot(state.document.getDocumentElement());
    canvas = new JSVGCanvas();
    canvas.setDocumentState(JSVGComponent.ALWAYS_DYNAMIC);
    JSVGScrollPane spane = new JSVGScrollPane(canvas);
    menu = new SVGContextMenu(spane, this);
    canvas.setSVGDocument(state.document);
    canvas.addMouseMotionListener(state.getMouseMotionListener(canvas));
    canvas.addMouseListener(state.getMouseClickedListener(menu, canvas));
    JFrame f = new JFrame();
    f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
// f.setLayout(new GridLayout());
    f.setContentPane(new JSVGScrollPane(canvas));
    f.setBounds(new Rectangle(0, 0, SVGProperties.MAP_WIDTH,
        SVGProperties.MAP_HEIGHT + 20));
    f.setVisible(true);
//    f.pack();
// f.setBounds(new Rectangle(0, 0, 300, 180));
  }

  public void prepareSave() {

    state.HOVERING_RECTANGLE_ID = -1;
    state.SELECTED_RECTANGLE_ID = -1;
    final Runnable selected = new SelectedRectangle(state);
    UpdateManager um = canvas.getUpdateManager();
    try {
      um.getUpdateRunnableQueue().invokeAndWait(selected);
    }
    catch (InterruptedException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}

class SVGState {

  final boolean CONTROL_HOVER = true;
  SVGDocument document = null;
  Node HOVERING_RECTANGLE = null;
  int HOVERING_RECTANGLE_ID = -1;
  Node SELECTED_RECTANGLE = null;
  int SELECTED_RECTANGLE_ID = -1;
  final List<Rectangle2D> tiles;
  private Runnable hovering;
  private Runnable selected;

  public SVGState(final SVGDocument doc, final List<Rectangle2D> rectangles) {

    document = doc;
    tiles = rectangles;
    hovering = new HoveringRectangle(this);
    selected = new SelectedRectangle(this);
  }

  private int getTile(Point p) {

    return (p.x / SVGProperties.CLUSTER_WIDTH)
        + (p.y / SVGProperties.CLUSTER_HEIGHT) * SVGProperties.CLUSTERS_ACROSS;
  }

  public MouseListener getMouseClickedListener(final SVGContextMenu menu,
      final JSVGCanvas canvas) {

    return new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

        if (e.getButton() == MouseEvent.BUTTON3) {
          menu.setLocation(e.getXOnScreen(), e.getYOnScreen());
          menu.setVisible(true);
        }
        else if (e.getButton() == MouseEvent.BUTTON1) {
          int newSelected = getTile(e.getPoint());
          if (newSelected >= tiles.size()) {
            return;
          }
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

  public MouseMotionListener getMouseMotionListener(final JSVGCanvas canvas) {

    return new MouseAdapter() {

      @Override
      public void mouseMoved(MouseEvent e) {

        if (CONTROL_HOVER && (e.getButton() == MouseEvent.NOBUTTON)) {
          int newHover = getTile(e.getPoint());
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
            canvas.getUpdateManager().getUpdateRunnableQueue().invokeLater(
                hovering);
          }
        }
      }
    };
  }
}

interface SVGProperties {

  static final int GUTTER = 10;
  static final int NODE_RADIUS = 12;
  static final int GRAPH_RADIUS = 60;
  static final int CLUSTERS_ACROSS = 2;
  static final int MAP_WIDTH = 300;
  static final int MAP_HEIGHT = 140;
  static final int CLUSTER_WIDTH = MAP_WIDTH / CLUSTERS_ACROSS;
  static final int CLUSTER_HEIGHT = CLUSTER_WIDTH;
  static final boolean STREAM_OUT = false;
  static final double INITIAL_ANGLE = Math.PI;
  static final Stroke SELECTED_STROKE = new BasicStroke(1.5f,
      BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,
      new float[] { 3, 3 }, 6);
  static final Color SELECTED_COLOR = Color.DARK_GRAY;
  static final Stroke HOVERING_STROKE = new BasicStroke(1.5f);
  static final Color HOVERING_COLOR = Color.MAGENTA;
}

class SVGDocumentBuilder {

  SVGDocument doc;
  List<Rectangle2D> tiles = new ArrayList<Rectangle2D>();

  private SVGDocumentBuilder() {

    // Create an SVG document.
    DOMImplementation impl = SVGDOMImplementation.getDOMImplementation();
    String svgNS = SVGDOMImplementation.SVG_NAMESPACE_URI;
    doc = (SVGDocument) impl.createDocument(svgNS, "svg", null);
    SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(doc);
    ctx.setComment("Generated by Etomica Cluster Viewer using Apache Batik.");
    ctx.setEmbeddedFontsOn(true);
    Element root = doc.getDocumentElement();
    root.setAttributeNS(null, "width", String.valueOf(SVGProperties.MAP_WIDTH));
    root.setAttributeNS(null, "height", String
        .valueOf(SVGProperties.MAP_HEIGHT));
    root
        .setAttributeNS(null, "style",
            "fill:white;fill-opacity:1;stroke:black;stroke-width:0.5;stroke-opacity:1");
    drawGraph(doc, Color.red, 20, 20);
    drawGraph(doc, Color.blue, 160, 20);
  }

  public static SVGState getState() {

    SVGDocumentBuilder b = new SVGDocumentBuilder();
    return new SVGState(b.doc, b.tiles);
  }

  private Node drawGraph(Document doc, Color color, int dX, int dY) {

    int size = 6;
    int r = SVGProperties.NODE_RADIUS;
    int R = SVGProperties.GRAPH_RADIUS;
    int gutter = SVGProperties.GUTTER;
    double[] x = new double[size];
    double[] y = new double[size];
    double minX = -1, minY = -1, maxX = -1, maxY = -1;
    double angle = 2.0 * Math.PI / size;
    Shape node = new Ellipse2D.Double(0, 0, r, r);
    // draw the cluster on an SVG canvas
    SVGGraphics2D g = new SVGGraphics2D(doc);
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON);
    g.setRenderingHint(RenderingHints.KEY_RENDERING,
        RenderingHints.VALUE_RENDER_QUALITY);
    // compute node centers
    for (int i = 0; i < size; i++) {
      x[i] = (R - r) * Math.cos(SVGProperties.INITIAL_ANGLE + i * angle)
          + (R - r + dX);
      y[i] = (R - r) * Math.sin(SVGProperties.INITIAL_ANGLE + i * angle)
          + (R - r + dY);
      if (minX == -1) {
        minX = x[i];
        maxX = x[i];
        minY = y[i];
        maxY = y[i];
      }
      if (minX > x[i]) {
        minX = x[i];
      }
      if (maxX < x[i]) {
        maxX = x[i];
      }
      if (minY > y[i]) {
        minY = y[i];
      }
      if (maxY < y[i]) {
        maxY = y[i];
      }
    }
    // compute the bounding rectangle
    tiles.add(new Rectangle2D.Double(minX - gutter, minY - gutter, maxX - minX
        + r + 2 * gutter, maxY - minY + r + 2 * gutter));
    // draw edges
    for (int i = 0; i < size; i++) {
      for (int j = i + 1; j < size; j++) {
        Point2D pi = new Point2D.Double(x[i] + r / 2, y[i] + r / 2);
        Point2D pj = new Point2D.Double(x[j] + r / 2, y[j] + r / 2);
        g.setPaint(Color.black);
        g.draw(new Line2D.Double(pi, pj));
      }
    }
    // draw nodes
    for (int i = 0; i < size; i++) {
      g.translate(x[i], y[i]);
      g.setPaint(Color.black);
      g.draw(node);
      g.setPaint(color);
      g.fill(node);
      g.translate(-x[i], -y[i]);
    }
    Element root = doc.getDocumentElement();
    Node child = g.getTopLevelGroup().getFirstChild();
    return root.appendChild(child);
  }
}

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

class SelectedRectangle extends PaintRunnable {

  SelectedRectangle(SVGState state) {

    super(state, SVGProperties.SELECTED_STROKE, SVGProperties.SELECTED_COLOR);
  }

  public void run() {

    clearHovering();
    clearSelection();
    if (state.SELECTED_RECTANGLE_ID != -1) {
      state.SELECTED_RECTANGLE = drawRectangle(state.SELECTED_RECTANGLE_ID);
    }
  }
}

class HoveringRectangle extends PaintRunnable {

  HoveringRectangle(SVGState state) {

    super(state, SVGProperties.HOVERING_STROKE, SVGProperties.HOVERING_COLOR);
  }

  public void run() {

    clearHovering();
    if (state.HOVERING_RECTANGLE_ID != -1) {
      state.HOVERING_RECTANGLE = drawRectangle(state.HOVERING_RECTANGLE_ID);
    }
  }
}