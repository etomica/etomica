package etomica.virial.cluster2.ui;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.*;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.List;
import java.util.ArrayList;

import javax.swing.*;

import org.apache.batik.swing.*;
import org.apache.batik.svggen.*;
import org.apache.batik.transcoder.Transcoder;
import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.svg2svg.SVGTranscoder;
import org.apache.batik.dom.svg.SVGDOMImplementation;

import org.w3c.dom.*;
import org.w3c.dom.svg.*;

public class SVGDraw {

  private static Node SELECTED_RECTANGLE = null;
  private static final double INITIAL_ANGLE = Math.PI;
  private static final int GUTTER = 10;
  private static final int GRAPH_RADIUS = 60;
  private static final int NODE_RADIUS = 12;
  private static final boolean STREAM_OUT = true;
  private static List<Rectangle2D> tiles = new ArrayList<Rectangle2D>();

  public static void main(String[] args) {

    // Create an SVG document.
    DOMImplementation impl = SVGDOMImplementation.getDOMImplementation();
    String svgNS = SVGDOMImplementation.SVG_NAMESPACE_URI;
    SVGDocument doc = (SVGDocument) impl.createDocument(svgNS, "svg", null);
    Element root = doc.getDocumentElement();
    root.setAttributeNS(svgNS, "width", "300");
    root.setAttributeNS(svgNS, "height", "140");
    root
        .setAttributeNS(svgNS, "style",
            "fill:white;fill-opacity:1;stroke:black;stroke-width:0.5;stroke-opacity:1");
    drawGraph(doc, Color.red, 20, 20);
    drawGraph(doc, Color.blue, 160, 20);
//    drawRectangle(doc, 0);
//    drawRectangle(doc, 1);
    if (STREAM_OUT) {
      streamOut(doc);
    }
    display(doc);
  }

  private static void display(final SVGDocument doc) {

    SVGGraphics2D g2 = new SVGGraphics2D(doc);
    // Populate the document root with the generated SVG content.
    g2.getRoot(doc.getDocumentElement());
    // Display the document.
    final JSVGScrollPane spane = new JSVGScrollPane(new JSVGCanvas());
    spane.getCanvas().setSVGDocument(doc);
    spane.getCanvas().addMouseListener(new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

        System.out.println("Mouse clicked @" + e.getPoint());
        for (int i = 0; i < tiles.size(); i++) {
          if (tiles.get(i).contains(e.getPoint())) {
            System.out.println("Belongs to rectangle #" + i);
            if (SELECTED_RECTANGLE != null) {
              Element root = doc.getDocumentElement();
              root.removeChild(SELECTED_RECTANGLE);
            }
            SELECTED_RECTANGLE = drawRectangle(doc, i);
            spane.getCanvas().setSVGDocument(doc);
          }
        }
      }
    });
    JFrame f = new JFrame();
    f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    f.getContentPane().add(spane);
    f.setVisible(true);
    f.setBounds(new Rectangle(0, 0, 300, 180));
  }

  private static void streamOut(SVGDocument doc) {

    // Stream out SVG to the standard output using UTF-8 encoding.
    try {
      TranscoderInput input = new TranscoderInput(doc);
      TranscoderOutput output = new TranscoderOutput(new OutputStreamWriter(
          System.out, "UTF-8"));
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

  private static Node drawRectangle(Document doc, int i) {

    // draw the cluster on an SVG canvas
    SVGGraphics2D g = new SVGGraphics2D(doc);
    g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON);
    g.setRenderingHint(RenderingHints.KEY_RENDERING,
        RenderingHints.VALUE_RENDER_QUALITY);
    Stroke stroke = new BasicStroke(1.5f, BasicStroke.CAP_BUTT,
        BasicStroke.JOIN_MITER, 1.0f, new float[] {4, 4}, 8);
    g.setStroke(stroke);
    g.setPaint(Color.lightGray);
    g.draw(tiles.get(i));
    Element root = doc.getDocumentElement();
    Node child = g.getTopLevelGroup().getFirstChild();
    return root.appendChild(child);
  }

  private static Node drawGraph(Document doc, Color color, int dX, int dY) {

    int size = 6;
    int r = NODE_RADIUS;
    int R = GRAPH_RADIUS;
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
      x[i] = (R - r) * Math.cos(INITIAL_ANGLE + i * angle) + (R - r + dX);
      y[i] = (R - r) * Math.sin(INITIAL_ANGLE + i * angle) + (R - r + dY);
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
    tiles.add(new Rectangle2D.Double(minX - GUTTER, minY - GUTTER, maxX - minX
        + r + 2 * GUTTER, maxY - minY + r + 2 * GUTTER));
    System.out.println("Rectangle #" + (tiles.size() - 1) + " = "
        + tiles.get(tiles.size() - 1));
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