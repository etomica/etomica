package etomica.virial.cluster2.mvc.view;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Stroke;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;

import javax.swing.JComponent;

import org.apache.batik.bridge.UpdateManager;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.swing.JSVGCanvas;
import org.apache.batik.swing.JSVGScrollPane;
import org.apache.batik.swing.svg.JSVGComponent;
import org.apache.batik.transcoder.Transcoder;
import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.svg2svg.SVGTranscoder;

import etomica.virial.cluster2.graph.GraphSet;
import etomica.virial.cluster2.ui.rasterizer.DestinationType;
import etomica.virial.cluster2.ui.rasterizer.Main;

public class SVGDraw {

  static final Stroke SELECTED_STROKE = new BasicStroke(1.5f,
      BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 1.0f,
      new float[] { 3, 5 }, 3);
  static final Color SELECTED_COLOR = Color.DARK_GRAY;
  static final Stroke HOVERING_STROKE = new BasicStroke(1.5f);
  static final Color HOVERING_COLOR = Color.MAGENTA;
  SVGState state;
  JSVGCanvas canvas;
  private JSVGScrollPane spane;

  public JComponent getPanel() {

    return spane;
  }

  public SVGDraw(GraphSet gs) {

    state = SVGClusterMap.getState(gs);
    SVGGraphics2D g2 = new SVGGraphics2D(state.document);
    g2.getRoot(state.document.getDocumentElement());
    canvas = new JSVGCanvas();
    canvas.setDocumentState(JSVGComponent.ALWAYS_DYNAMIC);
    spane = new JSVGScrollPane(canvas);
    SVGContextMenu menu = new SVGContextMenu(spane, this);
    canvas.setSVGDocument(state.document);
    canvas.addMouseMotionListener(state.getMouseMotionListener(menu, canvas));
    canvas.addMouseListener(state.getMouseClickedListener(menu, canvas));
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
      String[] options = new String[13];
      // destination
      options[0] = Main.CL_OPTION_OUTPUT;
      options[1] = fileName;
      // destination type
      options[2] = Main.CL_OPTION_MIME_TYPE;
      options[3] = getMimeType(fileName);
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