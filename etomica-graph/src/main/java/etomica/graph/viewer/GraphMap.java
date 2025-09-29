/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.viewer;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.HeadlessException;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;

import org.apache.batik.dom.svg.SAXSVGDocumentFactory;
import org.apache.batik.swing.JSVGCanvas;
import org.apache.batik.swing.JSVGScrollPane;
import org.apache.batik.transcoder.Transcoder;
import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.svg2svg.SVGTranscoder;
import org.apache.batik.util.XMLResourceDescriptor;
import org.w3c.dom.svg.SVGDocument;

import etomica.graph.model.Coefficient;
import etomica.graph.model.Graph;
import etomica.graph.view.rasterizer.DestinationType;
import etomica.graph.view.rasterizer.Main;

public class GraphMap {

  private static final int defaultBorder = 10;
  private int defaultGutter = 10;
  private int clusterDim = 200;
  private int clusterDimOut = clusterDim + 2 * defaultGutter;
  private boolean doIncludeCoefficients = false;

  private JSVGCanvas canvas;

  private int clustersAcross;
  private int clustersAlong;
  private int height;
  private JComponent container;
  private Set<Graph> graphs;

  private SVGDocument svg;

  private int width;

  public GraphMap(Set<Graph> graphs, JComponent container) {

    this.container = container;
    this.graphs = graphs;
    init();
  }
  
  public void setDoIncludeCoefficients(boolean doIncludeCoefficients) {
    this.doIncludeCoefficients = doIncludeCoefficients;
    defaultGutter = doIncludeCoefficients ? 25 : 10;
    clusterDimOut = clusterDim + 2 * defaultGutter;
    init();
  }
  
  protected void init() {
    clustersAcross = (container.getWidth() - 4 * defaultBorder) / (clusterDimOut);
    clustersAcross = clustersAcross > graphs.size() ? graphs.size() : clustersAcross;
    clustersAlong = clustersAcross > 0 ? (graphs.size() / clustersAcross + (graphs.size() % clustersAcross > 0 ? 1 : 0)) : 0;
    width = 4 * defaultBorder + clustersAcross * (clusterDimOut);
    height = 2 * defaultBorder + clustersAlong * (clusterDimOut);
  }

  public void draw() {

    String mapSVG = String
        .format(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" contentStyleType=\"text/css\" width=\"%dpx\" height=\"%dpx\" viewBox=\"0 0 %d %d\" preserveAspectRatio=\"xMinYMin slice\">\n",
            width, height, width, height);
    mapSVG += "<g>\n";
    mapSVG += String.format(
        "<rect style=\"fill: white; stroke: none\" x=\"0\" y=\"0\" width=\"%d\" height=\"%d\"/>\n", width,
        height);

    int col = 0;
    int row = 0;
    String gMapSVG = "";
    for (Graph g : graphs) {
      if (doIncludeCoefficients) {
        Coefficient c = g.coefficient();
        if (c.getDenominator() == 1) {
          String cStr = "";
          if (c.getNumerator() > 0) {
            if (row+col>0) {
              cStr = "+";
            }
            if (c.getNumerator() != 1) {
              cStr += Integer.toString(c.getNumerator());
            }
          }
          else {
            cStr = "-";
            if (c.getNumerator() != -1) {
              cStr = Integer.toString(c.getNumerator());
            }
          }
          gMapSVG += String.format("<text x=\"%d\" y=\"%d\">%s</text>", col*clusterDimOut+defaultBorder, row*clusterDimOut+clusterDimOut/2+defaultBorder, cStr);
        }
        else {
          if (row+col>0 || c.getNumerator() < 0) {
            String sStr = c.getNumerator() > 0 ? "+" : "-";
            gMapSVG += String.format("<text x=\"%d\" y=\"%d\">%s</text>", col*clusterDimOut+defaultBorder, row*clusterDimOut+clusterDimOut/2+defaultBorder, sStr);
          }
          gMapSVG += String.format("<text x=\"%d\" y=\"%d\">%d</text>", col*clusterDimOut+defaultBorder+15, row*clusterDimOut+clusterDimOut/2+defaultBorder-10, Math.abs(c.getNumerator()));
          gMapSVG += String.format("<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\"/>", col*clusterDimOut+defaultBorder+10, row*clusterDimOut+clusterDimOut/2+defaultBorder-5, col*clusterDimOut+defaultBorder+30, row*clusterDimOut+clusterDimOut/2+defaultBorder-5);
          gMapSVG += String.format("<text x=\"%d\" y=\"%d\">%d</text>", col*clusterDimOut+defaultBorder+15, row*clusterDimOut+clusterDimOut/2+defaultBorder+10, c.getDenominator());
        }
      }
      gMapSVG += String
          .format(
              "<g id=\"g%d\" style=\"fill:none; stroke: black; stroke-width: 1.1;\" transform=\"translate(%d,%d)\">\n",
              row * clustersAcross + col, col * clusterDimOut + defaultBorder + defaultGutter, row * clusterDimOut + defaultBorder + defaultGutter);
      gMapSVG += "<title>"+g.toString()+"</title>\n";
      gMapSVG += g.toSVG(clusterDim);
      gMapSVG += "</g>\n";
      col++;
      if (col % clustersAcross == 0) {
        col = 0;
        row++;
        mapSVG += gMapSVG;
        gMapSVG = "";
      }
    }
    mapSVG += gMapSVG;
    mapSVG += "</g></svg>";

    File f = null;
    try {
      f = File.createTempFile("tmp_", ".svg");
      FileOutputStream fos = new FileOutputStream(f);
      try {
        fos.write(mapSVG.getBytes());
        fos.flush();
      }
      finally {
        fos.close();
      }
    }
    catch (IOException e) {
      e.printStackTrace();
      return;
    }
    if (graphs.size() > 5000) {
      System.out.println("Too many graphs to open a display.  SVG File is "+f.getAbsolutePath());
      return;
    }

    try {
      String parser = XMLResourceDescriptor.getXMLParserClassName();
      SAXSVGDocumentFactory df = new SAXSVGDocumentFactory(parser);
      svg = df.createSVGDocument(f.toURI().toString());
    }
    catch (IOException ex) {
      // do your error handling here
    }

    canvas = new JSVGCanvas();
    // tooltips don't work with ALWAYS_STATIC...
//    canvas.setDocumentState(JSVGComponent.ALWAYS_STATIC);
    JSVGScrollPane spane = new JSVGScrollPane(canvas);
    SVGContextMenu menu = new SVGContextMenu(spane);
    canvas.setSVGDocument(svg);
    canvas.addMouseMotionListener(getMouseMotionListener(menu));
    canvas.addMouseListener(getMouseClickedListener(menu));
    // add the pane to the frame
    container.add(spane, BorderLayout.CENTER);
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
      e.printStackTrace();
    }
  }

  public void streamOut(OutputStream out) {

    // Stream out SVG to the standard output using UTF-8 encoding.
    try {
      TranscoderInput input = new TranscoderInput(svg);
      TranscoderOutput output = new TranscoderOutput(new BufferedWriter(new OutputStreamWriter(out, "UTF-8")));
      Transcoder t = new SVGTranscoder();
      t.transcode(input, output);
    }
    catch (TranscoderException e) {
      e.printStackTrace();
    }
    catch (UnsupportedEncodingException uee) {
      uee.printStackTrace();
    }
  }

  private MouseListener getMouseClickedListener(final SVGContextMenu menu) {

    return new MouseAdapter() {

      @Override
      public void mouseClicked(MouseEvent e) {

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
      }
    };
  }

  private MouseMotionListener getMouseMotionListener(final SVGContextMenu menu) {

    return new MouseMotionListener() {

      public void mouseMoved(MouseEvent e) {

        if (e.getButton() == MouseEvent.NOBUTTON) {
          menu.setVisible(false);
        }
      }

      public void mouseDragged(MouseEvent e) {

        // no-op
      }
    };
  }

  class SVGContextMenu extends JPopupMenu {

    private static final long serialVersionUID = -5605994099188052122L;
    public static final float EXPORT_WIDTH = 1024f;
    private Component owner;

    public SVGContextMenu(Component parent) {

      super();
      this.owner = parent;
      JMenuItem menuSaveAs = new JMenuItem("Save as...");
      add(menuSaveAs);
      MouseListener popupListener = new SaveAsListener();
      menuSaveAs.addMouseListener(popupListener);
    }

    /**
     * This listener exports the SVGCanvas to an appropriate file based on the selection
     * of destination file name defined by the user.
     * 
     * @author Demian Lessa
     */
    class SaveAsListener extends MouseAdapter {

      private void save(File file) throws Exception {

        String ext = SVGFileFilter.getExtension(file);
        if (ext == null) {
          return;
        }
        if (ext.equals(SVGFileFilter.EXT_svg)) {
          writeOut(file, new SVGTranscoder());
        }
        else if (ext.equals(SVGFileFilter.EXT_gif)) {
          File f = File.createTempFile("tmp_", ".png");
          imageOut(f.getAbsolutePath());
          // read the temporary PNG file into an image
          BufferedImage image = ImageIO.read(f);
          // write the image to disk
          ImageIO.write(image, ext, file);
          f.delete();
        }
        else {
          imageOut(file.getAbsolutePath());
        }
      }

      // for SVG only
      private void writeOut(File file, Transcoder t) throws Exception {

        // Set the transcoder input and output.
        TranscoderInput input = new TranscoderInput(svg);
        OutputStream ostream = new FileOutputStream(file);
        Writer w;
        if (t instanceof SVGTranscoder) {
          w = new OutputStreamWriter(ostream, "UTF-8");
        }
        else {
          w = new OutputStreamWriter(ostream);
        }
        TranscoderOutput output = new TranscoderOutput(w);
        // Transcode
        t.transcode(input, output);
        ostream.flush();
        ostream.close();
      }

      // make sure that after selecting the "save as" option, we clear all
      // adornments on the SVGCanvas- namely, the selection and hover rectangles
      public void mouseClicked(MouseEvent e) {

        canvas.setEnabled(false);
        setVisible(false);
        JFileChooser fc = new SVGFileChooser();
        if (fc.showSaveDialog(owner) == JFileChooser.APPROVE_OPTION) {
          File f = fc.getSelectedFile().getAbsoluteFile();
          if (f.exists()) {
            int actionDialog = JOptionPane.showConfirmDialog(owner, "Replace existing file?");
            if (actionDialog != JOptionPane.YES_OPTION) {
              return;
            }
          }
          try {
            save(fc.getSelectedFile().getAbsoluteFile());
            // System.out.println("Saved file: " +
            // fc.getSelectedFile().getAbsoluteFile());
          }
          catch (Exception e1) {
            // cowardly ignore
            e1.printStackTrace();
          }
        }
        canvas.setEnabled(true);
      }
    }
  }

  /**
   * Filter with all supported destination formats for exporting the image in the
   * SVGCanvas.
   * 
   * @author Demian Lessa
   */
  static class SVGFileFilter extends javax.swing.filechooser.FileFilter {

    public final static String EXT_gif = "gif";
    public final static String EXT_jpg = "jpg";
    public final static String EXT_pdf = "pdf";
    public final static String EXT_png = "png";
    public final static String EXT_svg = "svg";
    public final static String EXT_tif = "tif";

    private final static String DESC_gif = "GIF bitmap (*.gif)";
    private final static String DESC_jpg = "JPEG bitmap (*.jpeg)";
    private final static String DESC_pdf = "PDF document (*.pdf)";
    private final static String DESC_png = "PNG bitmap (*.png)";
    private final static String DESC_svg = "SVG vector graphic (*.svg)";
    private final static String DESC_tif = "TIFF bitmap (*.tif)";

    private final static SVGFileFilter FILTER_GIF = new SVGFileFilter(EXT_gif, DESC_gif);
    private final static SVGFileFilter FILTER_JPG = new SVGFileFilter(EXT_jpg, DESC_jpg);
    private final static SVGFileFilter FILTER_PDF = new SVGFileFilter(EXT_pdf, DESC_pdf);
    private final static SVGFileFilter FILTER_PNG = new SVGFileFilter(EXT_png, DESC_png);
    private final static SVGFileFilter FILTER_SVG = new SVGFileFilter(EXT_svg, DESC_svg);
    private final static SVGFileFilter FILTER_TIFF = new SVGFileFilter(EXT_tif, DESC_tif);

    private String extension;
    private String description;

    public SVGFileFilter(String ext, String desc) {

      super();
      extension = ext;
      description = desc;
    }

    public static String getExtension(File f) {

      String ext = null;
      String s = f.getName();
      int i = s.lastIndexOf('.');
      if (i > 0 && i < s.length() - 1) {
        ext = s.substring(i + 1).toLowerCase();
      }
      return ext;
    }

    @Override
    public boolean accept(File f) {

      if (f.isDirectory()) {
        return true;
      }
      String ext = getExtension(f);
      return (ext != null) && (ext.equals(extension));
    }

    @Override
    public String getDescription() {

      return description;
    }

    public String getExtension() {

      return extension;
    }
  }

  /**
   * Customized JFileChooser aware of the SVGCanvas export filters.
   * 
   * @author Demian Lessa
   */
  static class SVGFileChooser extends JFileChooser {

    private static final long serialVersionUID = -9020881522405108590L;

    public SVGFileChooser() {

      super();
      setDialogTitle("Export Cluster Diagrams");
      setAcceptAllFileFilterUsed(false);
      addChoosableFileFilter(SVGFileFilter.FILTER_GIF);
      addChoosableFileFilter(SVGFileFilter.FILTER_JPG);
      addChoosableFileFilter(SVGFileFilter.FILTER_PDF);
      addChoosableFileFilter(SVGFileFilter.FILTER_PNG);
      addChoosableFileFilter(SVGFileFilter.FILTER_SVG);
      addChoosableFileFilter(SVGFileFilter.FILTER_TIFF);
      setFileFilter(SVGFileFilter.FILTER_SVG);
    }

    @Override
    public int showSaveDialog(Component component) throws HeadlessException {

      int result = super.showSaveDialog(component);
      if (result == JFileChooser.APPROVE_OPTION) {
        File sf = getSelectedFile();
        String ext = SVGFileFilter.getExtension(sf);
        String fext = ((SVGFileFilter) getFileFilter()).getExtension();
        if ((ext == null) || (!ext.equals(fext))) {
          File nsf = new File(sf.getAbsolutePath() + "." + fext);
          setSelectedFile(nsf);
        }
      }
      return result;
    }
  }
}