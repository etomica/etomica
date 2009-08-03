package etomica.virial.cluster2.ui;

import java.awt.Color;
import java.awt.Component;
import java.awt.HeadlessException;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;

import org.apache.batik.swing.JSVGCanvas;
import org.apache.batik.swing.JSVGScrollPane;
import org.apache.batik.transcoder.Transcoder;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.image.PNGTranscoder;
import org.apache.batik.transcoder.image.TIFFTranscoder;
import org.apache.batik.transcoder.svg2svg.SVGTranscoder;
import org.apache.fop.svg.PDFTranscoder;
import org.w3c.dom.Document;
import org.w3c.dom.svg.SVGDocument;

import etomica.virial.cluster2.ui.rasterizer.DestinationType;

public class SVGContextMenu extends JPopupMenu {

  private static final long serialVersionUID = -5605994099188052122L;
  public static final float EXPORT_WIDTH = 1024f;
  private JMenuItem menuSaveAs;
  private Component owner;
// private JPopupMenu self;
  private SVGDraw svgDraw;

  public SVGContextMenu(Component parent, SVGDraw d) {

    super();
    owner = parent;
// self = this;
    svgDraw = d;
    menuSaveAs = new JMenuItem("Save as...");
    add(menuSaveAs);
    MouseListener popupListener = new SaveAsListener();
    menuSaveAs.addMouseListener(popupListener);
  }

  class SaveAsListener extends MouseAdapter {

    private void save(File file) throws Exception {

      String ext = SVGFileFilter.getExtension(file);
      if (ext == null) {
        return;
      }
      Transcoder t = getTranscoder(ext);
      if (t == null) {
        System.out.println("no transcoder found for " + file + " ("
            + SVGFileFilter.getExtension(file) + ")");
        return;
      }
      if (ext.equals(SVGFileFilter.EXT_svg)) {
        writeOut(file, t);
      }
      else if (ext.equals(SVGFileFilter.EXT_gif) || ext.equals(SVGFileFilter.EXT_bmp)) {
        File f = File.createTempFile("tmp_", ".png");
        svgDraw.imageOut(f.getAbsolutePath());
        //read the temporary PNG file into an image
        BufferedImage image = ImageIO.read(f);
        // write the image to disk
        ImageIO.write(image, ext, file);
        f.delete();
      }
      else {
        svgDraw.imageOut(file.getAbsolutePath());
      }
// if ((t instanceof PDFTranscoder) || (t instanceof TIFFTranscoder)
// || ext.equals(SVGFileFilter.EXT_png)) {
// ByteArrayOutputStream os = streamOut(t);
// FileOutputStream fos = new FileOutputStream(file);
// fos.write(os.toByteArray());
// fos.flush();
// fos.close();
// }
// else {
// ByteArrayOutputStream os = streamOut(t);
// // create a PNG image from the raw XML source
// BufferedImage image = ImageIO.read(new ByteArrayInputStream(os
// .toByteArray()));
// // write the image to disk
// ImageIO.write(image, ext, file);
// }
    }

    private ByteArrayOutputStream streamOut(Transcoder t) throws Exception {

      // Set the transcoder input and output.
      TranscoderInput input = new TranscoderInput(svgDraw.state.document);
      ByteArrayOutputStream ostream = new ByteArrayOutputStream();
      TranscoderOutput output = new TranscoderOutput(ostream);
      t.transcode(input, output);
      ostream.flush();
      ostream.close();
      return ostream;
    }

    // for SVG, ...
    private void writeOut(File file, Transcoder t) throws Exception {

      // Set the transcoder input and output.
      TranscoderInput input = new TranscoderInput(svgDraw.state.document);
      OutputStream ostream = new FileOutputStream(file);
      Writer w;
      if (t instanceof SVGTranscoder) {
        w = new OutputStreamWriter(ostream, "UTF-8");
      }
      else {
        w = new OutputStreamWriter(ostream);
      }
      TranscoderOutput output = new TranscoderOutput(w);
      // Perform the transcoding.
      t.transcode(input, output);
      ostream.flush();
      ostream.close();
    }

    private Transcoder getTranscoder(String ext) {

      Transcoder t = null;
      if (ext.equals(SVGFileFilter.EXT_pdf)) {
        PDFTranscoder pdf = new PDFTranscoder();
        pdf.addTranscodingHint(PDFTranscoder.KEY_WIDTH, EXPORT_WIDTH);
        t = pdf;
      }
      else if (ext.equals(SVGFileFilter.EXT_svg)) {
        SVGTranscoder svg = new SVGTranscoder();
        t = svg;
      }
      else if (ext.equals(SVGFileFilter.EXT_tif)) {
        TIFFTranscoder tiff = new TIFFTranscoder();
        tiff.addTranscodingHint(TIFFTranscoder.KEY_WIDTH, EXPORT_WIDTH);
        tiff.addTranscodingHint(TIFFTranscoder.KEY_BACKGROUND_COLOR,
            Color.white);
        t = tiff;
      }
      else {
        PNGTranscoder png = new PNGTranscoder();
        png.addTranscodingHint(PNGTranscoder.KEY_BACKGROUND_COLOR, Color.white);
        png.addTranscodingHint(PNGTranscoder.KEY_WIDTH, EXPORT_WIDTH);
        t = png;
      }
      return t;
    }

    public void mouseClicked(MouseEvent e) {

      svgDraw.canvas.setEnabled(false);
      svgDraw.prepareSave();
      setVisible(false);
      JFileChooser fc = new SVGFileChooser();
      if (fc.showSaveDialog(owner) == JFileChooser.APPROVE_OPTION) {
        File f = fc.getSelectedFile().getAbsoluteFile();
        if (f.exists()) {
          int actionDialog = JOptionPane.showConfirmDialog(owner,
              "Replace existing file?");
          if (actionDialog != JOptionPane.YES_OPTION) {
            return;
          }
        }
        try {
          save(fc.getSelectedFile().getAbsoluteFile());
          System.out.println("Saved file: "
              + fc.getSelectedFile().getAbsoluteFile());
        }
        catch (Exception e1) {
          // cowardly ignore
          e1.printStackTrace();
        }
      }
      svgDraw.canvas.setEnabled(true);
    }
  }
}

class SVGFileFilter extends javax.swing.filechooser.FileFilter {

  public final static String EXT_bmp = "bmp";
  public final static String EXT_gif = "gif";
  public final static String EXT_jpg = "jpg";
  public final static String EXT_pdf = "pdf";
  public final static String EXT_png = "png";
  public final static String EXT_svg = "svg";
  public final static String EXT_tif = "tif";
  private final static String DESC_bmp = "BMP bitmap (*.bmp)";
  private final static String DESC_gif = "GIF bitmap (*.gif)";
  private final static String DESC_jpg = "JPEG bitmap (*.jpeg)";
  private final static String DESC_pdf = "PDF document (*.pdf)";
  private final static String DESC_png = "PNG bitmap (*.png)";
  private final static String DESC_svg = "SVG vector graphic (*.svg)";
  private final static String DESC_tif = "TIFF bitmap (*.tif)";
  final static SVGFileFilter FILTER_BMP = new SVGFileFilter(EXT_bmp, DESC_bmp);
  final static SVGFileFilter FILTER_GIF = new SVGFileFilter(EXT_gif, DESC_gif);
  final static SVGFileFilter FILTER_JPG = new SVGFileFilter(EXT_jpg, DESC_jpg);
  final static SVGFileFilter FILTER_PDF = new SVGFileFilter(EXT_pdf, DESC_pdf);
  final static SVGFileFilter FILTER_PNG = new SVGFileFilter(EXT_png, DESC_png);
  final static SVGFileFilter FILTER_SVG = new SVGFileFilter(EXT_svg, DESC_svg);
  final static SVGFileFilter FILTER_TIFF = new SVGFileFilter(EXT_tif, DESC_tif);
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

class SVGFileChooser extends JFileChooser {

  private static final long serialVersionUID = -9020881522405108590L;

  public SVGFileChooser() {

    super();
    setDialogTitle("Export Cluster Diagrams");
    setAcceptAllFileFilterUsed(false);
    addChoosableFileFilter(SVGFileFilter.FILTER_BMP);
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
