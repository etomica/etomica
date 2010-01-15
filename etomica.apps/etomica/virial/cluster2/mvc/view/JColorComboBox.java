package etomica.virial.cluster2.mvc.view;

import java.util.List;
import java.awt.*;
import javax.swing.*;

public class JColorComboBox extends JComboBox {

  private static final long serialVersionUID = 2510594771259282093L;

  public JColorComboBox(List<ColorEntry> colors) {

    for (ColorEntry entry : colors) {
      addItem(entry);
    }
    setRenderer(new ColorListCellRenderer(new Dimension(40, 24)));
  }
}