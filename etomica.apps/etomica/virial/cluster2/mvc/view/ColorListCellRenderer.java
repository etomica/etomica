/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.*;

import javax.swing.*;
import javax.swing.border.*;

public class ColorListCellRenderer extends JLabel implements ListCellRenderer {

  private static final long serialVersionUID = 3711349950393326027L;

  private ColorEntry colorEntry;

  public ColorListCellRenderer(Dimension preferredSize) {

    super();
    setOpaque(true);
    setVerticalAlignment(CENTER);
    setHorizontalAlignment(CENTER);
    setFont(new Font("sans", Font.PLAIN, 14));
    setPreferredSize(preferredSize);
    setBorder(new CompoundBorder(new MatteBorder(2, 8, 2, 8, Color.white), new LineBorder(Color.black)));
  }

  public Component getListCellRendererComponent(JList list, Object obj, int row, boolean sel, boolean hasFocus) {

    if (obj instanceof ColorEntry) {
      colorEntry = ((ColorEntry) obj);
    }
    return this;
  }

  public void paint(Graphics g) {

    setForeground(Color.white);
    setText(String.format("(%s)", colorEntry.getText()));
    setBackground(colorEntry.getColor());
    super.paint(g);
  }
}