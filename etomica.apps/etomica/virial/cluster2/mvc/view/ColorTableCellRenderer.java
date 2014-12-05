/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.mvc.view;

import java.awt.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.TableCellRenderer;

public class ColorTableCellRenderer extends JLabel implements TableCellRenderer {

  private static final long serialVersionUID = 8243082920702301671L;

  private ColorEntry colorEntry;

  public ColorTableCellRenderer() {

    super();
    setOpaque(true);
    setVerticalAlignment(CENTER);
    setHorizontalAlignment(CENTER);
    setFont(new Font("sans", Font.PLAIN, 14));
    setBorder(new CompoundBorder(new MatteBorder(2, 8, 2, 8, Color.white), new LineBorder(Color.black)));
  }

  public Component getTableCellRendererComponent(JTable table, Object obj, boolean isSelected,
      boolean hasFocus, int row, int column) {

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