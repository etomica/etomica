import java.awt.*;

import java.net.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import java.awt.geom.*;

import javax.swing.text.View;

import java.util.*;

import java.awt.event.*;

import java.io.*;

import java.awt.print.*;

import java.beans.*;

import javax.swing.text.*;
import javax.swing.text.html.*;
import javax.swing.plaf.basic.*;
import javax.swing.plaf.metal.*;
import javax.swing.plaf.*;


/**
 @Linda
 @Experimental
 @09.09.2001
 */

class BlueRoseTheme extends DefaultMetalTheme {

    public String getName() {
        return "BlueRoses Theme...";
    }

    private final ColorUIResource   primary1            =
        new ColorUIResource(Colors.MBlue);
    private final ColorUIResource   primary2            =
        new ColorUIResource(Colors.lgray);    // highlights menuItems, too
    private final ColorUIResource   primary3            =
        new ColorUIResource(Colors.lgray);
    private final ColorUIResource   secondary1          =
        new ColorUIResource(102, 102, 102);
    private final ColorUIResource   secondary2          =
        new ColorUIResource(128, 128, 128);
    private final ColorUIResource   secondary3          =
        new ColorUIResource(Colors.MBlue);    //Menu, Panelbackground, Java-PrintDialog
    protected final ColorUIResource getControlHighlight =
        new ColorUIResource(240, 240, 240);
    protected final ColorUIResource text                =
        new ColorUIResource(Colors.lgray);
    protected final ColorUIResource light               =
        new ColorUIResource(Color.white);
    protected ColorUIResource       focus               =
        new ColorUIResource(Colors.orange);
    protected ColorUIResource       controlhigh         =
        new ColorUIResource(Colors.galadrelwhite);
    protected ColorUIResource       controlshade        =
        new ColorUIResource(Colors.lgray);
    private FontUIResource          systemtextfont      =
        new FontUIResource("Dialog", Font.BOLD, 14);
    private FontUIResource          menutextfont        =
        new FontUIResource("Dialog", Font.BOLD, 14);
    private FontUIResource          windowtitlefont     =
        new FontUIResource("Serif", Font.ITALIC, 14);
    private FontUIResource          controltextfont     =
        new FontUIResource("Serif", Font.ITALIC, 14);

    protected ColorUIResource getPrimary1() {
        return primary1;
    }

    protected ColorUIResource getPrimary2() {
        return primary2;
    }

    protected ColorUIResource getPrimary3() {
        return primary3;
    }

    protected ColorUIResource getSecondary1() {
        return secondary1;
    }

    protected ColorUIResource getSecondary2() {
        return secondary2;
    }

    protected ColorUIResource getSecondary3() {
        return secondary3;
    }

    public ColorUIResource getControlHighlight() {
        return light;
    }

    //Why the name getBlack? Does it make sense?
    protected ColorUIResource getBlack() {
        return text;
    }

    protected ColorUIResource getWhite() {
        return light;
    }

    public ColorUIResource getFocusColor() {
        return focus;
    }

    public ColorUIResource getPrimaryControlHighlight() {
        return controlhigh;
    }

    public ColorUIResource getPrimaryControlShadow() {
        return controlshade;
    }

    public FontUIResource getMenuTextFont() {
        return menutextfont;
    }

    public FontUIResource getSystemTextFont() {
        return systemtextfont;
    }

    public FontUIResource getWindowTitleFont() {    // JIFTitleFont
        return windowtitlefont;
    }

    public FontUIResource getControlTextFont() {    // e.g.ButtonText
        return controltextfont;
    }

    public void addCustomEntriesToTable(UIDefaults table) {

        super.addCustomEntriesToTable(table);
	table.put("control",new ColorUIResource(Colors.lgray));
        table.put("Button.focus", new ColorUIResource(Colors.lgray));
        table.put("Button.background", new ColorUIResource(Colors.lgray));
        table.put("Button.foreground", new ColorUIResource(Colors.drose));
        table.put("Button.border", BorderFactory.createRaisedBevelBorder());
        table.put("OptionPane.background", new ColorUIResource(Colors.dblue));
        table.put("OptionPane.buttonAreaBorder",
                  BorderFactory.createEmptyBorder(5, 5, 5, 5));
        table.put("OptionPane.messageAreaBorder",
                  BorderFactory.createEmptyBorder(5, 5, 5, 5));
        table.put("Button.font", new FontUIResource(Fonts.font));
        table.put("MenuItem.selectionBackground",
                  new ColorUIResource(Colors.lightrose));
        table.put("Menu.selectionBackground",
                  new ColorUIResource(Colors.lightrose));
    }
}

class SizeListener extends ComponentAdapter {

    static String nl = System.getProperty("line.separator");

    public void componentResized(ComponentEvent e) {
        System.out.println("Note: Frame has been resized: " + nl
                           + e.getComponent());
    }

    public void componentMoved(ComponentEvent e) {
        System.out.println("Note: Frame has been moved: " + nl
                           + e.getComponent());
    }
}

/**
 *My BorderClass
 */
class Borders {

    static Border       orange       =
        BorderFactory.createLineBorder(Colors.MGold);
    static Border       empty        = BorderFactory.createEmptyBorder(5, 5,
                                           5, 5);
    static Border       compound     =
        BorderFactory.createCompoundBorder(orange, empty);
    static Border       space        = BorderFactory.createEmptyBorder(5, 5,
                                           5, 10);
    static Border       rb           =
        BorderFactory.createRaisedBevelBorder();
    static Border       lb           = BorderFactory.createBevelBorder(1,
                                           Colors.drose, Colors.lgray);
    static Border       b2           = BorderFactory.createCompoundBorder(lb,
                                           rb);
    static Border       raisedbevel  =
        BorderFactory.createRaisedBevelBorder();
    static Border       loweredbevel =
        BorderFactory.createLoweredBevelBorder();
    static Border       status       =
        BorderFactory.createCompoundBorder(raisedbevel, loweredbevel);
    static Border       arealine     =
        BorderFactory.createLineBorder(Colors.highlighter);
    static TitledBorder areatitle    =
        BorderFactory.createTitledBorder(arealine, "Source");
    static TitledBorder buttontitle  =
        BorderFactory.createTitledBorder(arealine, "Info");
    static Border       b3           = BorderFactory.createCompoundBorder(rb,
                                           buttontitle);
}

/**
 * ColorClass
 */
class Colors {

    static final Color MRed           = new Color(102, 0, 51);
    static final Color lightrose      = new Color(200, 188, 212);
    static final Color MBlue          = new Color(24, 35, 87);
    static final Color MGold          = new Color(204, 180, 0);
    static final Color MWhite         = new Color(245, 245, 245);
    static final Color dblue          = new Color(63, 64, 124);
    static final Color drose          = new Color(159, 61, 100);
    static final Color lgray          = new Color(105, 120, 141);
    static final Color highlighter    = new Color(169, 31, 93);
    static final Color orange         = new Color(255, 153, 102);
    static final Color titlegold      = new Color(230, 184, 110);
    static final Color galadrelbronze = new Color(188, 154, 142);
    static final Color galadrelwhite  = new Color(193, 184, 180);
}

/**
 * GradientClass
 */
class Gradients {

    static final GradientPaint gradient        = new GradientPaint(0, 0,
                                                     Colors.lgray, 0, 50,
                                                     Colors.lightrose);
    static final GradientPaint gradientdesktop = new GradientPaint(0, 0,
                                                     Colors.lightrose, 0, 50,
                                                     Colors.lgray);
}

/**
 *FontClass
 */
class Fonts {
    static Font font = new Font("Serif", Font.ITALIC, 16);
}
