package etomica.graphics;
import java.awt.Color;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.plaf.*;
import javax.swing.plaf.metal.*;


/**
 * This class overrides the colors in the DefaultMetalTheme.
 * The fonts and other theme properties are left unaltered.
 * 
 * @author Andrew Walker
 **/
public class EtomicaTheme extends DefaultMetalTheme {

   public String getName() {return "EtomicaTheme";}

   static final Color blush = new Color(153,102,102);
   static final Color brightRed = new Color(153,0,0);
   static final Color darkRed = new Color(102,0,0);
   static final Color tan = new Color(204,204,153);
   static final Color darkKhaki = new Color(102,102,51);
   static final Color khaki = new Color(153,153,102);
   static final Color brightKhaki = khaki.brighter().brighter();


   //
   // Primary Colors (foreground)
   //
   private final ColorUIResource primary1 = new ColorUIResource(tan);
   private final ColorUIResource primary2 = super.getPrimary2();//new ColorUIResource(tan);
   private final ColorUIResource primary3 = super.getPrimary3();//new ColorUIResource(tan);

   protected ColorUIResource getPrimary1() {return primary1;}
   protected ColorUIResource getPrimary2() {return primary2;}
   protected ColorUIResource getPrimary3() {return primary3;}

   //
   // Secondary Colors (background)
   //
   private final ColorUIResource secondary1 = new ColorUIResource(darkKhaki);
   private final ColorUIResource secondary2 = new ColorUIResource(tan);
   private final ColorUIResource secondary3 = new ColorUIResource(brightKhaki);

   protected ColorUIResource getSecondary1() {return secondary1;}
   protected ColorUIResource getSecondary2() {return secondary2;}
   protected ColorUIResource getSecondary3() {return secondary3;}

   //
   // System Colors
   //
   protected final ColorUIResource black = new ColorUIResource(darkRed);
   protected final ColorUIResource white = new ColorUIResource(brightKhaki);

   protected ColorUIResource getBlack() {return black;}
   protected ColorUIResource getWhite() {return white;}

   public ColorUIResource getControl() {
      return new ColorUIResource(tan);}
   public ColorUIResource getControlShadow() {
      return new ColorUIResource(tan.darker());}
   public ColorUIResource getControlDarkShadow() {
      return new ColorUIResource(tan.darker().darker());}
   public ColorUIResource getControlHighlight() {
      return new ColorUIResource(tan.brighter());}
   public ColorUIResource getPrimaryControl() {
      return new ColorUIResource(tan);}
   public ColorUIResource getPrimaryControlShadow() {
      return new ColorUIResource(tan.darker());}
   public ColorUIResource getPrimaryControlHighlight() {
      return new ColorUIResource(tan.brighter());}


   //
   // Other Specific Customizations
   //

   BorderUIResource.CompoundBorderUIResource buttonBorder = new BorderUIResource.CompoundBorderUIResource(BorderFactory.createRaisedBevelBorder(),BorderFactory.createEmptyBorder(2,14,2,14));

   BorderUIResource.LineBorderUIResource titledBorder = new BorderUIResource.LineBorderUIResource(darkRed);

   BorderUIResource.BevelBorderUIResource fieldBorder = new BorderUIResource.BevelBorderUIResource(1,tan.brighter(),tan.darker());

   public void addCustomEntriesToTable(UIDefaults table) {

      super.addCustomEntriesToTable(table);
      table.put("Panel.background",new ColorUIResource(brightKhaki));
      table.put("Panel.foreground",new ColorUIResource(darkRed));
      table.put("Button.focus",new ColorUIResource(darkRed));
      table.put("Button.background",new ColorUIResource(darkRed));
      table.put("Button.foreground",new ColorUIResource(brightKhaki));
      table.put("Button.border",buttonBorder);
      table.put("ComboBox.background",new ColorUIResource(darkRed));
      table.put("ComboBox.foreground",new ColorUIResource(brightKhaki));

      table.put("TextField.inactiveBackground",new ColorUIResource(tan));
      table.put("TextField.border",fieldBorder);
      table.put("control",new ColorUIResource(tan));
      table.put("controlShadow",new ColorUIResource(tan.darker()));
      table.put("controlDkShadow",new ColorUIResource(tan.darker().darker()));
      table.put("controlHighlight",new ColorUIResource(tan.brighter()));
      table.put("controlLtLighligh",new ColorUIResource(tan.brighter().brighter()));
      table.put("TabbedPane.foreground",new ColorUIResource(brightKhaki));
      table.put("TabbedPane.background",new ColorUIResource(darkRed));
      table.put("TabbedPane.highlight",new ColorUIResource(darkRed.brighter()));
      table.put("TabbedPane.lightHighlight",new ColorUIResource(darkRed.brighter().brighter()));
      table.put("TabbedPane.shadow",new ColorUIResource(darkRed.darker()));
      table.put("TabbedPane.darkShadow",new ColorUIResource(darkRed.darker().darker()));
      table.put("TabbedPane.selected",new ColorUIResource(tan));
      table.put("TabbedPane.selectHighlight",new ColorUIResource(tan.brighter()));
      table.put("Slider.foreground",new ColorUIResource(darkRed));
      table.put("Slider.focus",new ColorUIResource(tan));
      table.put("Slider.background",new ColorUIResource(brightKhaki));
      //      table.put("Slider.highlight",new ColorUIResource(darkRed.brighter()));
      //      table.put("Slider.shadow",new ColorUIResource(darkRed.darker()));

      table.put("TitledBorder.border",titledBorder);

      /*      table.put("TitledBorder.border",border);
      table.put("Button.foreground",buttonTextColor1);
      table.put("Button.background",contrast1);
      table.put("ComboBox.foreground",buttonTextColor1);
      table.put("ComboBox.background",contrast1);
      table.put("ComboBox.selectionBackground",panelColor1);
      table.put("ComboBox.selectionForeground",contrast1);
      table.put("Panel.background",background1);
      table.put("Slider.background",background1);
      table.put("Slider.foreground",contrast1);
      table.put("Slider.highlight",contrast1.brighter());
      table.put("Slider.shadow",contrast1.darker());
      table.put("TitledBorder.titleColor",contrast1);
      table.put("ScrollBar.background",background1);
      table.put("ScrollBar.foreground",contrast1);
      table.put("ScrollBar.thumb",contrast1);
      table.put("ScrollBar.thumbHighlight",contrast1.brighter());
      table.put("ScrollBar.thumbShadow",contrast1.darker());
      table.put("TextField.inactiveBackground",panelColor1);
      table.put("control",panelColor1);
      table.put("controlHighlight",panelColor1.brighter());
      table.put("controlShadow",panelColor1.darker());
      table.put("LtHighlight",contrast1.brighter().brighter());
      table.put("DkShadow",contrast1.darker());
      table.put("TabbedPane.foreground",buttonTextColor1);
      table.put("TabbedPane.background",contrast1);
      table.put("TabbedPane.highlight",contrast1.brighter());
      table.put("TabbedPane.shadow",contrast1.darker());
      table.put("textText",contrast1);
      table.put("TextField.border",border1);
      table.put("Button.border",border2);
      table.put("ComboBox.border",border2);
      table.put("Slider.border",background1);
      */
   }	 

}  // EtomicaTheme
