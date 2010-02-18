package etomica.graph.traversal;

/*
 * BitmapUtils are used by the traversal algorithms.
 *
 * @author Demian Lessa
 */
public class BitmapUtils {

  protected static final byte SZ_INT = 32;

  public static int bitMask(byte bitSize) {

    return (0xFFFFFFFF >>> (SZ_INT - bitSize));
  }

  public static int bitOnMask(byte bit) {

    return (0x00000001 << bit);
  }

  public static int bitOffMask(byte bit) {

    return ~bitOnMask(bit);
  }

  public static byte leftmostBit(int map) {

    for (byte i = 0; i < SZ_INT; i++) {
      if ((map & bitOnMask((byte) i)) > 0) {
        return i;
      }
    }
    return -1;
  }
}