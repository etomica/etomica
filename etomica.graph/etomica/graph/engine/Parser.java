/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine;

import java.util.Map;

public interface Parser {

  // public command names
  public static final String COMMAND_DISPLAY = "display";
  public static final String COMMAND_CLEAR = "clear";
  public static final String COMMAND_DROP = "drop";
  public static final String COMMAND_LIST = "list";
  public static final String COMMAND_PRINT = "print";
  public static final String COMMAND_QUIT = "quit";
  public static final String COMMAND_READ = "read";
  public static final String COMMAND_READDB = "readdb";
  public static final String COMMAND_RUN = "run";
  public static final String COMMAND_SAVE = "save";
  public static final String COMMAND_SET = "set";
  public static final String COMMAND_SUMMARY = "summary";
  public static final String COMMAND_WRITE = "write";
  public static final String COMMAND_WRITEDB = "writedb";

  // public constructor names
  public static final String CONSTRUCTOR_MONO = "mono";
  public static final String CONSTRUCTOR_COLORED = "colored";
  public static final String CONSTRUCTOR_ISOMONO = "iso_mono";
  public static final String CONSTRUCTOR_ISOCOLORED = "iso_colored";

  // public filter names
  public static final String FILTER_HAS_ARTICULATION_POINT = "has_articulation_point";
  public static final String FILTER_HAS_ARTICULATION_PAIR = "has_articulation_pair";
  public static final String FILTER_HAS_NODAL_POINT = "has_nodal_point";
  public static final String FILTER_HAS_NO_ROOT_EDGE = "has_no_root_edge";
  public static final String FILTER_IS_BICONNECTED = "is_biconnected";
  public static final String FILTER_IS_CONNECTED = "is_connected";

  // public unary operation names
  public static final String UNARY_OP_NDIF = "ndif";
  public static final String UNARY_OP_EDIF = "edif";
  public static final String UNARY_OP_EXP = "exp";
  public static final String UNARY_OP_INT = "int";
  public static final String UNARY_OP_ISO = "iso";
  public static final String UNARY_OP_NCOPY = "ncopy";
  public static final String UNARY_OP_PCOPY = "pcopy";
  public static final String UNARY_OP_POWER = "pow";
  public static final String UNARY_OP_RELABEL = "relabel";
  public static final String UNARY_OP_SPLIT = "split";

  // public binary operation names
  public static final String BINARY_OP_CONV = "conv";
  public static final String BINARY_OP_DEL = "del";
  public static final String BINARY_OP_MUL = "mul";
  public static final String BINARY_OP_SUB = "sub";
  public static final String BINARY_OP_SUM = "sum";
  public static final String BINARY_OP_UNION = "union";

  // the parser parses a single source statement given as a string
  public Statement parse(String source) throws ParserException;

  // parser specific exception type
  public class ParserException extends Exception {

    private static final long serialVersionUID = 941457687073478706L;

    public ParserException(String message) {

      super(message);
    }

    public ParserException(String message, Throwable cause) {

      super(message, cause);
    }
  }

  // Statement ::= Assignment | Command
  public interface Statement {

  }

  // Assignment ::= Variable '=' StrictExpression
  public interface Assignment extends Statement {

    public StrictExpression getExpression();

    public Variable getVariable();
  }

  // Command
  public interface Command extends Statement {

    public String getCommand();
  }

  // CommandPropertyValue
  public interface CommandPropertyValue extends Command {

    public Property getProperty();

    public Value getValue();
  }

  // CommandValue
  public interface CommandValue extends Command {

    public Value getValue();
  }

  // CommandVariableValue
  public interface CommandVariable extends Command {

    public Variable getVariable();
  }

  // CommandVariableValue
  public interface CommandVariableValue extends Command {

    public Value getValue();

    public Variable getVariable();
  }

  // Expression ::= Variable | StrictExpression
  public interface Expression {

  }

  // Variable
  public interface Variable extends Expression {

    public String getName();
  }

  // StrictExpression ::= Constructor ['>' filter_name]* | UnaryOp | BinaryOp
  public interface StrictExpression extends Expression {

  }

  // Constructor
  public interface Constructor extends StrictExpression {

    public void addFilter(String filterName);

    public boolean hasFilter(String filterName);
  }

  // ConstructorMono
  public interface ConstructorMono extends Constructor {

    public byte getFieldNodes();

    public byte getRootNodes();

    public boolean isIsoFree();
  }

  // ConstructorColored
  public interface ConstructorColored extends Constructor {

    public Map<Character, Byte> getFieldColorMap();

    public Map<Character, Byte> getRootColorMap();

    public boolean isIsoFree();
  }

  // UnaryOp
  public interface UnaryOp extends StrictExpression {

    public Expression getExpression();

    public String getOperation();
  }

  // UnaryOpColor
  public interface UnaryOpColor extends UnaryOp {

    public char getColor1();
  }

  // UnaryOpByte
  public interface UnaryOpByte extends UnaryOp {

    public byte getByte1();
  }

  // UnaryOpByteMap
  public interface UnaryOpByteMap extends UnaryOp {

    public Map<Byte, Byte> getColorMap();
  }

  // UnaryOpThreeColor
  public interface UnaryOpThreeColor extends UnaryOp {

    public char getColor1();

    public char getColor2();

    public char getColor3();
  }

  // UnaryOpTwoByte
  public interface UnaryOpTwoByte extends UnaryOp {

    public byte getByte1();

    public byte getByte2();
  }

  // BinaryOp
  public interface BinaryOp extends StrictExpression {

    public Expression getExpression1();

    public Expression getExpression2();

    public String getOperation();
  }

  // BinaryOpByte
  public interface BinaryOpByte extends BinaryOp {

    public byte getByte1();
  }

  // Property
  public interface Property {

    public String getName();
  }

  // Value
  public interface Value {

    public String getValue();
  }
}