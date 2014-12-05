/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.engine.impl;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import etomica.graph.engine.Parser;

public class ParserImpl implements Parser {

  public static final boolean DEBUG_MODE = false;

  // parsing symbols consumed during parsing
  private static final String SYMBOL_LPARENS = "(";
  private static final String SYMBOL_RPARENS = ")";
  private static final String SYMBOL_LBRACKET = "{";
  private static final String SYMBOL_RBRACKET = "}";
  private static final String SYMBOL_AT = "@";
  private static final String SYMBOL_DOLLAR = "$";
  private static final String SYMBOL_COLUMN = ":";
  private static final String SYMBOL_COMMA = ",";
  private static final String SYMBOL_EQUAL = "=";
  private static final String SYMBOL_GT = ">";
  private static final String SYMBOL_QUOTE = "'";

  // parsing context: used only for parsing
  private static final String CONTEXT_STATEMENT = "STATEMENT/";
  private static final String CONTEXT_ASSIGNMENT = "ASSIGNMENT/";
  private static final String CONTEXT_COMMAND = "COMMAND/";
  private static final String CONTEXT_STRICT_EXPRESSION = "STRICT_EXPRESSION/";
  private static final String CONTEXT_CONSTRUCTOR = "CONSTRUCTOR/";
  private static final String CONTEXT_CONSTRUCTOR_FILTER = "CONSTRUCTOR_FILTER/";
  private static final String CONTEXT_BINARY_OP = "BINARY_OP/";
  private static final String CONTEXT_UNARY_OP = "UNARY_OP/";
  private static final String CONTEXT_BYTE_MAP = "BYTE_MAP/";
  private static final String CONTEXT_CHAR_MAP = "CHAR_MAP/";
  private static final String CONTEXT_ROOT = "/";

  // delimiter list: used only for parsing
  private static final String WHITESPACE = " \t\f\r\n";
  private static final String DELIMITERS = WHITESPACE + SYMBOL_LPARENS + SYMBOL_RPARENS + SYMBOL_LBRACKET
      + SYMBOL_RBRACKET + SYMBOL_AT + SYMBOL_DOLLAR + SYMBOL_COLUMN + SYMBOL_COMMA + SYMBOL_EQUAL + SYMBOL_GT
      + SYMBOL_QUOTE;

  // command lists: used only for parsing
  private static final List<String> IDENT_COMMANDS = Arrays.asList(new String[] { COMMAND_DISPLAY,
      COMMAND_DROP, COMMAND_PRINT, COMMAND_SUMMARY });
  private static final List<String> ARITY0_COMMANDS = Arrays
      .asList(new String[] { COMMAND_CLEAR, COMMAND_LIST, COMMAND_QUIT });
  private static final List<String> ARITY1_COMMANDS = Arrays.asList(new String[] { COMMAND_DISPLAY,
      COMMAND_DROP, COMMAND_PRINT, COMMAND_READDB, COMMAND_RUN, COMMAND_SAVE, COMMAND_SUMMARY,
      COMMAND_WRITEDB });
  private static final List<String> ARITY2_COMMANDS = Arrays.asList(new String[] { COMMAND_READ,
      COMMAND_WRITE, COMMAND_SET });

  // constructor list: used only for parsing
  private static final List<String> CONSTRUCTORS = Arrays.asList(new String[] { CONSTRUCTOR_MONO,
      CONSTRUCTOR_COLORED, CONSTRUCTOR_ISOMONO, CONSTRUCTOR_ISOCOLORED });

  // filter list: used only for parsing
  private static final List<String> FILTERS = Arrays.asList(new String[] { FILTER_HAS_ARTICULATION_PAIR,
      FILTER_HAS_ARTICULATION_POINT, FILTER_HAS_NO_ROOT_EDGE, FILTER_HAS_NODAL_POINT, FILTER_IS_BICONNECTED,
      FILTER_IS_CONNECTED });

  // unary operation list: used only for parsing
  private static final List<String> UNARY_OPS = Arrays.asList(new String[] { UNARY_OP_NDIF, UNARY_OP_EDIF,
      UNARY_OP_EXP, UNARY_OP_INT, UNARY_OP_ISO, UNARY_OP_NCOPY, UNARY_OP_PCOPY, UNARY_OP_POWER,
      UNARY_OP_RELABEL, UNARY_OP_SPLIT });

  // binary operation list: used only for parsing
  private static final List<String> BINARY_OPS = Arrays.asList(new String[] { BINARY_OP_CONV, BINARY_OP_DEL,
      BINARY_OP_MUL, BINARY_OP_SUB, BINARY_OP_SUM, BINARY_OP_UNION });

  private String token;
  private String context;
  private StringTokenizer tokenizer;
  private boolean isQString;

  public Statement parse(String source) throws ParserException {

    isQString = false;
    tokenizer = new StringTokenizer(source, DELIMITERS, true);
    context = CONTEXT_ROOT;
    Statement result = statement();
    if (DEBUG_MODE) {
      System.out.println(context);
    }
    return result;
  }

  private String scanQString() {

    String next = "";
    String result = "";
    while (tokenizer.hasMoreTokens() && !next.equals(SYMBOL_QUOTE)) {
      next = tokenizer.nextToken();
      if (!next.equals(SYMBOL_QUOTE)) {
        result += next;
      }
    }
    if (tokenizer.hasMoreTokens()) {
      return result;
    }
    else {
      return null;
    }
  }

  // advances token to the next non-whitespace token and returns true only if such a token
  // exists; if it does not exist, token is set to null
  private boolean next() {

    if (isQString) {
      token = scanQString();
      return token != null;
    }
    if (tokenizer.hasMoreTokens()) {
      token = tokenizer.nextToken();
      if (WHITESPACE.contains(token)) {
        return next();
      }
      context += token + "/";
      return true;
    }
    token = null;
    return false;
  }

  // advances only if the next token is not null
  private boolean advance() throws ParserException {

    check(next());
    return true;
  }

  // advances only if the next token is the expected string
  private boolean advance(String expected) throws ParserException {

    check(next(), expected);
    return true;
  }

  // advances only if the next token is either of the strings
  private boolean advance(String one, String other) throws ParserException {

    check(next(), one, other);
    return true;
  }

  // advances only if the next token is one of the strings in the expected set
  private boolean advance(List<String> expected) throws ParserException {

    check(next(), expected);
    return true;
  }

  private boolean check(boolean expression) throws ParserException {

    return expression || die(String.format("Syntax error. Invalid token: %s.", token));
  }

  private boolean check(boolean expression, String expected) throws ParserException {

    return (expression && token.equals(expected))
        || die(String.format("Syntax error. Expected token: %s. Found token: %s.", expected, token));
  }

  private boolean check(boolean expression, String one, String other) throws ParserException {

    return (expression && (token.equals(one) || token.equals(other)))
        || die(String.format("Syntax error. Expected one of the tokens: %s or %s. Found token: %s.", one,
            other, token));
  }

  private boolean check(boolean expression, List<String> expected) throws ParserException {

    return (expression && expected.contains(token))
        || die(String.format("Syntax error. Expected one of the tokens: %s. Found token: %s.", expected
            .toString(), token));
  }

  private boolean die(String message) throws ParserException {

    if (DEBUG_MODE) {
      System.out.println("Context: " + context);
    }
    throw new ParserException(message);
  }

  private boolean die(Throwable cause) throws ParserException {

    if (DEBUG_MODE) {
      System.out.println("Context: " + context);
    }
    throw new ParserException(String.format(
        "Syntax error. Invalid token (%s) while parsing statement at context %s.", token, context), cause);
  }

  private Statement statement() throws ParserException {

    context += CONTEXT_STATEMENT;
    advance();
    if (token.equals(SYMBOL_DOLLAR)) {
      return assignment();
    }
    else {
      return command();
    }
  }

  private Assignment assignment() throws ParserException {

    context += "/" + CONTEXT_ASSIGNMENT;
    advance();
    Variable variable = variable();
    advance(SYMBOL_EQUAL);
    advance();
    return new AssignmentImpl(variable, strictExpression());
  }

  private StrictExpression strictExpression() throws ParserException {

    context += CONTEXT_STRICT_EXPRESSION;
    StrictExpression result = null;
    String strict = token.toLowerCase();
    if (CONSTRUCTORS.contains(strict)) {
      result = constructor(strict);
    }
    else if (UNARY_OPS.contains(strict)) {
      result = unaryOperator(strict);
    }
    else if (BINARY_OPS.contains(strict)) {
      result = binaryOperation(strict);
    }
    if (result != null) {
      return result;
    }
    die(String.format("Syntax error. StrictExpression expected but found '%s' instead.", strict));
    return null;
  }

  private StrictExpression binaryOperation(String strict) throws ParserException {

    context += CONTEXT_BINARY_OP;
    BinaryOp operation;
    advance(SYMBOL_LPARENS);
    advance();
    Expression expr1 = expression(token.equals(SYMBOL_DOLLAR));
    check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
    advance();
    Expression expr2 = expression(token.equals(SYMBOL_DOLLAR));
    // no additional parameters
    if (!strict.equals(BINARY_OP_CONV)) {
      operation = new BinaryOpImpl(strict, expr1, expr2);
    }
    // single byte parameter
    else {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      byte param1 = typeByte();
      next();
      operation = new BinaryOpByteImpl(strict, expr1, expr2, param1);
    }
    check(token.equals(SYMBOL_RPARENS), SYMBOL_RPARENS);
    // move passed the last token of the binary operator
    next();
    return operation;
  }

  private StrictExpression unaryOperator(String strict) throws ParserException {

    context += CONTEXT_UNARY_OP;
    UnaryOp operation;
    advance(SYMBOL_LPARENS);
    advance();
    Expression expression = expression(token.equals(SYMBOL_DOLLAR));
    // no additional parameters
    if (strict.equals(UNARY_OP_ISO) || strict.equals(UNARY_OP_NCOPY) || strict.equals(UNARY_OP_PCOPY)) {
      operation = new UnaryOpImpl(strict, expression);
    }
    // single char parameter
    else if (strict.equals(UNARY_OP_NDIF) || strict.equals(UNARY_OP_EDIF)) {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      char param1 = typeColor();
      next();
      operation = new UnaryOpColorImpl(strict, expression, param1);
    }
    // single byte parameter
    else if (strict.equals(UNARY_OP_INT) || strict.equals(UNARY_OP_POWER)) {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      byte param1 = typeByte();
      next();
      operation = new UnaryOpByteImpl(strict, expression, param1);
    }
    // single byte map parameter
    else if (strict.equals(UNARY_OP_RELABEL)) {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      Map<Byte, Byte> param1 = typeByteMap();
      next();
      operation = new UnaryOpByteMapImpl(strict, expression, param1);
    }
    // two byte parameters
    else if (strict.equals(UNARY_OP_EXP)) {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      byte param1 = typeByte();
      advance(SYMBOL_COMMA);
      advance();
      byte param2 = typeByte();
      next();
      operation = new UnaryOpTwoByteImpl(strict, expression, param1, param2);
    }
    // three char parameters
    else {
      check(token.equals(SYMBOL_COMMA), SYMBOL_COMMA);
      advance();
      char param1 = typeColor();
      advance(SYMBOL_COMMA);
      advance();
      char param2 = typeColor();
      advance(SYMBOL_COMMA);
      advance();
      char param3 = typeColor();
      next();
      operation = new UnaryOpThreeColorImpl(strict, expression, param1, param2, param3);
    }
    check(token.equals(SYMBOL_RPARENS), SYMBOL_RPARENS);
    // move passed the last token of the unary operator
    next();
    return operation;
  }

  private StrictExpression constructor(String strict) throws ParserException {

    context += CONTEXT_CONSTRUCTOR;
    Constructor constructor;
    if (strict.equals(CONSTRUCTOR_MONO) || strict.equals(CONSTRUCTOR_ISOMONO)) {
      advance(SYMBOL_LPARENS);
      advance();
      byte param1 = typeByte();
      advance(SYMBOL_COMMA);
      advance();
      byte param2 = typeByte();
      advance(SYMBOL_RPARENS);
      next();
      constructor = new ConstructorMonoImpl(param1, param2, strict.equals(CONSTRUCTOR_ISOMONO));
    }
    else {
      advance(SYMBOL_LPARENS);
      advance();
      Map<Character, Byte> param1 = typeColorMap();
      advance(SYMBOL_COMMA);
      advance();
      Map<Character, Byte> param2 = typeColorMap();
      advance(SYMBOL_RPARENS);
      next();
      constructor = new ConstructorColoredImpl(param1, param2, strict.equals(CONSTRUCTOR_ISOCOLORED));
    }
    if (token != null && token.equals(SYMBOL_GT)) {
      context += CONTEXT_CONSTRUCTOR_FILTER;
      do {
        check(token.equals(SYMBOL_GT), SYMBOL_GT);
        advance(FILTERS);
        constructor.addFilter(token);
        next();
      } while (token != null && token.equals(SYMBOL_GT));
    }
    // move passed the last token of the constructor
    return constructor;
  }

  private Expression expression(boolean isVariable) throws ParserException {

    Expression result;
    if (isVariable) {
      // consume the $
      advance();
      result = variable();
      // move to the next unprocessed token
      next();
    }
    else {
      result = strictExpression();
    }
    return result;
  }

  private Command command() throws ParserException {

    context += CONTEXT_COMMAND;
    Command result = null;
    String command = token.toLowerCase();
    if (ARITY0_COMMANDS.contains(command)) {
      result = new CommandImpl(command);
    }
    if (ARITY1_COMMANDS.contains(command)) {
      advance(SYMBOL_LPARENS);
      if (IDENT_COMMANDS.contains(command)) {
        advance(SYMBOL_DOLLAR);
        advance();
        Variable variable = variable();
        advance(SYMBOL_RPARENS);
        result = new CommandVariableImpl(command, variable);
      }
      else {
        Value qstring = qstring();
        advance(SYMBOL_RPARENS);
        result = new CommandValueImpl(command, qstring);
      }
    }
    if (ARITY2_COMMANDS.contains(command)) {
      advance(SYMBOL_LPARENS);
      advance(SYMBOL_DOLLAR, SYMBOL_AT);
      Variable variable = null;
      Property property = null;
      if (token.equals(SYMBOL_DOLLAR)) {
        advance();
        variable = variable();
      }
      else {
        advance();
        property = property();
      }
      advance(SYMBOL_COMMA);
      Value qstring = qstring();
      advance(SYMBOL_RPARENS);
      if (variable != null) {
        result = new CommandVariableValueImpl(command, variable, qstring);
      }
      else {
        result = new CommandPropertyValueImpl(command, property, qstring);
      }
    }
    next();
    if (token == null && result != null) {
      return result;
    }
    if (token != null) {
      die(String.format("Syntax error. Invalid token after Command: '%s'.", token));
    }
    else {
      die(String.format("Syntax error. Command expected but found '%s' instead.", command));
    }
    return null;
  }

  private Variable variable() throws ParserException {

    return new VariableImpl(token);
  }

  private Property property() throws ParserException {

    return new PropertyImpl(token);
  }

  private Value qstring() throws ParserException {

    advance(SYMBOL_QUOTE);
    isQString = true;
    advance();
    Value result = new ValueImpl(token);
    isQString = false;
    return result;
  }

  private Byte typeByte() throws ParserException {

    try {
      return Byte.parseByte(token);
    }
    catch (NumberFormatException nfe) {
      die(nfe);
    }
    return null;
  }

  private Character typeColor() throws ParserException {

    check(token.length() == 1 && Character.isLetter(token.charAt(0)));
    return token.charAt(0);
  }

  private Map<Byte, Byte> typeByteMap() throws ParserException {

    byte lhs;
    byte rhs;
    Map<Byte, Byte> result = new HashMap<Byte, Byte>();
    context += CONTEXT_BYTE_MAP;
    check(token.equals(SYMBOL_LBRACKET), SYMBOL_LBRACKET);
    do {
      advance();
      lhs = typeByte();
      advance(SYMBOL_COLUMN);
      advance();
      rhs = typeByte();
      result.put(lhs, rhs);
      advance();
    } while (token.equals(SYMBOL_COMMA));
    check(token.equals(SYMBOL_RBRACKET), SYMBOL_RBRACKET);
    return result;
  }

  private Map<Character, Byte> typeColorMap() throws ParserException {

    char lhs;
    byte rhs;
    Map<Character, Byte> result = new HashMap<Character, Byte>();
    context += CONTEXT_CHAR_MAP;
    check(token.equals(SYMBOL_LBRACKET), SYMBOL_LBRACKET);
    do {
      advance();
      lhs = typeColor();
      advance(SYMBOL_COLUMN);
      advance();
      rhs = typeByte();
      result.put(lhs, rhs);
      advance();
    } while (token.equals(SYMBOL_COMMA));
    check(token.equals(SYMBOL_RBRACKET), SYMBOL_RBRACKET);
    return result;
  }
}