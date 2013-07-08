import scala.util.parsing.combinator.syntactical._


object GLParser extends StandardTokenParsers {
  val keywords = List("update","reduce_edges","reduce_vertices","print")
  val types = List("bool","int","float")

  lexical.reserved ++= keywords
  lexical.reserved ++= types

  lexical.delimiters ++= List("{","}","(",")",".",",","+","-","*","/","<",">","=")

  def numLit:Parser[Int] = opt("-") ~ numericLit ^^ {
    case None ~ n => n.toInt
    case Some(_) ~ n => (-1)*n.toInt
  }



}
