import java.lang.Math._
import scala.collection.mutable.ArrayBuffer

/**
  * Created by ypu on 16/03/2017.
  */
object wowa {

  private def owa_gaussian(len: Int, mu: Double, sig: Double): Array[Double] = {
    val alpha = 0.5
    val c = 1.0 / (sig * sqrt(2 * PI))
//    (1 to len).map(x => c * exp(-alpha * pow((x - mu) / sig, 2))).toArray
    // Normalization
    val weights_raw = (1 to len).map(x => c * exp(-alpha * pow((x - mu) / sig, 2))).toArray
    val sum = weights_raw.sum
    val weights = weights_raw.map(_ / sum)
    print("OWA_Gaussian: ")
    weights.indices foreach (i => print(weights(i) + ", "))
    println()
    weights
  }

  /**
    * Compute the OWA of values
    *
    * @param values the values to weight as OWA
    * @param mu     the mean of the gaussians to apply (a.k.a. atLeast)
    * @return the OWA score
    */
  def owa(values: Array[Double], mu: Int): Double = {
    val weights = owa_gaussian(values length, mu, sig = 1.0)
    val valsorted = values sortWith (_ > _)
    (valsorted zip weights).map(Function.tupled(_ * _)) sum // dot product
  }

  /**
    * Score an array of measures with the WOWA algorithm
    * Required : all arrays have the same length
    *
    * @param wa       weights of the WA
    * @param owamu    gaussian mean of the OWA (a.k.a.atLeast)
    * @param measures measures to score
    * @return the WOWA scoring value
    */
  def wowa(wa: Array[Double], owamu: Int, measures: Array[Double]): Double = {
    val owa = owa_gaussian(measures length, owamu, sig = 1.0)
    wowa(wa, owa, measures)
  }

  /**
    * Score an array of measures with the WOWA algorithm
    * Required : all arrays have the same length
    *
    * @param wa       weights of the WA
    * @param owa      weights of the OWA
    * @param measures measures to score
    * @return the WOWA scoring value
    */
  def wowa(wa: Array[Double], owa: Array[Double], measures: Array[Double]): Double = {
    // Compute the WOWA interpolation function (depends only on OWA coefs and is a straight line passing by the origin)
    // points = Set{(i/n, sum(j<=i, owa(i))|i = 1,...,n} and {(0,0)}
    println("OWA weights: ")
    owa.indices foreach(i => print(owa(i) + ", "))
    print(owa.sum + "\n")
    val len = owa.length
    val points = ArrayBuffer[(Double, Double)]((0, 0))
    owa.indices foreach (i => points += (( (i+1).toDouble / len.toDouble, points(i)._2 + owa(i)))) //points: (i/n, partial sum of owa)
//    points += ((1.0, 1.0))
    println("Points: ")
    points.indices foreach(i => print(points(i) + ", "))
    print("\n")
    val slopes = wowa_slopes(points.toArray)
    println("Slopes: ")
    slopes.indices foreach(i => print(slopes(i) + ", "))
    print("\n")
    val knots = insertKnots(slopes, points.toArray)
    println("Knot points: ")
    knots.O.indices foreach(i => print("V: " + knots.V(i) + " O: " + knots.O(i) + " W: " + knots.W(i) + "; "))
    print("\n")
    // Sort the measures and apply the same permutation to the WA coefs
    val sorted = (measures zip wa) sortWith (_._1 > _._1)
    val measures_sorted = sorted map (_._1)
    val wa_sorted = sorted map (_._2)
    // Compute the incremental sum of WA sorted
    val wa_sum = ArrayBuffer[Double](wa_sorted(0))
    (1 until len) foreach (i => wa_sum += wa_sum(i - 1) + wa_sorted(i))
    // Compute the wowa coefs
    val wowa_weights = ArrayBuffer[Double](wowa_interpolation(wa_sum(0), knots, points.toArray))
    (1 until len) foreach (i => wowa_weights += wowa_interpolation(wa_sum(i), knots, points.toArray)
      - wowa_interpolation(wa_sum(i - 1), knots, points.toArray))

    println("WOWA weights: ")
    (0 until len) foreach(i => print(wowa_weights(i) + ", "))
    print("\n")
    println("Sum of WOWA weights: " + wowa_weights.sum)
    // Finally, compute the WOWA score
    (wowa_weights zip measures_sorted).map(Function.tupled(_ * _)) sum
  }

  /**
    * Compute the slope of the wowa coefficients
    * Implement the interpolation function of WOWA following Chen and Otto's interpolation
    * TODO: May have the infinity problems
    * @param points points which are used for interpolation
    * @return the corresponding slope of each point
    */
  private def wowa_slopes(points: Array[(Double, Double)]): Array[Double] = {
    // Compute the slope for each segment between pairs of points
    val len = points.length
    val slopes = ArrayBuffer[Double](0) // slopes(0) is useless
    for(i <- 1 until len) slopes += (points(i)._2 - points(i-1)._2) / (points(i)._1 - points(i-1)._1)
    // Compute the corresponding slope of each point
    var slopesOfPoints = new Array[Double](len)
    for(i <- 1 until len-1) {
      if(slopes(i) * slopes(i+1) <= 0) slopesOfPoints(i) = 0
      else if(slopes(i) == slopes(i+1)) slopesOfPoints(i) = slopes(i)
      else if(math.abs(slopes(i)) > 0 && math.abs(slopes(i+1)) > 0 && slopes(i) * slopes(i+1) > 0) {
        if (math.abs(slopes(i)) > math.abs(slopes(i+1))) {
          /*
            Extend a line through points(i) with slope slopes(i) until it intersects the horizontal line
            through points(i+1) at point b=(b, Y_i+1).
            slopes(i) = (Y_i+1 - Y_i)/(c - X_i), where c = (b + X_i+1)/2, b = (Y_i+1 - Y_i) / S_i + X_i
           */
          val b = (points(i+1)._2 - points(i)._2) / slopes(i) + points(i)._1
          val c = (b + points(i+1)._1) / 2
          slopesOfPoints(i) = (points(i+1)._2 - points(i)._2) / (c - points(i)._1)
        }
        else if(math.abs(slopes(i)) < math.abs(slopes(i+1))) {
          /*
            Extend a line through points(i) with slope slopes(i+1) until it intersects the horizontal line
            through points(i-1) at point b=(b, y_i).
            slopes(i) = (Y_i - Y_i-1)/(X_i - c), where c = (b + X_i-1)/2, b = (Y_i-1 - Y_i) / S_i+1 + X_i
           */
          val b = (points(i-1)._2 - points(i)._2) / slopes(i+1) + points(i)._1
          val c = (b + points(i-1)._1) / 2
          slopesOfPoints(i) = (points(i)._2 - points(i-1)._2) / (points(i)._1 - c)
        }
      }
    }
    slopesOfPoints(0) = {
      if(slopesOfPoints(1) == 0 && slopes(1) == 0) 0
      else pow(slopes(1), 2) / slopesOfPoints(1)
    }
    slopesOfPoints(len-1) = {
      if(slopesOfPoints(len-2) == 0 && slopes(len-1) == 0) 0
      else pow(slopes(len-1), 2) / slopesOfPoints(len - 2)
    }
    slopesOfPoints
  }

  /**
    *
    * @param V
    * @param O Knot points
    * @param W
    */
  private case class Knots(V: Array[(Double, Double)], O: Array[(Double, Double)], W: Array[(Double, Double)])

  /**
    *
    * @param slopes Slopes of all the points to interpolate
    * @param points Points which are used for interpolation
    * @return Knot points between contiguous pairs of points
    */
  private def insertKnots(slopes: Array[Double], points: Array[(Double, Double)]): Knots = {
    /**
      * Li is the line that passes through Di(xi, yi) with slope mi
      * Zi is the intersection point of Li and Li+1, Zi = (ti, zi)
      * Vi = ((xi + ti)/2, Li((xi + ti)/2)), Wi = ((xi+1 + ti)/2, Li+1((xi+1 + ti)/2))
      * L is the line joining Vi and Wi, knot point Oi is defined asfollows:
      * Oi = (ti, L(ti))
      */
    val len = points.length
    val O, V, W = new Array[(Double, Double)](len-1)
    for(i <- 0 until len-1) {
      // Point Z is within the rectangle R determined by the points Di and Di+1
      var t = ((points(i+1)._2 - slopes(i+1)*points(i+1)._1) - (points(i)._2 - slopes(i)*points(i)._1)) / (slopes(i) - slopes(i+1))
      var z = slopes(i)*(t - points(i)._1) + points(i)._2
      // Point Z is outside R or it is Di or Di+1 or Li coincide with Li+1
      if(!(t < points(i+1)._1 && t > points(i)._1 && z < points(i+1)._2 && z > points(i)._2)
        || (slopes(i) == slopes(i+1) && (points(i+1)._2 - points(i)._2)/(points(i+1)._1 - points(i)._1) == slopes(i))) {
        t = (points(i)._1 + points(i+1)._1)/2
        z = slopes(i)*(t - points(i)._1) + points(i)._2
      }
      val v = ((points(i)._1 + t)/2, slopes(i)*(t - points(i)._1)/2 + points(i)._2)
      val w = ((points(i+1)._1 + t)/2, slopes(i+1)*(t - points(i+1)._1)/2 + points(i+1)._2)
      val slope_L = (w._2 - v._2)/(w._1 - v._1) // could be Infinity?
//      println("V: " + v)
//      println("W: " + w)
//      println("Slope_L: " + slope_L)
      val o = (t, slope_L*(t - v._1) + v._2)
      V(i) = v; W(i) = w; O(i) = o
    }
    Knots(V, O, W)
  }

  /**
    *
    * @param x Input of the interpolation function
    * @param knots Knot points which are between contiguous pairs of points
    * @param points Points which are used for interpolation
    * @return Value of interpolation
    */
  private def wowa_interpolation(x: Double, knots: Knots, points: Array[(Double, Double)]): Double = {
    /**
      * Interpolation method based on the second-degree Bernstein polynomial
      * Y(x) = B[pi, vi, oi](x) on [xi, ti] or B[oi, wi, pi+1](x) on [ti, xi+1]
      * B[P1, w, P2](x) = B(g)(x) = (g(x1)*pow(x2-x, 2) + 2b(x-x1)(x2-x) + g(x2)*pow(x-x1, 2)) / pow(x2-x1, 2)
      */
    val i = points.sliding(2).indexWhere{case Array(a, b) => a._1 <= x && b._1 >= x}
    i match {
      case -1 => 0
      case _ =>
        val t = knots.O(i)._1
        val bernstein: ((Double, Double), (Double, Double), (Double, Double), Double) => Double = (P, O, D, x) =>
          (P._2*pow(D._1 - x, 2) + 2*O._2*(x - P._1)*(D._1 - x) + D._2*pow(x - P._1, 2)) / pow(D._1 - P._1, 2)
        if(x >= points(i)._1 && x <= t) bernstein(points(i), knots.V(i), knots.O(i), x)
        else bernstein(knots.O(i), knots.W(i), points(i+1), x)
    }
  }


  // The code below is for WOWA testing
  case class Session(seq: Double, OF: Double, US: Double, DC: Double, SP:Double, FA:Double,
                OS: Double, DV: Double, AG:Double, IP:Double, TS:Double, DI:Double, WD:Double) {
  }

  private def wowa_test(session: Session): Unit = {
    val wa_raw = Array[Double](2, 1.5, 1, 0.8, 0.7, 0.5, 0.3, 0.2)
    val sum = wa_raw.sum
    val wa = wa_raw.map(_ / sum)

    val s1 = owa(Array(session.OF, session.US, session.DC), 1) // office, user-sign, duty_code
    val s2 = owa(Array(session.SP, session.FA), 1) // sap, farm
    val s3 = owa(Array(session.OS, session.DV, session.AG), 1) // os, device, agent
    val seq = session.seq //last action's score
    println(s"OWA scores are: $s1, $s2, $s3 ")
    // Array(seq, ip, timestamp, ip_dist, week_day)
//    println("Result of OWA is(compared to WOWA): " + owa(mu = 2, values = Array(seq, session.IP, session.TS, session.DI, session.WD, s1, s2, s3)) + "\n")
    println("Result of WOWA is: " + wowa(wa = wa, owamu = 2, measures = Array(seq, session.IP, session.TS, session.DI, session.WD, s1, s2, s3)) + "\n")
  }

  def main(args: Array[String]): Unit = {

    //score of {Seq, office, user_sign, duty_code, sap, farm, os, device, agent, ip, timestamp, ip_dist, week_day}
    val session1 = Session(seq = 0.982,
      OF = 0, US = 0, DC = 0,
      SP = 1, FA = 0,
      OS = 0, DV = 0, AG = 1,
      IP = 1, TS = 0.799, DI = 0, WD = 0.819)

    val session2 = Session(seq = 0.984,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0.026,
      IP = 0, TS = 0.44, DI = 0, WD = 0.657)

    val session3 = Session(seq = 0.824,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 0.07, TS = 0.108, DI = 0, WD = 0.696)

    val session4 = Session(seq = 0.753,
      OF = 0, US = 0, DC = 0,
      SP = 1, FA = 0,
      OS = 0, DV = 0, AG = 1,
      IP = 1, TS = 0.799, DI = 0, WD = 0.819)

    val session5 = Session(seq = 0.853,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 1, TS = 0.454, DI = 0, WD = 1)

    val session6 = Session(seq = 0.408,
      OF = 0, US = 0, DC = 0,
      SP = 1, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 0, TS = 0.475, DI = 0, WD = 1)

    val session7 = Session(seq = 0.146,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 1,
      IP = 0, TS = 0.319, DI = 0, WD = 0.696)

    val session8 = Session(seq = 1,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 1, TS = 0.454, DI = 0, WD = 1)

    val session9 = Session(seq = 0.696,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 1, TS = 0.594, DI = 0, WD = 1)

    val session10 = Session(seq = 0,
      OF = 0, US = 0, DC = 0,
      SP = 0, FA = 0,
      OS = 0, DV = 0, AG = 0,
      IP = 0, TS = 0.006, DI = 0, WD = 1)


    println("\n Test 1 \n")
    wowa_test(session1)

    println("\n Test 2 \n")
    wowa_test(session2)

    println("\n Test 3 \n")
    wowa_test(session3)

    println("\n Test 4 \n")
    wowa_test(session4)

    println("\n Test 5 \n")
    wowa_test(session5)

    println("\n Test 6 \n")
    wowa_test(session6)

    println("\n Test 7 \n")
    wowa_test(session7)

    println("\n Test 8 \n")
    wowa_test(session8)

    println("\n Test 9 \n")
    wowa_test(session9)

    println("\n Test 10 \n")
    wowa_test(session10)
  }
}
