#ifndef _THRESHOLD_H
#define _THRESHOLD_H

template<class X>
bool IsBetween(X x, X a, X b) {
  return ((a <= x) && (x < b)) || ((b < x) && (x <= a));
}


template<class Number>
Number Threshold2(Number density,
                  Number threshold_01_a,
                  Number threshold_01_b)
{
  // Then the thresholding function is either:
  //
  //           if thresh_01_a < thresh_01_b
  //
  // output
  // density
  //  /|\                              _________________\
  // 1 |                           _.-'                 /
  //   |                       _,-'                 
  //   |                   _,-'            
  // 0 |________________,-'                     ________\ input
  //                 thresh         thresh              / density
  //                  01_a           01_b 
  //
  //
  //     OR    if thresh_01_b < thresh_01_a
  //
  //
  // output
  // density
  //
  //  /|\
  //   |________________
  // 1 |                `-._
  //   |                    `-._ 
  //   |                        '-._            
  // 0 |                            `-.___________________\ input
  //                 thresh         thresh                / density
  //                  01_b           01_a 
        

  Number g; // "g" is the new intensity after threshold filter.

  if (IsBetween(density,
                threshold_01_a,
                threshold_01_b))
    g =
      (density-threshold_01_a) /
      (threshold_01_b-threshold_01_a);
  else if ((density-threshold_01_a)
           *
           (threshold_01_b-threshold_01_a)
           > 0.0)
    g = 1.0;
  else
    g = 0.0;

  return g;
} // Threshold2()





template<class Number>
Number Threshold4(Number density,
                  Number threshold_01_a,
                  Number threshold_01_b,
                  Number threshold_10_a,
                  Number threshold_10_b)
{
  // In that case, the thresholding function is:
  //
  // output
  // density
  //
  //  /|\
  // 1 |                 ________________                
  //   |             _,-'                `-._
  //   |         _,-'                        `-._
  // 0 |______,-'                                `-._______\ input
  //        thresh    thresh          thresh    thresh     / density
  //         01_a      01_b            10_a      10_b
  //
  //  Or, if the user selects thresh_10_b < thresh_01_a, then use:
  //
  //  /|\                                                   
  // 1 |_____                                       _______\ input
  //   |     `-._                               _.-'       / density
  //   |         `-._                       _,-'
  // 0 |             `-._________________,-'
  //        thresh    thresh          thresh    thresh
  //         10_a      10_b            01_a      01_b
  //
  //
  
  Number g;
  g = Threshold2(density,
                 threshold_01_a,
                 threshold_01_b);

  if ((threshold_01_b == threshold_10_a) &&
      (threshold_01_b == threshold_10_b))
    return g; // Default: if 2 extra arguments don't make sense, ignore them

  if (IsBetween(density,
                threshold_01_a,
                threshold_01_b))
    g =
      (density-threshold_01_a) /
      (threshold_01_b-threshold_01_a);
  else if (IsBetween(density,
                     threshold_10_a,
                     threshold_10_b))
    g =
      (density-threshold_10_b) /
      (threshold_10_a-threshold_10_b);
  else if (IsBetween(density,
                     threshold_01_b,
                     threshold_10_a))
    g = 1.0;
  else if (IsBetween(density,
                     threshold_10_b,
                     threshold_01_a))
    g = 0.0;
  else if (threshold_01_b <
           threshold_10_a) //upper picture
    g = 0.0;
  else if (threshold_10_b <
           threshold_01_a) //lower picture
    g = 1.0;
  else
    throw "Invalid threshold settings (order violation).\n";

  return g;

} // Threshold4()




#endif //#ifndef _THRESHOLD_H
