"use strict";

const EPS = 3e-16;
const BIG = 1.0 / EPS;
const PI = 3.141592653589793238462643383279502884197;
const PID2 = 1.570796326794896619231321691639751442099;
const RPI = 0.3183098861837906715377675267450287240689;
const RPISQ = RPI * RPI;
const LARGEF = RPI / EPS;
const LARGEG = (RPI * LARGEF) ** (1.0 / 3.0);
const LARGEX = 1.0/Math.sqrt(Math.sqrt(2.225E-308)); // DMACH(1) constant

function DCOSPX (X) {
  /*
  > 1996-01-29 DCOSPX WVS JPL Add better acknowledgement for origins.
  > 1994-10-20 DCOSPX Krogh  Changes to use M77CON
  > 1994-04-22 DCOSPX CLL Made SP and DP codes similar.
  > 1993-05-06 DCOSPX WVS JPL Convert from NSWC to Math 77
  ----------------------------------------------------------------------
  This procedure was originally procedure COS1 from the Naval Surface
  Warfare Center library.
  ----------------------------------------------------------------------
 
                          EVALUATION OF COS(PI*X)
 
                              --------------
 
      THE EXPANSION FOR SIN(PI*A) (ABS(A) .LE. 1/4) USING A1,...,A13
      IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
      THE EXPANSION FOR COS(PI*A) (ABS(A) .LE. 1/4) USING B1,...,B13
      IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.
 
      The polynomials using coefficients SA0, ... SA6 and SB1, ..., SB6
      give approximations whose largest observed relative error in the
      relevant intervals is 0.129e-14.
      We will use this latter approximation when the machine epsilon
      is larger than 0.2d-14.
  ----------------------------------------------------------------------
  */

  const CUTOFF = 0.2e-14;

  const SA0 = .314159265358979e1, SA1 = -.516771278004995e1,
        SA2 = .255016403987327e1, SA3 = -.599264528932149,
        SA4 = .821458689493251e-1, SA5 = -.737001831310553e-2,
        SA6 = .461514425296398e-3,
        SB1 = -.493480220054460e1, SB2 = .405871212639605e1,
        SB3 = -.133526276691575e1, SB4 = .235330543508553,
        SB5 = -.258048861575714e-1, SB6 = .190653140279462e-2,
        A1 = -.1028083791780141522795259479153765743002,
        A2 =  .3170868848763100170457042079710451905600e-2,
        A3 = -.4657026956105571623449026167864697920000e-4,
        A4 =  .3989844942879455643410226655783424000000e-6,
        A5 = -.2237397227721999776371894030796800000000e-8,
        A6 =  .8847045483056962709715066675200000000000e-11,
        A7 = -.2598715447506450292885585920000000000000e-13,
        A8 =  .5893449774331011070033920000000000000000e-16,
        A9 = -.1062975472045522550784000000000000000000e-18,
        A10 =  .1561182648301780992000000000000000000000e-21,
        A11 = -.1903193516670976000000000000000000000000e-24,
        A12 =  .1956617650176000000000000000000000000000e-27,
        A13 = -.1711276032000000000000000000000000000000e-30,
        B1 = -.3084251375340424568385778437461297229882,
        B2 =  .1585434424381550085228521039855226435920e-1,
        B3 = -.3259918869273900136414318317506279360000e-3,
        B4 =  .3590860448591510079069203991239232000000e-5,
        B5 = -.2461136950494199754009084061808640000000e-7,
        B6 =  .1150115912797405152263195572224000000000e-9,
        B7 = -.3898073171259675439899172864000000000000e-12,
        B8 =  .1001886461636271969091584000000000000000e-14,
        B9 = -.2019653396886572027084800000000000000000e-17,
        B10 =  .3278483561466560512000000000000000000000e-20,
        B11 = -.4377345082051788800000000000000000000000e-23,
        B12 =  .4891532381388800000000000000000000000000e-26,
        B13 = -.4617089843200000000000000000000000000000e-29;

  let A = Math.abs(X);
  if (A >= BIG) throw 'ABS(X) is too large';

  let N = Math.floor(A);
  A = A - N;
  let DCOSPX;
  if (0.25 <= A && A <= 0.75) {
    A = 0.25 + (0.25 - A);
    if(EPS < CUTOFF) {
       const T = 16.0*A*A;
       const W = (((((((((((((A13*T + A12)*T + A11)*T + A10)*T + A9)*T +
                A8)*T + A7)*T + A6)*T + A5)*T + A4)*T + A3)*T +
                A2)*T + A1)*T + 0.5) + 0.5;
       DCOSPX = PI*A*W;
    } else {
       T = A*A;
       DCOSPX = ((((((SA6*T + SA5)*T + SA4)*T + SA3)*T + SA2)*T
                     + SA1)*T + SA0)*A;
    }
  } else {
    // A < 0.25 or A > 0.75
    if (A > 0.75) {
      A = 0.25 + (0.75 - A);
      N = N - 1;
    }

    if (EPS < CUTOFF) {
       const T = 16.0*A*A;
       DCOSPX = (((((((((((((B13*T + B12)*T + B11)*T + B10)*T + B9)*T+
                    B8)*T + B7)*T + B6)*T + B5)*T + B4)*T + B3)*T +
                    B2)*T + B1)*T + 0.5) + 0.5;
    } else {
       const T = A*A;
       DCOSPX = ((((((SB6*T + SB5)*T + SB4)*T + SB3)*T + SB2)*T
                      + SB1)*T + 0.5) + 0.5;
    }
  }

  if ((N % 2) != 0) DCOSPX = -DCOSPX;
  return DCOSPX;
}

function DCSPXX(X) {
  /*
  1996-01-08 DCSPXX WV Snyder Original code
  COS(PI * X * X / 2) carefully to avoid loss of precision for large X
  DCOSPX is used to compute COS(PI * X)

  BIG = 1 / round-off = biggest integer exactly representable by F.P.
     If X > BIG then to the working precision x**2 is an integer (which
     we assume to be a multiple of four), so cos(pi/2 * x**2) = 1.
  N = [X], and later [K F]
  F = X - N = fractional part of X
  K = [ N / 2 ]
  J = N mod 2
  M = [K F]
  G = K F - M = fractional part of K F
  */

  let f = Math.abs(X);
  if (f > BIG) {
    // Assume x is an even integer.
    return 1.0;
  }

  let n = Math.floor(f);
  f = f - n;
  let k = Math.floor(n / 2);
  let j = n % 2;
  let g = k * f;
  n = Math.floor(g);
  g = g - n;
  if (j == 0)
    return DCOSPX(0.5*f*f + g + g);
  else
    return -DSINPX(0.5*f*f + f + g + g);
}

function DSINPX (X) {
  /*
  > 1996-01-29 SCOSPX WVS JPL Add better acknowledgement for origins.
  > 1994-10-20 DSINPX Krogh  Changes to use M77CON
  > 1994-04-22 DSINPX CLL Made SP and DP codes similar.
  > 1993-05-06 DSINPX WVS JPL Convert from NSWC to Math 77

  ----------------------------------------------------------------------
  This procedure was originally procedure SIN1 from the Naval Surface
  Warfare Center library.
  ----------------------------------------------------------------------
 
                         EVALUATION OF SIN(PI*X)
 
                              --------------
 
      THE EXPANSION FOR SIN(PI*A) (ABS(A) .LE. 1/4) USING A1,...,A13
      IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
      THE EXPANSION FOR COS(PI*A) (ABS(A) .LE. 1/4) USING B1,...,B13
      IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.
 
      The polynomials using coefficients SA0, ... SA6 and SB1, ..., SB6
      give approximations whose largest observed relative error in the
      relevant intervals is 0.129d-14.
      We will use this latter approximation when the machine epsilon
      is larger than 0.2d-14.
 -----------------------------------------------------------------------
 */
  const CUTOFF = 0.2e-14;

  const SA0 = .314159265358979e1, SA1 = -.516771278004995e1,
        SA2 = .255016403987327e1, SA3 = -.599264528932149e0,
        SA4 = .821458689493251e-1, SA5 = -.737001831310553e-2,
        SA6 = .461514425296398e-3,
        SB1 = -.493480220054460e1, SB2 = .405871212639605e1,
        SB3 = -.133526276691575e1, SB4 = .235330543508553e0,
        SB5 = -.258048861575714e-1, SB6 = .190653140279462e-2,
        A1 = -.1028083791780141522795259479153765743002E0,
        A2 =  .3170868848763100170457042079710451905600E-2,
        A3 = -.4657026956105571623449026167864697920000E-4,
        A4 =  .3989844942879455643410226655783424000000E-6,
        A5 = -.2237397227721999776371894030796800000000E-8,
        A6 =  .8847045483056962709715066675200000000000E-11,
        A7 = -.2598715447506450292885585920000000000000E-13,
        A8 =  .5893449774331011070033920000000000000000E-16,
        A9 = -.1062975472045522550784000000000000000000E-18,
        A10 =  .1561182648301780992000000000000000000000e-21,
        A11 = -.1903193516670976000000000000000000000000e-24,
        A12 =  .1956617650176000000000000000000000000000e-27,
        A13 = -.1711276032000000000000000000000000000000e-30,
        B1 = -.3084251375340424568385778437461297229882E0,
        B2 =  .1585434424381550085228521039855226435920E-1,
        B3 = -.3259918869273900136414318317506279360000E-3,
        B4 =  .3590860448591510079069203991239232000000E-5,
        B5 = -.2461136950494199754009084061808640000000E-7,
        B6 =  .1150115912797405152263195572224000000000E-9,
        B7 = -.3898073171259675439899172864000000000000E-12,
        B8 =  .1001886461636271969091584000000000000000E-14,
        B9 = -.2019653396886572027084800000000000000000E-17,
        B10 =  .3278483561466560512000000000000000000000e-20,
        B11 = -.4377345082051788800000000000000000000000e-23,
        B12 =  .4891532381388800000000000000000000000000e-26,
        B13 = -.4617089843200000000000000000000000000000e-29;

  let A = Math.abs(X);
  if (A >= BIG) throw "ABS(X) is too large";

  let N = Math.floor(A);
  A = A - N;
  let DSINPX;
  if (0.25 <= A && A <= 0.75) {
    A = 0.25 + (0.25 - A);
    if (EPS < CUTOFF) {
      let T = 16.0*A*A;
      DSINPX = (((((((((((((B13*T + B12)*T + B11)*T + B10)*T + B9)*T+
                  B8)*T + B7)*T + B6)*T + B5)*T + B4)*T + B3)*T +
                  B2)*T + B1)*T + 0.5) + 0.5;
    } else {
      let T = A*A;
      DSINPX = ((((((SB6*T + SB5)*T + SB4)*T + SB3)*T + SB2)*T
                    + SB1)*T + 0.5) + 0.5;
    }
  } else {
    // A < 0.25 or A > 0.75
    if (A > 0.75) {
      A = 0.25 + (0.75 - A);
    }
   
    if(EPS < CUTOFF) {
      let T = 16.0*A*A;
      let W = (((((((((((((A13*T + A12)*T + A11)*T + A10)*T + A9)*T +
                A8)*T + A7)*T + A6)*T + A5)*T + A4)*T + A3)*T +
                A2)*T + A1)*T + 0.5) + 0.5;
      DSINPX = PI*A*W;
    } else {
      let T = A*A;
      DSINPX = ((((((SA6*T + SA5)*T + SA4)*T + SA3)*T + SA2)*T
                      + SA1)*T + SA0)*A;
    }
  }

  if (X < 0.0) DSINPX = -DSINPX;
  if ((N % 2) != 0) DSINPX = -DSINPX;
    
  return DSINPX;
}

function DSNPXX(X) {
  /*
  1996-01-08 DSNPXX WV Snyder Original code
  SIN(PI * X * X / 2) carefully to avoid loss of precision for large X
  DSINPX is used to compute SIN(PI * X)
  BIGX = 1 / round-off = biggest integer exactly representable by F.P.
     If X > BIGX then to the working precision x**2 is an integer (which
     we assume to be a multiple of four), so sin(pi/2 * x**2) = 0.
  N = [X], and later [K F]
  F = X - N = fractional part of X
  K = [ N / 2 ]
  J = N mod 2
  G = K F - M = fractional part of K F
  */

  let f = Math.abs(X);
  if (f > BIG) {
    // Assume x is an even integer.
    return 0.0;
  }
  
  let n = Math.floor(f);
  f = f - n;
  let k = Math.floor(n / 2);
  let j = n % 2;
  let g = k * f;
  n = Math.floor(g);
  g = g - n;
  if (j == 0)
    return DSINPX(0.5*f*f + g + g)
  else
    return DCOSPX(0.5*f*f + f + g + g)
}

function FRESNEL(x) {
/*
      .  Copyright (C) 1992, California Institute of Technology.
      .  U. S. Government sponsorship under
      .  NASA contract NAS7-918 is acknowledged.
  > 1996-01-08 DFRENL WV Snyder Use DCSPXX for cos(Pi/2 x**2), etc.
  > 1995-11-03 DFRENL Krogh  Removed blanks in numbers for C conversion.
  > 1994-11-02 DFRENL Krogh  Changes to use M77CON
  > 1994-10-18 DFRENL WV Snyder More specializing instructions
  > 1993-02-25 DFRENL CLL. Edited to eliminate ENTRY and EQUIVALENCE.
  > 1992-09-15 DFRENL WV Snyder Specializing instructions
  > 1992-04-13 DFRENL WV Snyder Declare DFRENF, DFRENG, DFRENS
  > 1992-03-18 DFRENL WV Snyder Move declarations for coefficient arrays
  > 1992-01-24 DFRENL WV Snyder Original code

  Compute the Fresnel Cosine and Sine integrals C(x) and S(x) and the
  auxiliary functions f(x) and g(x), for any X.
 
  Developed by W. V. Snyder, Jet Propulsion Laboratory, 24 January 1992.
 
  Ref: W. J. Cody, "Chebyshev Approximations for the Fresnel Integrals",
  Mathematics of Computation, 1968, pp 450-453 plus Microfiche Suppl.
  W. V. Snyder, "Algorithm 723: Fresnel Integrals," ACM Trans. Math.
  Softw. 19, 4 (December 1993) 452-456.
  Accuracies of highest order formulae, where E is relative error:
 
  Range           Function   -log10(E)   Function   -log10(E)
  |X|<=1.2          C(x)       16.24       S(x)       17.26
  1.2<|X|<=1.6      C(x)       17.47       S(x)       18.66
  1.6<|X|<=1.9      f(x)       17.13       g(x)       16.25
  1.9<|X|<=2.4      f(x)       16.64       g(x)       15.65
  2.4<|X|           f(x)       16.89       g(x)       15.58
 
  Refer to Cody for accuracy of other approximations.


                         Internal variables.
 
  PID2 is pi / 2.
  RPI is the reciprocal of PI.
  RPISQ is the reciprocal of PI squared.
  AX is abs(x).
  BIGX is 1/sqrt(round-off).  If X > BIGX then to the working
          precision x**2 is an integer (which we assume to be a multiple
          of four), so cos(pi/2 * x**2) = 1, and sin(pi/2 * x**2) = 0.
  C and S are values of C(x) and S(x), respectively.
  CX and SX are cos(pi/2 * ax**2) and sin(pi/2 * ax**2), respectively.
  F and G are used to compute f(x) and g(x) when X > 1.6.
  LARGEF is 1/(pi * underflow).  If X > LARGEF then f ~ 0.
  LARGEG is cbrt(1/(pi**2 * underflow)).  If X > LARGEG then g ~ 0.
  LARGEX is 1/sqrt(sqrt(underflow)).  If X > LARGEX then f ~ 1/(pi * x)
          and g ~ 1/(pi**2 * x**3).
  MODE indicates the function to be computed: 1 = C(x), 2 = S(x),
          3 = f(x), 4 = g(x).
  X4 is either X ** 4 or (1.0/X) ** 4.
      If you change the order of approximation, you must change the
      declarations and DATA statements for the coefficient arrays,
      and the executable statements that evaluate the approximations.
  */

  const pc1 = [9.999999999999999421e-1, -1.994608988261842706e-1, 1.761939525434914045e-2,
               -5.280796513726226960e-4, 5.477113856826871660e-6];
  const qc1 = [null, 4.727921120104532689e-2, 1.099572150256418851e-3, 1.552378852769941331e-5, 1.189389014228757184e-7];
  const pc2 = [1.00000000000111043640e0, -2.07073360335323894245e-1, 1.91870279431746926505e-2,
               -6.71376034694922109230e-4, 1.02365435056105864908e-5, -5.68293310121870728343e-8];
  const qc2 = [null, 3.96667496952323433510e-2, 7.88905245052359907842e-4, 1.01344630866749406081e-5,
               8.77945377892369265356e-8, 4.41701374065009620393e-10];
  const ps1 = [5.2359877559829887021e-1, -7.0748991514452302596e-2, 3.8778212346368287939e-3,
               -8.4555728435277680591e-5, 6.7174846662514086196e-7];
  const qs1 = [null, 4.1122315114238422205e-2, 8.1709194215213447204e-4, 9.6269087593903403370e-6, 5.9528122767840998345e-8];
  const ps2 = [5.23598775598344165913e-1, -7.37766914010191323867e-2, 4.30730526504366510217e-3,
               -1.09540023911434994566e-4, 1.28531043742724820610e-6, -5.76765815593088804567e-9];
  const qs2 = [null, 3.53398342767472162540e-2, 6.18224620195473216538e-4, 6.87086265718620117905e-6,
               5.03090581246612375866e-8, 2.05539124458579596075e-10];
  const pf1 = [3.1830975293580985290e-1, 1.2226000551672961219e1, 1.2924886131901657025e2,
               4.3886367156695547655e2, 4.1466722177958961672e2, 5.6771463664185116454e1];
  const qf1 = [null, 3.8713003365583442831e1, 4.1674359830705629745e2, 1.4740030733966610568e3,
               1.5371675584895759916e3, 2.9113088788847831515e2];
  const pf2 = [3.183098818220169217e-1, 1.958839410219691002e1, 3.398371349269842400e2,
               1.930076407867157531e3, 3.091451615744296552e3, 7.177032493651399590e2];
  const qf2 = [null, 6.184271381728873709e1, 1.085350675006501251e3, 6.337471558511437898e3,
               1.093342489888087888e4, 3.361216991805511494e3];
  const pf3 = [-9.675460329952532343e-2, -2.431275407194161683e1, -1.947621998306889176e3,
               -6.059852197160773639e4, -7.076806952837779823e5, -2.417656749061154155e6, -7.834914590078317336e5];
  const qf3 = [null, 2.548289012949732752e2, 2.099761536857815105e4, 6.924122509827708985e5,
               9.178823229918143780e6, 4.292733255630186679e7, 4.803294784260528342e7];
  const pg1 = [1.013206188102747985e-1, 4.445338275505123778e0, 5.311228134809894481e1,
               1.991828186789025318e2, 1.962320379716626191e2, 2.054214324985006303e1];
  const qg1 = [null, 4.539250196736893605e1, 5.835905757164290666e2, 2.544731331818221034e3,
               3.481121478565452837e3, 1.013794833960028555e3];
  const pg2 = [1.01321161761804586e-1, 7.11205001789782823e0, 1.40959617911315524e2,
               9.08311749529593938e2, 1.59268006085353864e3, 3.13330163068755950e2];
  const qg2 = [null, 7.17128596939302198e1, 1.49051922797329229e3, 1.06729678030580897e4,
               2.41315567213369742e4, 1.15149832376260604e4];
  const pg3 = [-1.53989733819769316e-1, -4.31710157823357568e1, -3.87754141746378493e3, -1.35678867813756347e5,
               -1.77758950838029676e6, -6.66907061668636416e6, -1.72590224654836845e6];
  const qg3 = [null, 2.86733194975899483e2, 2.69183180396242536e4, 1.02878693056687506e6,
               1.62095600500231646e7, 9.38695862531635179e7, 1.40622441123580005e8];

  let ax = Math.abs(x);
  let c, s, f, g;
  if (ax <= 1.6) {
    let x4 = ax**4;

    if (ax <= 1.2)
       c = x * ((((pc1[4]*x4+pc1[3])*x4+pc1[2])*x4+pc1[1])*x4+
                      pc1[0])
            / ((((qc1[4]*x4+qc1[3])*x4+qc1[2])*x4+qc1[1])*x4+1.0);
    else
       c = x * (((((pc2[5]*x4+pc2[4])*x4+pc2[3])*x4+pc2[2])*x4+
                      pc2[1])*x4+pc2[0])
                / (((((qc2[5]*x4+qc2[4])*x4+qc2[3])*x4+qc2[2])*x4+
                    qc2[1])*x4+1.0);

    if (ax <= 1.2)
       s = x**3*((((ps1[4]*x4+ps1[3])*x4+ps1[2])*x4+ps1[1])*x4+
                       ps1[0])
            / ((((qs1[4]*x4+qs1[3])*x4+qs1[2])*x4+qs1[1])*x4+1.0);
    else
       s = x**3*(((((ps2[5]*x4+ps2[4])*x4+ps2[3])*x4+ps2[2])*x4+
                       ps2[1])*x4+ps2[0])
               /   (((((qs2[5]*x4+qs2[4])*x4+qs2[3])*x4+qs2[2])*x4+
                       qs2[1])*x4+1.0);

    let cx = DCSPXX(ax);
    let sx = DSNPXX(ax);
    f = (0.5 - s) * cx - (0.5 - c) * sx;
    g = (0.5 - c) * cx + (0.50 - s) * sx;
  } else {
    if (ax <= LARGEX) {
      let x4 = (1.0 / ax) ** 4;

      if (ax <= 1.9)
         f = (((((pf1[5]*x4+pf1[4])*x4+pf1[3])*x4+pf1[2])*x4+
                 pf1[1])*x4+pf1[0]) /
          ((((((qf1[5]*x4+qf1[4])*x4+qf1[3])*x4+qf1[2])*x4+
                 qf1[1])*x4+1.0) * ax);
      else if (ax <= 2.4)
         f = (((((pf2[5]*x4+pf2[4])*x4+pf2[3])*x4+pf2[2])*x4+
                 pf2[1])*x4+pf2[0]) /
          ((((((qf2[5]*x4+qf2[4])*x4+qf2[3])*x4+qf2[2])*x4+
                 qf2[1])*x4+1.0) * ax);
      else
         f = (RPI +
           x4*((((((pf3[6]*x4+pf3[5])*x4+pf3[4])*x4+pf3[3])*x4+
                pf3[2])*x4+pf3[1])*x4+pf3[0]) /
              ((((((qf3[6]*x4+qf3[5])*x4+qf3[4])*x4+qf3[3])*x4+
                qf3[2])*x4+qf3[1])*x4+1.0)) / ax;

      if (ax <= 1.9)
         g = (((((pg1[5]*x4+pg1[4])*x4+pg1[3])*x4+pg1[2])*x4+
                    pg1[1])*x4+pg1[0]) /
             ((((((qg1[5]*x4+qg1[4])*x4+qg1[3])*x4+qg1[2])*x4+
                    qg1[1])*x4+1.0) * ax**3);
      else if (ax <= 2.4)
         g = (((((pg2[5]*x4+pg2[4])*x4+pg2[3])*x4+pg2[2])*x4+
                     pg2[1])*x4+pg2[0]) /
             ((((((qg2[5]*x4+qg2[4])*x4+qg2[3])*x4+qg2[2])*x4+
                    qg2[1])*x4+1.0) * ax**3);
      else
         g = (RPISQ +
              x4*((((((pg3[6]*x4+pg3[5])*x4+pg3[4])*x4+pg3[3])*x4+
                   pg3[2])*x4+pg3[1])*x4+pg3[0]) /
                ((((((qg3[6]*x4+qg3[5])*x4+qg3[4])*x4+qg3[3])*x4+
                   qg3[2])*x4+qg3[1])*x4+1.0)) / ax**3;
    } else {
      if (x <= LARGEF)
        f = RPI / ax;
      else
        f = 0.0;

      if (x <= LARGEG)
        g = RPISQ / ax**3;
      else
        g = 0.0;
    }

    let cx = DCSPXX(ax);
    let sx = DSNPXX(ax);
    c = 0.5 + f*sx - g*cx;
    s = 0.5 - f*cx - g*sx;
    if (x < 0.0) {
      c = -c
      s = -s
    }

    if (x < 0.0) {
      // We COULD do the following before the preceeding, and then
      // not put in a test in the preceeding for x .lt. 0, but
      // even though the results are mathematically identical, we
      // would have some cancellation above if we did so.
      g = cx + sx - g;
      f = cx - sx - f;
    }
  }

  return {
    x: c,
    y: s,
    f,
    g
  };
}

module.exports = FRESNEL;

