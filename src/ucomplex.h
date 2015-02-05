typedef struct COMPLEX {double r,i;} complex;

/* complex functions */
complex Cadd(complex z1, complex z2);  /* addition */
complex Csub(complex z1, complex z2);  /* substraction */
complex Cmul(complex z1, complex z2);  /* multiplication */
complex RCmul(double r, complex z);    /* double*complex */
complex Cdiv(complex z1, complex z2);  /* division */
complex Complex(double r, double i);   /* convert to double */

complex Conjg (complex z);    /* conjuge */

/* inverse function 1/z */
complex Cinv (complex z);

/* complex functions with double return values */
double Cabs (complex z);        /* module */
double Carc (complex z);        /* argument : z1 / z = p.e^ia */

/* elementary functions */
complex Cexp (complex z);    /* exponential */
complex Clog (complex z);    /* natural logarithm */
complex Csqrt (complex z);   /* square root */

/* complex trigonometric functions  */
complex Ccos (complex z);    /* cosinus */
complex Csin (complex z);    /* sinus */
complex Ctan (complex z);    /* tangent */

/* inverse complex trigonometric functions */
complex Carc_cos (complex z);   /* arc cosinus */
complex Carc_sin (complex z);   /* arc sinus */
complex Carc_tan (complex z);   /* arc tangent */

/* hyperbolic complex functions */
complex Cch (complex z);     /* hyperbolic cosinus */
complex Csh (complex z);     /* hyperbolic sinus */
complex Cth (complex z);     /* hyperbolic tangent */

/* inverse hyperbolic complex functions */
complex Carc_ch (complex z); /* hyperbolic arc cosinus */
complex Carc_sh (complex z); /* hyperbolic arc sinus */
complex Carc_th (complex z); /* hyperbolic arc tangente */

