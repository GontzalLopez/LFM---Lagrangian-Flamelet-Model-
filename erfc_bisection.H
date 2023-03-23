double erfc_inv(double x)
{
    const double eps = 1e-8; // tolerance
    double a = -10, b = 10; // initial bracketing interval
    double fa = std::erfc(a) - x; // function value at a
    double fb = std::erfc(b) - x; // function value at b

    while (b - a > eps) {
        double c = (a + b) / 2.0; // midpoint of the interval
        double fc = std::erfc(c) - x; // function value at c

        if (fa * fc > 0) {
            a = c;
            fa = fc;
        } else {
            b = c;
            fb = fc;
        }
    }

    return (a + b) / 2.0;
}
