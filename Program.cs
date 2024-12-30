using System;
using System.Collections.Generic;
using System.IO;

class LotkaVolterra
{
    private double alpha;
    private double beta;
    private double gamma;
    private double delta;

    private double x0;
    private double y0;

    public LotkaVolterra(double alpha, double beta, double gamma, double delta, double x0, double y0)
    {
        this.alpha = alpha;
        this.beta = beta;
        this.gamma = gamma;
        this.delta = delta;
        this.x0 = x0;
        this.y0 = y0;
    }

    private (double, double) SolveNewton(double xPrev, double yPrev, double h)
    {
        double xK = xPrev;
        double yK = yPrev;

        
        for (int i = 0; i < 100; ++i)
        {
            double f1 = xK - xPrev - 0.5 * h * (alpha * xK - beta * xK * yK + alpha * xPrev - beta * xPrev * yPrev);
            double f2 = yK - yPrev - 0.5 * h * (delta * xK * yK - gamma * yK + delta * xPrev * yPrev - gamma * yPrev);

            double df1_dx = 1 - 0.5 * h * (alpha - beta * yK);
            double df1_dy = 0.5 * h * beta * xK;
            double df2_dx = 0.5 * h * delta * yK;
            double df2_dy = 1 - 0.5 * h * (-gamma + delta * xK);

            double det = df1_dx * df2_dy - df1_dy * df2_dx;

            double dx = (f1 * df2_dy - f2 * df1_dy) / det;
            double dy = (df1_dx * f2 - df2_dx * f1) / det;

            xK -= dx;
            yK -= dy;

            if (Math.Abs(dx) < 1e-6 && Math.Abs(dy) < 1e-6)
            {
                break;
            }
        }

        return (xK, yK);
    }

    private (double, double) Dt(double x, double y)
    {
        double dxdt = alpha * x - beta * x * y;
        double dydt = delta * x * y - gamma * y;

        return (dxdt, dydt);
    }

    public void TrapezoidMethod(double h, double tEnd)
    {
        List<(double, double)> results = new List<(double, double)>();

        double t = 0.0;
        double x = x0;
        double y = y0;

        results.Add((t, x));
        results.Add((t, y));

        while (t < tEnd)
        {
            var nextPoint = SolveNewton(x, y, h);

            t += h;
            x = nextPoint.Item1;
            y = nextPoint.Item2;

            results.Add((t, x));
            results.Add((t, y));
        }

        SaveResults("xy_slv.txt", results);
        SaveResults("ty_slv.txt", results, true);
        SaveResults("tx_slv.txt", results, false);
    }

    public void RungeKutta4(double h, double tEnd)
    {
        double t = 0.0;
        double x = x0;
        double y = y0;

        List<(double, double)> results = new List<(double, double)>();

        while (t < tEnd)
        {
            results.Add((t, x));
            results.Add((t, y));

            var k1 = Dt(x, y);
            var k2 = Dt(x + h * k1.Item1 / 2, y + h * k1.Item2 / 2);
            var k3 = Dt(x + h * k2.Item1 / 2, y + h * k2.Item2 / 2);
            var k4 = Dt(x + h * k3.Item1, y + h * k3.Item2);

            x += h / 6 * (k1.Item1 + 2 * k2.Item1 + 2 * k3.Item1 + k4.Item1);
            y += h / 6 * (k1.Item2 + 2 * k2.Item2 + 2 * k3.Item2 + k4.Item2);

            t += h;
        }

        SaveResults("xy_rk4.txt", results);
        SaveResults("ty_rk4.txt", results, true);
        SaveResults("tx_rk4.txt", results, false);
    }

    private void SaveResults(string fileName, List<(double, double)> results, bool isY = false)
    {
        using (var writer = new StreamWriter(fileName))
        {
            foreach (var result in results)
            {
                if (isY)
                    writer.WriteLine($"{result.Item1} {result.Item2}");
                else
                    writer.WriteLine($"{result.Item2} {result.Item1}");
            }
        }
    }
}

class Program
{
    static void Main()
    {
        double alpha = 4.66;
        double beta = 0.01;
        double gamma = 2.354;
        double delta = 0.002;

        double x0 = 2303;
        double y0 = 993;

        double h = 0.1;
        double tEnd = 365; 

        var lotkaVolterra = new LotkaVolterra(alpha, beta, gamma, delta, x0, y0);

        lotkaVolterra.RungeKutta4(h, tEnd);
        lotkaVolterra.TrapezoidMethod(h, tEnd);
    }
}
