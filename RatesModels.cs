using System;
using System.Collections.Generic;

namespace Statistics
{
    public class GaussianFunctions
    {
        static public double spdf(double x)
        {
            double p = Math.Exp(-(x * x) / 2) / (Math.Sqrt(2 * Math.PI));

            return p;
        }

        static public double scdf(double x)
        {
            double x_min = -5.0;
            double x_max = 5.0;
            double intervals = 10000;
            double increase = (x_max - x_min) / intervals;

            if (x <= x_min)
            {
                return 0;
            }

            if (x >= x_max)
            {
                return 1;
            }

            else
            {
                double sum = 0;
                double counter = x_min;
                double var = 0;
                while (counter < x)
                {
                    var = counter + increase / 2;
                    sum = sum + GaussianFunctions.spdf(var) * (increase);
                    counter = counter + increase;
                }
                return sum;
            }
        }
    }

    public class NonCentralChiSquared
    {
        public static double NCX_cdf(double x, double l, double k) // l = nc, k = df
        {
            double h = 1 - 2 / 3 * (k + l) * (k + 3 * l) / Math.Pow(k + 2 * l, 2);
            double p = (k + 2 * l) / Math.Pow(k + l, 2);
            double m = (h - 1) * (1 - 3 * h);

            double denominator = h * Math.Sqrt(2 * p) * (1 + 0.5 * m * p);
            double nom1term = Math.Pow(x / (k + l), h);
            double nom2term = 1 + h * p * (h - 1 - 0.5 * (2 - h) * m * p);
            double nominator = nom1term - nom2term;

            double value = nominator / denominator;

            return GaussianFunctions.scdf(value);
        }
    }
}

namespace RatesModels
{
    using Statistics;

    public interface IBondModel
    {
        double P(double t, double s);
        double R(double t, double s);
        double YieldVolatility(double t, double s);
    }

    public abstract class BondModel : IBondModel
    {
        // SDE Data
        public double kappa; // Speed of mean-reversion
        public double theta; // Long term mean leversion level
        public double sigma; // Volatility of the short rate
        public double r0;  // Current Rate

        // Inherited from interface
        public abstract double P(double t, double s);
        public abstract double R(double t, double s);
        public abstract double YieldVolatility(double t, double s);

        //Main Function
        public BondModel(double kappa, double theta, double vol, double r)
        {
            this.kappa = kappa;
            this.theta = theta;
            this.sigma = vol;
            this.r0 = r;
        }

        //Visitor Pattern
        public abstract void Accept(BondVisitor visitor);
    }

    // Class for Vasicek model
    public class VasicekModel : BondModel
    {
        private double LTRate;

        public VasicekModel(double kappa, double theta, double vol, double r) : base(kappa, theta, vol, r)
        {
            LTRate = theta - 0.5 * vol * vol / (kappa * kappa);
        }

        private double A(double t, double s)
        {
            double R = LTRate;
            double smt = s - t;
            double exp = 1 - Math.Exp(-kappa * smt);

            double result = R * exp / kappa - smt * R - 0.25 * sigma * sigma * exp * exp / (kappa * kappa * kappa);

            return Math.Exp(result);
        }

        private double B(double t, double s)
        {
            return (1 - Math.Exp(-kappa * (s - t))) / kappa;
        }

        public override double P(double t, double s)
        {
            return A(t, s) * Math.Exp(-r0 * B(t, s));
        }

        public override double R(double t, double s)
        {
            return (-Math.Log(A(t, s)) + B(t, s) * r0) / (s - t);
        }

        public override double YieldVolatility(double t, double s)
        {
            return sigma * (1.0 - Math.Exp(-kappa * (s - t))) / (kappa * (s - t));
        }

        public override void Accept(BondVisitor visitor)
        {
            visitor.Visit(this);
        }
    }

    // Class for CIR Model
    public class CIRModel : BondModel
    {

        public CIRModel(double kappa, double theta, double vol, double r) : base(kappa, theta, vol, r)
        {

        }
        public double A(double t, double s)
        {
            double h = Math.Sqrt(kappa * kappa + 2 * sigma * sigma);
            double smt = s - t;
            double power = 2 * kappa * theta / (sigma * sigma);
            double nominator = 2 * h * Math.Exp((kappa + h) * smt / 2);
            double denominator = (Math.Exp(smt * h) - 1) * (kappa + h) + 2 * h;

            return Math.Pow(nominator / denominator, power);
        }

        public double B(double t, double s)
        {
            double h = Math.Sqrt(kappa + 2 * sigma * sigma);
            double smt = s - t;
            double nominator = 2 * (Math.Exp(smt * h) - 1);
            double denominator = (kappa + h) * (Math.Exp(smt * h) - 1) + 2 * h;

            return nominator / denominator;
        }

        public override double P(double t, double s)
        {
            return A(t, s) * Math.Exp(-r0 * B(t, s));
        }

        public override double R(double t, double s)
        {
            return (-Math.Log(A(t, s)) + B(t, s) * r0) / (s - t);
        }

        public override double YieldVolatility(double t, double s)
        {
            return sigma * (1.0 - Math.Exp(-kappa * (s - t))) / (kappa * (s - t));
        }

        public override void Accept(BondVisitor visitor)
        {
            visitor.Visit(this);
        }

    }

    // Visitor Class
    public abstract class BondVisitor
    {
        public BondVisitor()
        {
        }

        public BondVisitor(BondVisitor sourse)
        {
        }

        public abstract void Visit(VasicekModel model);
        public abstract void Visit(CIRModel model);
    }

    public enum OptionType : uint
    {
        Call = 0,
        Put = 1
    }

    public class BondOptionPricing : BondVisitor
    {
        // Base Parameters
        private double t;
        private double T;
        private double s;
        private double K;

        // Type of the options
        private OptionType type;

        // Computed Value
        public double OptionPrice;

        //Constructor
        public BondOptionPricing(double start, double optmat, double bondmat, double strike, OptionType otype)
        {
            t = start;
            T = optmat;
            s = bondmat;
            K = strike;
            type = otype;
        }

        //Methods
        public override void Visit(VasicekModel model)
        {
            if (type == OptionType.Call)
            {
                OptionPrice = CallPrice(model);
            }

            else
            {
                OptionPrice = PutPrice(model);
            }
        }

        public override void Visit(CIRModel model)
        {
            if (type == OptionType.Call)
            {
                OptionPrice = CallPrice(model);
            }

            else
            {
                OptionPrice = PutPrice(model);
            }
        }

        private double CallPrice(VasicekModel model)
        {
            double nu = 0.5 * model.sigma * model.sigma *
                (1.0 - Math.Exp(-2.0 * model.kappa * (T - t))) / model.kappa;

            nu = Math.Sqrt(nu);

            double SigP = nu * (1.0 - Math.Exp(-model.kappa * (s - T))) / model.kappa;

            double d1 = (Math.Log(model.P(t, s) / (model.P(t, T) * K)) / SigP) + 0.5 * SigP;

            double d2 = d1 - SigP;

            return model.P(t, s) * GaussianFunctions.scdf(d1) - K * model.P(t, s) * GaussianFunctions.scdf(d2);
        }

        private double PutPrice(VasicekModel model)
        {
            double nu = 0.5 * model.sigma * model.sigma *
                (1.0 - Math.Exp(-2.0 * model.kappa * (T - t))) / model.kappa;

            nu = Math.Sqrt(nu);

            double SigP = nu * (1.0 - Math.Exp(-model.kappa * (s - T))) / model.kappa;

            double d1 = (Math.Log(model.P(t, s) / (model.P(t, T) * K)) / SigP) + 0.5 * SigP;

            double d2 = d1 - SigP;

            return K * model.P(t, s) * GaussianFunctions.scdf(-d2) - model.P(t, s) * GaussianFunctions.scdf(-d1);
        }

        private double CallPrice(CIRModel model)
        {
            double df = 4.0 * model.kappa * model.theta / (model.sigma * model.sigma);
            double thet = Math.Sqrt(model.kappa * model.kappa + 2 * model.sigma * model.sigma);
            double phi = 2 * thet / (model.sigma * model.sigma * (Math.Exp(thet * (T - t)) - 1));
            double epsi = (model.kappa + thet) / (model.sigma * model.sigma);

            double rStar = Math.Log(model.A(T, s) / K) / model.B(T, s);
            double ncParam = 2 * phi * phi * model.r0 * Math.Exp(thet * (T - t)) / (phi + epsi + model.B(T, s));
            double x = 2 * rStar * (phi + epsi + model.B(T, s));

            double ncParam2 = 2 * phi * phi * model.r0 * Math.Exp(thet * (T - t)) / (phi + epsi);
            double x2 = 2 * rStar * (phi + epsi);

            return model.P(t, s) * NonCentralChiSquared.NCX_cdf(x, ncParam, df) -
                K * model.P(t, s) * NonCentralChiSquared.NCX_cdf(x2, ncParam2, df);
        }

        private double PutPrice(CIRModel model)
        {
            double Call = CallPrice(model);

            return Call - model.P(t, s) + K * model.P(t, T);
        }

    }
}

// Testing
namespace TestRatesModel
{
    using RatesModels;
    using Statistics;

    class Test
    {
        static void Main()
        {

            // Model Parameters
            double kappa = .15;
            double r0 = .05;
            double sigma = .01;
            double theta = r0;

            //Option Parameters
            OptionType type = OptionType.Call;
            double t = 0.0;
            double T = 1.0;
            double s = 5.0;
            double K = 0.67;

            // Model Calibration
            VasicekModel vasicek = new VasicekModel(kappa, theta, sigma, r0);
            CIRModel cir = new CIRModel(kappa, theta, sigma, r0);
            BondOptionPricing bv = new BondOptionPricing(t, T, s, K, type);

            //Vasicek output
            Console.WriteLine("Vasicek ZCB Price is: {0}", vasicek.P(.0, 1));
            Console.WriteLine("Vasicek yield rate is: {0}", vasicek.R(.0, 1));
            Console.WriteLine("Vasicek Yield Volatility is: {0}", vasicek.YieldVolatility(.0, 1));
            bv.Visit(vasicek);
            Console.WriteLine("Vasicek option price is: {0}\n", bv.OptionPrice);

            // CIR Output
            Console.WriteLine("CIR ZCB Price is: {0}", cir.P(.0, 1));
            Console.WriteLine("CIR yield rate is: {0}", cir.R(.0, 1));
            Console.WriteLine("CIR Yield Volatility is: {0}", cir.YieldVolatility(.0, 1));
            bv.Visit(cir);
            Console.WriteLine("CIR option price is: {0}\n", bv.OptionPrice);
        }
    }
}
