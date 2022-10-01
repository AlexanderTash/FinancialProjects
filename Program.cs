using System;

// Main Option pricing engine

namespace OptionPricingBS
{

    //Class for gaussian methods
    public class StandardGaussian
    {
        // PDF of the standard normal distribution
        static public double pdf(double x)
        {
            return Math.Exp(-Math.Pow(x, 2) / 2) / Math.Sqrt(2 * Math.PI);
        }

        //CDF of the standard normal distribution
        static public double cdf(double x)
        {
            double x_min = -5;
            double x_max = 5;
            double interval = 1e5;
            double increase = (x_max - x_min) / interval;

            if(x <= x_min)
            {
                return 0;
            }

            if(x >= x_max)
            {
                return 1;
            }

            else
            {
                double current = x_min;
                double sum = 0;
                while(current < x)
                {
                    sum = sum + StandardGaussian.pdf(current + increase / 2) * increase;
                    current = current + increase;
                }
                return sum;
            }
        }

    }


    // Class for creating options
    public class Option
    {
        public double r;          //Interest rate
        public double sigma;      //Implied Volatility
        public double K;          //Strike Price
        public double T;          //Maturity
        public string type;       //Type: Call or Put

        // Default constructor
        public Option()
        {

        }

        // Custom Constructor
        public Option(double rate, double vol, double strike, double maturity, string otype)
        {
            r = rate;
            sigma = vol;
            K = strike;
            T = maturity;
            type = otype;
        }

        // Private methods

        // Pricees for call and put options
        private double CallPrice(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return S * StandardGaussian.cdf(d1) - K * Math.Exp(-r * T) * StandardGaussian.cdf(d2);
        }

        private double PutPrice(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return K * Math.Exp(-r * T) * StandardGaussian.cdf(-d2) - S * StandardGaussian.cdf(-d1);
        }


        // Delta for call and put options
        private double CallDelta(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return StandardGaussian.cdf(d1);
        }

        private double PutDelta(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return -StandardGaussian.cdf(-d1);
        }


        // Gamma for call and put options
        private double CallGamma(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return StandardGaussian.pdf(d1) / (S * sigma * Math.Sqrt(T));
        }

        private double PutGamma(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return StandardGaussian.pdf(d1) / (S * sigma * Math.Sqrt(T));
        }


        // Vega for call and put options
        private double CallVega(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return StandardGaussian.pdf(d1) * S * Math.Sqrt(T);
        }

        private double PutVega(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            return StandardGaussian.pdf(d1) * S * Math.Sqrt(T);
        }


        // Theta for call and put options
        private double CallTheta(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return -StandardGaussian.pdf(d1) * S * sigma / (2 * Math.Sqrt(T))
                - r * K * Math.Exp(-r * T) * StandardGaussian.cdf(d2);
        }

        private double PutTheta(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return StandardGaussian.pdf(d1) * S * sigma / (2 * Math.Sqrt(T))
                + r * K * Math.Exp(-r * T) * StandardGaussian.cdf(-d2);
        }


        // Rho for call and put options
        private double CallRho(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return T * K * Math.Exp(-r * T) * StandardGaussian.cdf(d2);
        }

        private double PutRho(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return  -K * T * Math.Exp(-r * T) * StandardGaussian.cdf(-d2);
        }


        // Vanna for call and put options
        private double CallVanna(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return StandardGaussian.pdf(d1) * d2 / sigma;
        }

        private double PutVanna(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return StandardGaussian.pdf(d1) * d2 / sigma;
        }


        // Vomma for call and put options
        private double CallVomma(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return CallVega(S) * d1 * d2 / sigma;
        }

        private double PutVomma(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));

            return PutVega(S) * d1 * d2 / sigma;
        }


        // Charm for call and put options
        private double CallCharm(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double dec = ((2 * r) * T - d2 * sigma * Math.Sqrt(T)) / (2 * T * sigma * Math.Sqrt(T));
            return -StandardGaussian.pdf(d1) * dec;
        }

        private double PutCharm(double S)
        {
            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double d2 = (Math.Log(S / K) + (r - Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double dec = ((2 * r) * T - d2 * sigma * Math.Sqrt(T)) / (2 * T * sigma * Math.Sqrt(T));
            return -StandardGaussian.pdf(d1) * dec;
        }


        //Public functions
        //Price 
        public double price(double S)
        {
            if(type == "C")
            {
                return CallPrice(S);
            }
            else
            {
                return PutPrice(S);
            }
        }

        //Delta
        public double delta(double S)
        {
            if (type == "C")
            {
                return CallDelta(S);
            }
            else
            {
                return PutDelta(S);
            }
        }

        //Gamma
        public double gamma(double S)
        {
            if (type == "C")
            {
                return CallGamma(S);
            }
            else
            {
                return PutGamma(S);
            }
        }

        //Vega
        public double vega(double S)
        {
            if (type == "C")
            {
                return CallVega(S);
            }
            else
            {
                return PutVega(S);
            }
        }

        //Theta
        public double theta(double S)
        {
            if (type == "C")
            {
                return CallTheta(S);
            }
            else
            {
                return PutTheta(S);
            }
        }

        //Rho
        public double rho(double S)
        {
            if (type == "C")
            {
                return CallRho(S);
            }
            else
            {
                return PutRho(S);
            }
        }

        //Vanna
        public double vanna(double S)
        {
            if (type == "C")
            {
                return CallVanna(S);
            }
            else
            {
                return PutVanna(S);
            }
        }

        //Vomma
        public double vomma(double S)
        {
            if (type == "C")
            {
                return CallVomma(S);
            }
            else
            {
                return PutVomma(S);
            }
        }

        //Charm
        public double charm(double S)
        {
            if (type == "C")
            {
                return CallCharm(S);
            }
            else
            {
                return PutCharm(S);
            }
        }
    }
}

// Example of Inheritance

namespace OptionExtension
{
    using OptionPricingBS;

    public static class OptionMixins
    {
        public static void Display(this Option option, double S)
        {
            Console.WriteLine("Price: {0}", option.price(S));
            Console.WriteLine("Delta: {0}", option.delta(S));
            Console.WriteLine("Gamma: {0}", option.gamma(S));
            Console.WriteLine("Vega: {0}", option.vega(S));
            Console.WriteLine("Theta: {0}", option.theta(S));
            Console.WriteLine("Rho: {0}", option.rho(S));
            Console.WriteLine("Vanna: {0}", option.vanna(S));
            Console.WriteLine("Vomma: {0}", option.vomma(S));
            Console.WriteLine("Charm: {0}", option.charm(S));
        }

        public static double[] PriceSpots(this Option option, double low, double high, int steps)
        {
            double h = (high - low) / steps;
            double S = low;
            double[] output = new double[steps + 1];
            for(int i = 1; i <= steps; i++)
            {
                output[i] = option.price(S);
                S = S + h;
            }
            return output;
        }
    }
}

//Testing all features

namespace Testing
{
    using OptionPricingBS;
    using OptionExtension;
    using Microsoft.VisualBasic.FileIO;

    public class Program
    {
        static void Main(string[] args)
        {
            // Set option
            Option callOption1 = new Option();
            callOption1.type = "C";
            callOption1.K = 100.0;
            callOption1.T = 1.0;
            callOption1.r = 0.0;
            callOption1.sigma = 0.10;

            Option callOption2 = new Option(0, .1, 100, 1, "C");

            double S = 100.0;

            //Testing main pricing function
            Console.WriteLine("Option 1 price: {0}", callOption1.price(S));
            Console.WriteLine("Option 2 price: {0}", callOption2.price(S));

            //Testing extensions
            callOption1.Display(S);

            double low = 50;
            double up = 150;
            int step = 100;
            double[] prices = callOption1.PriceSpots(low, up, step);

            for (int i = 0; i < prices.Length; i++)
            {
                Console.WriteLine("Option price is: {0}", prices[i]);
            }
        }
    }
}