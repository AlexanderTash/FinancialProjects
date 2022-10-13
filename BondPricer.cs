using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;

namespace BondPricing
{
    public interface IBondPricer
    {
        public double CleanPrice(double rate);
        public double DirtyPrice(double rate);
        public double AccuredInterest();
        public double Yield(double price);
    }

    public class InterestRateCalculator
    {
        private double r;
        private int periods;

        //Constructors
        public InterestRateCalculator(int NumberPeriods, double Rate)
        {
            r = Rate;
            periods = NumberPeriods;
        }

        public InterestRateCalculator(DateTime start, DateTime end, double Rate)
        {
            r = Rate;
            periods = Convert.ToInt32((end - start).Days) / 365;
        }

        // Methods

        // Future Value Functions
        public double FutureValue(double P0, int m) //P0 - payment, m - number of payments py
        {
            double factor = 1.0 + r / m;
            return P0 * Math.Pow(factor, periods * m);
        }

        // Present Value Functions
        public double PresentValue(double P0, int m)
        {
            double factor = 1.0 + r / m;
            return P0 / Math.Pow(factor, periods * m);
        }

        public double PresentValueFlows(double[] cashflows, int m)
        {
            double factor = 1.0 + r / m;
            double PV = 0;
            for(int i = 0; i < cashflows.Length; i++)
            {
                PV = PV + cashflows[i] * 1 / Math.Pow(factor, i + 1);
            }
            return PV;
        }

        public double PresentValueCoupon(double Coupon, int m)
        {
            double factor = 1.0 / (1.0 + r / m);
            double PV = 0;

            for (int i = 1; i <= m * periods; i++)
            {
                PV = PV + Math.Pow(factor, i);
            }

            return PV * Coupon;
        }

        public double TWPVC(double Coupon, int m)
        {
            double factor = 1.0 / (1.0 + r / m);
            double PV = 0;

            for (int i = 1; i <= m * periods; i++)
            {
                PV = PV + Math.Pow(factor, i) * i;
            }

            return PV * Coupon;
        }

        public int NumberOfPeriods()
        {
            return periods;
        }

        public double Rate()
        {
            return r;
        }
    }

    public class ScheduleGenerator
    {

        public DateTime CheckDate(DateTime Date)
        {
            if (Date.DayOfWeek == DayOfWeek.Saturday)
            {
                Date = Date.AddDays(2);
            }

            if (Date.DayOfWeek == DayOfWeek.Sunday)
            {
                Date = Date.AddDays(1);
            }

            return Date;
        }

        // Fixed Constructor
        public List<DateTime> Schedule(DateTime StartDate, DateTime EndDate, int Frequency)
        {
            List<DateTime> schedule = new List<DateTime>();

            EndDate = CheckDate(EndDate);

            DateTime temp = EndDate;

            while (temp >= CheckDate(StartDate))
            {
                int month = 12 / Frequency;
                temp = temp.AddMonths(-month);
            }

            DateTime FirstDate = temp;
            temp = CheckDate(FirstDate);

            FirstDate = CheckDate(FirstDate);

            if (FirstDate < StartDate)
            {
                FirstDate = CheckDate(temp.AddMonths(12 / Frequency));
            }
            
            schedule.Add(CheckDate(FirstDate).Date);

            while (temp <= CheckDate(EndDate))
            {
                int month = 12 / Frequency;
                temp = temp.AddMonths(month);
                schedule.Add(CheckDate(temp).Date);
            }

            return schedule;
        }

        //Copy constructor
        public List<DateTime> Schedule(List<DateTime> ScheduleInput)
        {
            return ScheduleInput;
        }

    }

    public struct Bond
    {
        //Attributes (Assume act / 365 convention)
        public DateTime Start;
        public DateTime NextCouponDate;
        public DateTime Maturity;

        public double AnnualCouponRate;
        public double FaceValue;
        public int PaymentFrequency;

        public int YearsToMaturity;

        public double YearsToCoupon;


        // Constructors
        public Bond(double coupon, int frequency, DateTime start, DateTime end, double FV = 1000)
        {
            Start = start;

            FaceValue = FV;

            AnnualCouponRate = coupon;
            PaymentFrequency = frequency;

            ScheduleGenerator BondSchedule = new ScheduleGenerator();

            List<DateTime> schedule = BondSchedule.Schedule(start, end, frequency);

            NextCouponDate = schedule[0];
            Maturity = schedule.Last();

            YearsToMaturity = Convert.ToInt32((this.Maturity - this.NextCouponDate).Days) / 365;
            YearsToCoupon = Convert.ToDouble((this.NextCouponDate - this.Start).Days) / 365;
        }
    }


    public class BondPricer : IBondPricer
    {
        Bond mybond;

        // Activate Engine
        public InterestRateCalculator eng;

        // Constructor
        public BondPricer(Bond CustomBond)
        {
            mybond = CustomBond;
        }

        // Methods for Fixed rate constructor
        public double CleanPrice(double rate)
        {
            double YieldToMaturity = rate;
            eng = new InterestRateCalculator(mybond.YearsToMaturity, YieldToMaturity);

            double coupons = mybond.AnnualCouponRate * mybond.FaceValue / mybond.PaymentFrequency;
            double PVCoupons = eng.PresentValueCoupon(Coupon: coupons, m: mybond.PaymentFrequency);
            double PVFaceValue = eng.PresentValue(P0: mybond.FaceValue, m: mybond.PaymentFrequency);

            return PVCoupons + PVFaceValue;
        }

        public double AccuredInterest()
        {
            double accured_coupon = mybond.AnnualCouponRate * mybond.FaceValue / mybond.PaymentFrequency;
            double time = 1 - mybond.YearsToCoupon * 2;
            return accured_coupon * time;
        }

        public double DirtyPrice(double rate)
        {
            return CleanPrice(rate) + AccuredInterest();
        }

        public double Yield(double price)
        {
            double ErrorTolerance = Math.Pow(10, -6);
            double InitialGuess = 0.1;
            double PriceGuessed;

            eng = new InterestRateCalculator(mybond.YearsToMaturity, InitialGuess);

            PriceGuessed = CleanPrice(InitialGuess);

            while (Math.Abs(PriceGuessed - price) > ErrorTolerance)
            {
                if (PriceGuessed - price > 0)
                {
                    InitialGuess = InitialGuess + 0.001;
                    PriceGuessed = CleanPrice(InitialGuess);
                }

                else
                {
                    InitialGuess = InitialGuess - 0.001;
                    PriceGuessed = CleanPrice(InitialGuess);
                }
            }

            return InitialGuess;
        }

        public double MacauleyDuration(double rate)
        {

            eng = new InterestRateCalculator(mybond.YearsToMaturity, rate);

            double price = CleanPrice(rate);
            double coupons = mybond.AnnualCouponRate * mybond.FaceValue / mybond.PaymentFrequency;
            double nom = eng.TWPVC(coupons, mybond.PaymentFrequency);

            double fv = mybond.YearsToMaturity * mybond.PaymentFrequency * eng.PresentValue(mybond.FaceValue, mybond.PaymentFrequency);

            return (nom + fv) / price / mybond.PaymentFrequency;
        }

        public double ModifiedDuration(double rate)
        {
            double macd = MacauleyDuration(rate);

            return macd / (1 + rate);
        }

        public double DV01(double rate)
        {
            double modd = ModifiedDuration(rate);

            return modd * mybond.FaceValue / 10000;
        }

        public double BondConvexity(double rate)
        {
            double actual_diff = CleanPrice(rate) - CleanPrice(rate + 0.0001);
            double estimated_diff = DV01(rate);

            return (actual_diff - estimated_diff);
        }
    }
}



namespace TestBondPricer
{
    using BondPricing;

    class Test
    {
        static void Main(string[] args)
        {
            double r = .092;
            double coup = 0.1;
            int f = 2;

            DateTime st = new DateTime(2000, 1, 1);
            DateTime et = new DateTime(2010, 1, 1);

            Bond mybond1 = new Bond(coup, f, st, et);

            BondPricer bp = new BondPricer(mybond1);

            Console.WriteLine("Bond Clean Price is: {0}", bp.CleanPrice(r));
            Console.WriteLine("Bond Dirty Price is: {0}", bp.DirtyPrice(r));
            Console.WriteLine("Bond YTM is: {0}", bp.Yield(bp.CleanPrice(r)));
            Console.WriteLine("Bond Accured Interest is: {0}", bp.AccuredInterest());
            Console.WriteLine("Bond Macauley Duration is: {0}", bp.MacauleyDuration(r));
            Console.WriteLine("Bond Modified Duration is: {0}", bp.ModifiedDuration(r));
            Console.WriteLine("Bond DV01 is: {0}", bp.DV01(r));
            Console.WriteLine("Bond Convexity is: {0}", bp.BondConvexity(r));
        }
    }
}